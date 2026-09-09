# Ordinal (cumulative) MAIHDA.
#
# The frequentist path fits a cumulative link mixed model via ordinal::clmm()
# under a dedicated 'ordinal' engine; the Bayesian path uses brms::cumulative().
# A cumulative model has no usable predict() method (predict.clmm does not
# exist), so -- exactly as the wemix engine does for WeMixResults -- predictions
# are built manually from the stored location coefficients (beta), thresholds
# (alpha) and stratum conditional modes: P(Y <= k) = linkinv(alpha_k - eta) with
# eta = x'beta + u. The MAIHDA variance summaries live on the latent scale, where
# the level-1 variance is the standard pi^2/3 (logit) or 1 (probit) -- the same
# latent treatment the package applies to binomial models -- so
# VPC = sigma^2_u / (sigma^2_u + pi^2/3).
#
# Empirical notes on the clmm object (ordinal 2025.12.29), which the accessors
# below rely on: $alpha (the FREE threshold parameters -- named "1|2", ... and
# equal to the cut points only under the default threshold = "flexible"; a
# structured threshold request stores a shorter reparameterisation such as
# threshold.1/spacing, so the cut points come from $Theta / $tJac via
# maihda_clmm_cutpoints()), $beta (named location
# coefficients, NO intercept -- it is absorbed by the thresholds), $model (the
# model frame), $link, $xlevels, $terms (fixed-effects-only terms), and
# $optRes$convergence (0 = converged). ordinal exports VarCorr(), ranef() and
# condVar() for clmm (condVar returns conditional VARIANCES); nobs() and vcov()
# (which includes the threshold rows) dispatch off the loaded namespace;
# stats::family() is undefined for clmm, so the family the wrapper records at
# fit time is the source of truth downstream.

# Links for which the latent-scale level-1 variance (and hence the VPC) is
# defined; matches the binomial latent treatment elsewhere in the package.
.maihda_ordinal_links <- c("logit", "probit")

# Stop early with an installation hint when the ordinal package is unavailable.
maihda_require_ordinal <- function() {
  if (!requireNamespace("ordinal", quietly = TRUE)) {
    stop("Package 'ordinal' is required for the cumulative (engine = \"ordinal\") ",
         "fit. Please install it with: install.packages('ordinal') -- or use ",
         "engine = \"brms\" for the Bayesian cumulative model.", call. = FALSE)
  }
  invisible(TRUE)
}

#' Cumulative (ordinal) family marker for MAIHDA models
#'
#' @description
#' Specifies a cumulative (proportional-odds) model for an ordinal outcome in
#' \code{\link{fit_maihda}} / \code{\link{maihda}}, with a choice of link:
#' \code{maihda_cumulative("logit")} (the default, equivalent to
#' \code{family = "ordinal"}) or \code{maihda_cumulative("probit")}. It plays the
#' role a \code{stats} family object plays for the other families -- there is no
#' cumulative family constructor in \code{stats}, and using
#' \code{brms::cumulative()} would require brms for a frequentist fit.
#'
#' @param link The cumulative link: \code{"logit"} (default) or \code{"probit"}.
#'   These are the links for which the latent-scale VPC is defined
#'   (level-1 variance \eqn{\pi^2/3} and 1 respectively).
#' @return A family marker list with elements \code{family = "cumulative"} and
#'   \code{link}.
#' @examples
#' maihda_cumulative()
#' maihda_cumulative("probit")
#' @seealso \code{\link{fit_maihda}}
#' @export
maihda_cumulative <- function(link = c("logit", "probit")) {
  link <- match.arg(link)
  list(family = "cumulative", link = link)
}

# TRUE when a (resolved or raw) family specification requests a cumulative
# model: the strings "ordinal"/"cumulative", a marker list from
# maihda_cumulative(), or a brms::cumulative() family object. A bare function
# (e.g. brms::cumulative) is resolved by fit_maihda() before this is consulted.
maihda_family_is_ordinal <- function(family) {
  if (is.character(family) && length(family) == 1) {
    return(maihda_normalize_family_name(family) == "cumulative")
  }
  if (is.list(family) && !is.null(family$family)) {
    return(identical(maihda_normalize_family_name(family$family), "cumulative"))
  }
  FALSE
}

# The MAIHDA summaries for a cumulative model are defined for the logit and
# probit links only (the latent level-1 variance is pi^2/3 / 1); reject other
# cumulative links (cloglog, cauchit, ...) up front.
maihda_ordinal_check_family <- function(family) {
  if (!family$link %in% .maihda_ordinal_links) {
    stop("The cumulative (ordinal) MAIHDA model supports the ",
         paste(.maihda_ordinal_links, collapse = " and "), " links, for which ",
         "the latent-scale VPC is defined; this model uses link = '",
         family$link, "'.", call. = FALSE)
  }
  invisible(TRUE)
}

# The ordinal engine fits the canonical MAIHDA structure only: one
# intercept-only (1 | stratum) random effect. clmm() itself can fit more, but
# the variance/ranef/prediction helpers below (and the MAIHDA VPC) assume the
# single stratum effect -- the same restriction the wemix engine makes.
maihda_ordinal_check_formula <- function(formula) {
  re_terms <- reformulas::findbars(formula)
  ok <- length(re_terms) == 1 &&
    identical(paste(deparse(re_terms[[1]][[2]]), collapse = " "), "1") &&
    identical(all.vars(re_terms[[1]][[3]]), "stratum")
  if (!ok) {
    stop("engine = \"ordinal\" supports the canonical MAIHDA structure only: a ",
         "single intercept-only random effect (1 | stratum) (or the (1 | var1:var2) ",
         "shorthand that resolves to it). For crossed or additional random effects ",
         "(context =, decomposition = \"crossed-dimensions\", extra (1 | g) terms), ",
         "use engine = \"brms\" with family = \"ordinal\".", call. = FALSE)
  }
  invisible(TRUE)
}

# Validate / coerce the response of a cumulative model. clmm() and brms
# cumulative() need an (ordered) factor; a numeric response is rejected with a
# conversion hint rather than silently treated as interval-scaled, and an
# unordered factor is coerced to ordered in its declared level order with a
# message (the order is load-bearing for a cumulative model).
maihda_ordinal_prepare_response <- function(data, formula) {
  resp_expr <- formula[[2]]
  resp_vars <- all.vars(resp_expr)
  if (length(resp_vars) != 1 || !is.symbol(resp_expr)) {
    stop("The cumulative (ordinal) MAIHDA model needs a single outcome column ",
         "as the response (no cbind()/addition terms).", call. = FALSE)
  }
  resp_name <- resp_vars[1]
  if (!resp_name %in% names(data)) {
    stop("Response variable not found in data: ", resp_name, call. = FALSE)
  }
  y <- data[[resp_name]]
  if (!is.factor(y)) {
    stop("The cumulative (ordinal) MAIHDA model needs an ordered-factor ",
         "response; '", resp_name, "' is ", class(y)[1], ". Convert it first, ",
         "e.g. ", resp_name, " = factor(", resp_name, ", levels = ..., ",
         "ordered = TRUE), so the category order is explicit.", call. = FALSE)
  }
  y <- droplevels(y)
  if (nlevels(y) < 3) {
    stop("The cumulative (ordinal) MAIHDA model needs at least 3 response ",
         "categories; '", resp_name, "' has ", nlevels(y), ". A two-level ",
         "outcome is a binomial model (family = \"binomial\").", call. = FALSE)
  }
  if (!is.ordered(y)) {
    message("fit_maihda(): coercing the response '", resp_name, "' to an ",
            "ordered factor in its declared level order (",
            paste(levels(y), collapse = " < "), "). Set the levels explicitly ",
            "if this order is wrong.")
    y <- factor(y, levels = levels(y), ordered = TRUE)
  }
  data[[resp_name]] <- y
  data
}

# Assert a cumulative response still has the >= 3 categories a cumulative model
# needs AFTER analytic-sample filtering (the rows the fit actually uses). This is
# deliberately separate from maihda_ordinal_prepare_response()'s full-data check:
# that check runs before rows with missing predictors/outcomes are dropped, so a
# category present only on excluded rows passes it, yet clmm()/brms fit only the
# observed categories WITHOUT complaint -- silently turning an explicit 3+-level
# ordinal request into a binomial model whose response predictions use the wrong
# 1..K scale. Both the clmm and brms ordinal paths re-run this on exactly the rows
# their fit uses. Expects `y` already droplevels()'d so an empty category counts as
# absent. `resp_name` names the outcome for the message only.
maihda_ordinal_assert_min_levels <- function(y, resp_name) {
  nl <- nlevels(y)
  if (nl < 3) {
    stop("The cumulative (ordinal) MAIHDA model needs at least 3 response ",
         "categories, but only ", nl, " remain in the analytic sample -- the rows ",
         "actually fitted, after dropping rows with a missing outcome or predictor. ",
         "A category of '", resp_name, "' occurs only on excluded rows, so the fit ",
         "would silently collapse to a ", nl, "-category (binomial) model. ",
         "Investigate why that category coincides with missing data, or collapse '",
         resp_name, "' to the categories you can actually model.", call. = FALSE)
  }
  invisible(TRUE)
}

#' Fit a cumulative MAIHDA model via ordinal::clmm
#'
#' Internal engine call for \code{fit_maihda(engine = "ordinal")}. Builds the
#' analytic sample (complete cases on the model variables) so the stored
#' \code{data} matches the rows clmm fits, then calls \code{ordinal::clmm()}
#' with \code{Hess = TRUE} (needed for the threshold standard errors). The
#' analytic frame is passed by NAME (bound in a private environment) so the call
#' clmm stores stays one line: ordinal's \code{print}/\code{summary} methods
#' deparse \code{call$data}, and embedding the frame there made printing an
#' ordinal fit dump the whole data set.
#'
#' @param formula The resolved model formula (with \code{(1 | stratum)}).
#' @param data The data (after strata creation / response preparation).
#' @param family The cumulative family marker (link "logit" or "probit").
#' @param dot_vals Named list of evaluated \code{...} arguments forwarded to
#'   \code{ordinal::clmm()} (e.g. \code{nAGQ}, \code{control}).
#' @return A list with \code{model} (the \code{clmm} fit) and \code{data} (the
#'   analytic data frame actually fitted).
#' @keywords internal
maihda_fit_clmm <- function(formula, data, family, dot_vals) {
  # Keep exactly the rows clmm will fit: the evaluated analytic frame (response and
  # fixed-effect transformations applied, rows missing AFTER those transformations
  # dropped). A raw complete.cases() over the formula's columns misses NAs introduced
  # by a transformed term (e.g. log(x) of x <= 0), which left the stored wrapper data
  # holding more rows than clmm actually fit and broke downstream row alignment. Fall
  # back to the raw check only when the analytic frame cannot be built.
  keep <- maihda_analytic_keep_mask(formula, data)
  if (is.null(keep)) {
    model_vars <- intersect(all.vars(formula), names(data))
    keep <- stats::complete.cases(data[, model_vars, drop = FALSE])
  }
  if (!any(keep)) {
    stop("No usable rows remain for the ordinal fit after dropping rows with ",
         "missing model variables.", call. = FALSE)
  }
  if (sum(!keep) > 0) {
    warning(sprintf(paste0("fit_maihda(): dropped %d row(s) with missing model ",
                           "variables before the ordinal fit."), sum(!keep)),
            call. = FALSE)
    data <- data[keep, , drop = FALSE]
  }

  # Re-validate the response category count on this analytic sample and drop any now-
  # empty category before the fit: a level present only on the rows removed above
  # would otherwise let the explicit ordinal request silently collapse to a binomial
  # clmm fit (see maihda_ordinal_assert_min_levels). droplevels() is a no-op when
  # every declared category is still observed.
  resp_name <- all.vars(formula[[2]])[1]
  data[[resp_name]] <- droplevels(data[[resp_name]])
  maihda_ordinal_assert_min_levels(data[[resp_name]], resp_name)

  # Build the clmm call so it REFERENCES the analytic frame by name instead of
  # embedding it. do.call(ordinal::clmm, list(data = data, ...)) substitutes the
  # values it is handed into the call clmm records, so model$call$data became the
  # whole data frame -- and ordinal's print.clmm/summary.clmm deparse exactly that
  # element ("data:    ..."), turning a routine print() of an ordinal fit into
  # thousands of lines of dumped data (the frame also reached the call head, so
  # deparse(getCall(fit)) echoed clmm's entire source). Binding `data` in a private
  # environment and passing the symbol is the same machinery fit_maihda() uses for
  # the lme4/brms engines, and the reason those fits already print a one-line
  # "Data: data". The formula's environment is pointed at the same env so clmm
  # resolves the symbol however it chooses to evaluate the model frame.
  fit_env <- new.env(parent = environment(formula))
  fit_env$data <- data
  environment(formula) <- fit_env
  args <- list(
    formula = formula,
    data = quote(data),
    link = family$link,
    Hess = TRUE
  )
  fit_call <- as.call(c(list(quote(ordinal::clmm)), args, dot_vals))
  model <- eval(fit_call, fit_env)

  list(model = model, data = data)
}

#' Variance components of a cumulative (clmm) MAIHDA fit
#'
#' Reads the between-stratum variance from \code{ordinal::VarCorr()} and pairs
#' it with the latent-scale level-1 variance (\eqn{\pi^2/3} for logit, 1 for
#' probit), matching the latent treatment of binomial models in the other
#' engines.
#'
#' @param object A \code{maihda_model} with engine \code{"ordinal"}.
#' @return A list with \code{stratum} and \code{residual} variances.
#' @keywords internal
maihda_clmm_variances <- function(object) {
  vc <- tryCatch(ordinal::VarCorr(object$model), error = function(e) NULL)
  if (is.null(vc) || !"stratum" %in% names(vc)) {
    stop("Could not read the 'stratum' random-effect variance from the clmm fit.",
         call. = FALSE)
  }
  var_stratum <- as.numeric(vc[["stratum"]][1, 1])

  var_residual <- if (identical(object$family$link, "probit")) 1 else (pi^2) / 3

  list(stratum = var_stratum, residual = var_residual)
}

#' Expanded cut points of a cumulative (clmm) fit
#'
#' The \eqn{K-1} cut points of a \eqn{K}-category cumulative model, whatever
#' threshold structure was fitted.
#'
#' \code{clmm}'s \code{$alpha} holds the \emph{free} threshold parameters, which
#' equal the cut points only under the default \code{threshold = "flexible"}. A
#' structured threshold request (\code{"equidistant"}, \code{"symmetric"},
#' \code{"symmetric2"}, all part of \code{clmm}'s documented API and reachable
#' through \code{fit_maihda}'s \code{...}) leaves \code{$alpha} holding a shorter,
#' differently-parameterised vector -- \code{threshold.1} and \code{spacing} for
#' the equidistant case -- so reading it as the cut points understates the number
#' of categories and mislocates every one of them. The expanded cut points are
#' stored on the fit as \code{$Theta} (a \eqn{1 \times (K-1)} matrix), equivalently
#' \code{$tJac \%*\% $alpha}; \code{$tJac} is the identity under
#' \code{"flexible"}, so this returns \code{$alpha} unchanged there.
#'
#' @param model A fitted \code{clmm}.
#' @return A numeric vector of \eqn{K-1} increasing cut points, named
#'   \code{"1|2"}, \code{"2|3"}, ... when the fit supplies those labels.
#' @keywords internal
maihda_clmm_cutpoints <- function(model) {
  theta <- model$Theta
  if (!is.null(theta) && length(theta) > 0) {
    # A 1 x (K-1) matrix whose column names are the "1|2", "2|3", ... labels.
    out <- as.numeric(theta)
    nms <- if (is.matrix(theta)) colnames(theta) else names(theta)
  } else {
    # Older/leaner fits may carry only the free parameters and the Jacobian that
    # expands them; fall back to alpha itself when neither is present.
    alpha <- model$alpha
    if (is.null(alpha) || length(alpha) == 0) {
      stop("No thresholds found on the clmm fit.", call. = FALSE)
    }
    tJac <- model$tJac
    if (!is.null(tJac) && is.matrix(tJac) && ncol(tJac) == length(alpha)) {
      out <- as.numeric(tJac %*% as.numeric(alpha))
      nms <- rownames(tJac)
    } else {
      out <- as.numeric(alpha)
      nms <- names(alpha)
    }
  }
  # NA-safe: a label vector carrying NA would make all(nms == "") return NA and
  # turn the guard itself into an error.
  if (!is.null(nms) && length(nms) == length(out) &&
      any(!is.na(nms) & nzchar(nms))) {
    names(out) <- nms
  }

  # Invariant: K - 1 cut points for the K categories actually fitted, so every
  # downstream probability matrix has exactly one column per fitted category. This
  # is the guard that would have caught the free-parameter/cut-point confusion --
  # an equidistant 5-category fit stores 2 free parameters, and reading those as
  # cut points silently produced a 3-column probability matrix.
  y_levels <- model$y.levels
  if (!is.null(y_levels) && length(y_levels) > 0 &&
      length(out) != length(y_levels) - 1L) {
    stop(sprintf(paste0("Cumulative fit inconsistency: %d cut point(s) recovered ",
                        "for %d fitted response categories (expected %d)."),
                 length(out), length(y_levels), length(y_levels) - 1L),
         call. = FALSE)
  }
  out
}

#' Threshold (cut-point) estimates of a cumulative (clmm) MAIHDA fit
#'
#' The thresholds \eqn{\alpha_k} take the place of the intercept in a cumulative
#' model: \eqn{P(Y \le k) = g^{-1}(\alpha_k - \eta)}. Reported as the \eqn{K-1}
#' expanded cut points (\code{\link{maihda_clmm_cutpoints}}), so the table means
#' the same thing under every threshold structure rather than echoing the free
#' parameters of a constrained fit. Standard errors come from the Hessian-based
#' \code{vcov()} (hence \code{Hess = TRUE} at fit time) and degrade to \code{NA}
#' when unavailable; under a structured threshold they are the delta-method SEs of
#' the cut points, \eqn{J V J'} for the fit's threshold Jacobian \code{$tJac}.
#' That Jacobian is the identity under the default \code{threshold = "flexible"},
#' where the transform reproduces the direct \code{vcov()} SEs exactly.
#'
#' @param object A \code{maihda_model} with engine \code{"ordinal"}.
#' @return A data frame with \code{term}, \code{estimate}, \code{se}.
#' @keywords internal
maihda_clmm_thresholds <- function(object) {
  model <- object$model
  theta <- maihda_clmm_cutpoints(model)
  alpha <- model$alpha
  V <- tryCatch(stats::vcov(model), error = function(e) NULL)

  se <- rep(NA_real_, length(theta))
  # length(names(alpha)) == length(alpha) is load-bearing, not belt-and-braces:
  # names(alpha) is character(0) on an unnamed vector, %in% then gives logical(0),
  # and all(logical(0)) is TRUE -- so the row selection below would pass the guard
  # and silently produce a 0 x 0 covariance block. Every clmm fit does name alpha
  # (all four threshold structures), so this only ever fires on a malformed fit;
  # it degrades to NA standard errors instead of a non-conformable error that the
  # caller's tryCatch would swallow, dropping the whole table without a word.
  if (!is.null(V) && !is.null(alpha) && length(alpha) > 0 &&
      length(names(alpha)) == length(alpha) &&
      all(names(alpha) %in% rownames(V))) {
    V_alpha <- V[names(alpha), names(alpha), drop = FALSE]
    tJac <- model$tJac
    if (!is.null(tJac) && is.matrix(tJac) &&
        nrow(tJac) == length(theta) && ncol(tJac) == length(alpha)) {
      # Delta method for the cut points, which are an exact linear function of the
      # free parameters (Theta = tJac %*% alpha). Under "flexible" tJac is the
      # identity and this returns the untransformed vcov() diagonal unchanged.
      se <- sqrt(pmax(diag(tJac %*% V_alpha %*% t(tJac)), 0))
    } else if (identical(as.numeric(alpha), as.numeric(theta))) {
      # No reparameterisation happened (flexible, or a fit carrying only alpha),
      # so the free-parameter SEs ARE the cut-point SEs. Equal LENGTH is not
      # enough: a 3-category equidistant fit has 2 free parameters and 2 cut
      # points, and pairing those SEs with these estimates would be silently
      # wrong. Without a usable Jacobian the SEs stay NA rather than misreport.
      se <- sqrt(pmax(diag(V_alpha), 0))
    }
  }

  labels <- names(theta)
  if (is.null(labels) || length(labels) != length(theta)) {
    labels <- paste0(seq_len(length(theta)), "|", seq_len(length(theta)) + 1L)
  }
  data.frame(
    term = labels,
    estimate = as.numeric(theta),
    se = as.numeric(se),
    row.names = NULL,
    stringsAsFactors = FALSE
  )
}

#' Stratum random-effect table for a cumulative (clmm) fit
#'
#' Mirrors \code{maihda_stratum_ranef_lme4()}: one row per stratum with the
#' conditional mode, its conditional standard error (from
#' \code{ordinal::condVar()}, which returns conditional \emph{variances}), and a
#' 95\% interval. At a boundary fit (zero between-stratum variance) the
#' conditional distribution collapses on 0, so the SE is 0.
#'
#' @param object A \code{maihda_model} with engine \code{"ordinal"}.
#' @return A data frame with \code{stratum}, \code{stratum_id},
#'   \code{random_effect}, \code{se}, \code{lower_95}, \code{upper_95}.
#' @keywords internal
maihda_clmm_stratum_ranef <- function(object) {
  re_list <- tryCatch(ordinal::ranef(object$model), error = function(e) NULL)
  if (is.null(re_list) || !"stratum" %in% names(re_list)) {
    stop("No 'stratum' random effects found in the clmm fit.", call. = FALSE)
  }
  tab <- re_list[["stratum"]]
  cols <- intersect(c("(Intercept)", "Intercept"), colnames(tab))
  if (length(cols) == 0) {
    stop("The 'stratum' random effect must include an intercept for MAIHDA ",
         "stratum estimates.", call. = FALSE)
  }
  re <- stats::setNames(as.numeric(tab[[cols[1]]]), rownames(tab))

  tau2 <- tryCatch(maihda_clmm_variances(object)$stratum,
                   error = function(e) NA_real_)
  if (is.finite(tau2) && tau2 < 1e-8) {
    se <- rep(0, length(re))
  } else {
    cv_list <- tryCatch(ordinal::condVar(object$model), error = function(e) NULL)
    se <- rep(NA_real_, length(re))
    if (!is.null(cv_list) && "stratum" %in% names(cv_list)) {
      cv_tab <- cv_list[["stratum"]]
      cv_col <- intersect(c("(Intercept)", "Intercept"), colnames(cv_tab))
      if (length(cv_col) > 0) {
        cv <- stats::setNames(as.numeric(cv_tab[[cv_col[1]]]), rownames(cv_tab))
        se <- sqrt(pmax(cv[names(re)], 0))
      }
    }
  }

  data.frame(
    stratum = names(re),
    stratum_id = suppressWarnings(as.integer(names(re))),
    random_effect = unname(re),
    se = unname(se),
    lower_95 = unname(re - 1.96 * se),
    upper_95 = unname(re + 1.96 * se),
    stringsAsFactors = FALSE
  )
}

#' Location linear predictor of a cumulative (clmm) fit
#'
#' \code{predict.clmm} does not exist, so the location part
#' \eqn{\eta = x'\beta (+ u)} is built directly: the fixed design matrix is
#' constructed with the training data's factor levels AND transformation basis (so a
#' data-dependent term such as \code{scale(x)} uses the fit's centre and scale rather
#' than recomputing them from \code{newdata}) and multiplied by the
#' location coefficients \code{beta} (a clmm has \emph{no} intercept column --
#' it is absorbed by the thresholds -- so \code{beta}'s names select the right
#' columns), any formula offset term is evaluated on \code{newdata} and added
#' (so an offset-only null model still predicts its offset), and
#' \code{include_re} adds each row's stratum conditional mode (an
#' unseen stratum contributes 0 -- the population-average fallback that
#' \code{\link{predict_maihda}} only reaches when \code{allow_new_levels = TRUE},
#' having otherwise rejected unseen strata upstream). Everything is on the latent
#' (link) scale; map through \code{\link{maihda_ordinal_eta_to_score}} for the
#' response-scale expected category score.
#'
#' @param object A \code{maihda_model} with engine \code{"ordinal"}.
#' @param newdata Data to predict for; defaults to the analytic data.
#' @param include_re Add the stratum random effect (conditional mode)?
#' @return A numeric vector of latent-scale location predictions.
#' @keywords internal
maihda_clmm_linpred <- function(object, newdata = NULL, include_re = TRUE) {
  if (is.null(newdata)) {
    newdata <- object$data
  }
  beta <- object$model$beta

  # Build the fixed-effect model frame once. It supplies both the design matrix
  # (when there are location coefficients) and any formula offset term, so it is
  # needed even for a null (thresholds-only) model that carries only an offset.
  # Terms rebuilt from the FITTED data, so a data-dependent transformation such as
  # scale(x) / poly(x, 2) / ns(x, 3) evaluates on the fit's basis instead of being
  # recomputed from the prediction batch (see maihda_fitted_predict_terms()).
  basis <- maihda_fitted_predict_terms(object$formula, object$data)
  tt <- basis$terms
  mf <- stats::model.frame(tt, newdata, xlev = basis$xlev,
                           na.action = stats::na.pass)

  if (is.null(beta) || length(beta) == 0) {
    # Null (thresholds-only) model: the location fixed part is identically 0
    # (the offset, if any, is added below).
    eta <- rep(0, nrow(newdata))
  } else {
    X <- stats::model.matrix(tt, mf)
    missing_cols <- setdiff(names(beta), colnames(X))
    if (length(missing_cols) > 0) {
      stop("Could not rebuild the clmm design matrix; missing column(s): ",
           paste(missing_cols, collapse = ", "), call. = FALSE)
    }
    eta <- drop(X[, names(beta), drop = FALSE] %*% beta)
  }

  # A formula offset term (offset(.) in the model formula) is part of the latent
  # location clmm fits but is NOT a column of X, so add it explicitly -- including
  # for an offset-only null model, whose location is otherwise identically 0.
  off <- stats::model.offset(mf)
  if (!is.null(off)) {
    eta <- eta + off
  }

  if (include_re) {
    re_tab <- maihda_clmm_stratum_ranef(object)
    re <- stats::setNames(re_tab$random_effect, re_tab$stratum)
    u <- re[as.character(newdata$stratum)]
    u[is.na(u)] <- 0
    eta <- eta + unname(u)
  }
  eta
}

# ---- pure cumulative-probability helpers (shared by the clmm and brms paths) --

#' Category probabilities of a cumulative model
#'
#' Pure function (no fit object): given latent locations \code{eta}, ordered
#' thresholds \code{alpha} and the link, returns the category probability matrix
#' via \eqn{P(Y \le k) = g^{-1}(\alpha_k - \eta)} and differencing. Rows are
#' observations, columns categories \code{1..K} (\code{K = length(alpha) + 1}).
#'
#' @param eta Numeric vector of latent locations.
#' @param thresholds Numeric vector of increasing thresholds \eqn{\alpha_k}.
#' @param link \code{"logit"} or \code{"probit"}.
#' @return A numeric matrix with \code{length(eta)} rows that sum to 1.
#' @keywords internal
maihda_ordinal_category_probs <- function(eta, thresholds, link = "logit") {
  if (!link %in% .maihda_ordinal_links) {
    stop("Unsupported cumulative link: ", link, call. = FALSE)
  }
  linkinv <- if (identical(link, "probit")) stats::pnorm else stats::plogis
  thresholds <- as.numeric(thresholds)
  # Strictly increasing, as the @param doc states ("increasing thresholds"). Equal
  # adjacent thresholds imply a zero-probability category (a degenerate cut point),
  # so reject them rather than silently returning an empty category.
  if (is.unsorted(thresholds, strictly = TRUE)) {
    stop("Cumulative thresholds must be strictly increasing.", call. = FALSE)
  }
  cum <- vapply(thresholds, function(a) linkinv(a - eta),
                numeric(length(eta)))
  cum <- matrix(cum, nrow = length(eta))
  full <- cbind(cum, 1)
  probs <- full - cbind(0, cum)
  # Numerical guard: differencing can leave tiny negatives.
  probs[probs < 0] <- 0
  colnames(probs) <- as.character(seq_len(ncol(probs)))
  probs
}

#' Expected category score from a probability matrix
#'
#' The response-scale summary of a cumulative model used throughout the package
#' (the plot layer's "Average Expected Category Score"): \eqn{\sum_k k\, p_k},
#' with categories scored 1..K in order.
#'
#' @param probs A category-probability matrix (rows = observations).
#' @return A numeric vector of expected scores in \eqn{[1, K]}.
#' @keywords internal
maihda_ordinal_expected_score <- function(probs) {
  drop(probs %*% seq_len(ncol(probs)))
}

#' Latent location to expected category score
#'
#' Convenience composition of \code{\link{maihda_ordinal_category_probs}} and
#' \code{\link{maihda_ordinal_expected_score}}.
#'
#' @param eta Numeric vector of latent locations.
#' @param thresholds Numeric vector of increasing thresholds.
#' @param link \code{"logit"} or \code{"probit"}.
#' @return A numeric vector of expected category scores.
#' @keywords internal
maihda_ordinal_eta_to_score <- function(eta, thresholds, link = "logit") {
  maihda_ordinal_expected_score(
    maihda_ordinal_category_probs(eta, thresholds, link)
  )
}

#' Cumulative brmsfamily for the brms ordinal engine
#'
#' Turns whatever the caller supplied for an ordinal \code{engine = "brms"} fit
#' into a proper \code{brms::cumulative()} family, preserving the
#' cumulative-specific options a user set.
#'
#' A rebuild is needed because \code{family} may still be the plain
#' \code{list(family = "cumulative", link = )} marker that the family-string path
#' and \code{\link{maihda_cumulative}} produce, which \code{brm()} cannot use. It
#' must be a rebuild rather than a pass-through for the same reason. Rebuilding
#' from the link ALONE, however, silently reset every other option to brms's
#' defaults: \code{brms::cumulative(threshold = "equidistant")} fitted the
#' flexible model, with no warning that a different model had been fitted from
#' the one asked for. \code{threshold} and \code{link_disc} are therefore carried
#' across when the supplied family actually has them; the marker lists carry
#' neither, so those paths still get brms's defaults.
#'
#' Note the contrast with the \code{clmm} engine, where a structured threshold
#' also has to be handled on the way OUT
#' (\code{\link{maihda_clmm_cutpoints}}). brms needs no such treatment: its
#' generated quantities expand the constrained thresholds into the full
#' \code{b_Intercept[1..K-1]} vector, so \code{brms::fixef()} reports all
#' \eqn{K-1} cut points under every threshold structure and
#' \code{\link{maihda_brms_ordinal_thresholds}} reads them unchanged.
#'
#' @param family The family carried into the brms branch of \code{fit_maihda()}:
#'   a \code{brmsfamily}, or the plain cumulative marker list.
#' @return A \code{brms::cumulative()} family object.
#' @keywords internal
maihda_brms_cumulative_family <- function(family) {
  args <- list(link = family$link)
  for (nm in c("threshold", "link_disc")) {
    val <- family[[nm]]
    # Absent or explicitly unset: leave it to brms's default.
    if (is.null(val) || length(val) == 0L) next
    if (length(val) == 1L && (is.na(val) || !nzchar(as.character(val)))) next
    # Anything else is a setting the caller made, so hand it to brms and let brms
    # judge it. A malformed value must ERROR rather than quietly fall back to the
    # default -- silently substituting a different model is the defect this
    # function exists to fix, and that argument applies to bad input too.
    args[[nm]] <- as.character(val)
  }
  do.call(brms::cumulative, args)
}

#' Posterior-mean cumulative thresholds of a brms cumulative fit
#'
#' The brms analogue of \code{clmm}'s \code{object$model$alpha}: the ordered cut
#' points \eqn{\alpha_k} of a \code{brms::cumulative()} fit, read as posterior
#' means from \code{brms::fixef()}, where they appear as \code{Intercept[1]},
#' \code{Intercept[2]}, ... (any location predictors are separate rows and are
#' dropped here). Returned in threshold order so they pair with the
#' \code{brms::posterior_linpred(re_formula = NA)} location -- which excludes the
#' thresholds, exactly the latent \eqn{\eta} that
#' \code{\link{maihda_ordinal_category_probs}} expects (\eqn{P(Y \le k) =
#' g^{-1}(\alpha_k - \eta)}).
#'
#' Unlike \code{clmm}'s \code{$alpha}, this needs no expansion under a structured
#' \code{threshold}. brms keeps the constrained parameterisation in the model
#' block only -- an equidistant fit samples \code{first_Intercept} and
#' \code{delta} -- and its generated quantities write the full
#' \code{b_Intercept[1..K-1]} vector, which is what \code{fixef()} reports. Verified
#' on a fitted equidistant model: \code{fixef()} carries all \eqn{K-1}
#' \code{Intercept[k]} rows (the raw \code{delta} is not among them) and the cut
#' points come back equally spaced. Contrast
#' \code{\link{maihda_clmm_cutpoints}}, which must do that expansion itself.
#'
#' @param model A fitted \code{brmsfit} from \code{brms::cumulative()}.
#' @return A numeric vector of thresholds (length \eqn{K-1} for \eqn{K} categories).
#' @keywords internal
maihda_brms_ordinal_thresholds <- function(model) {
  fx <- tryCatch(brms::fixef(model), error = function(e) NULL)
  if (is.null(fx) || is.null(dim(fx)) || is.null(rownames(fx)) ||
      !"Estimate" %in% colnames(fx)) {
    stop("Could not read the cumulative thresholds from the brms fit ",
         "(brms::fixef() returned no usable population-level summary).",
         call. = FALSE)
  }
  rn <- rownames(fx)
  thr_rows <- grep("^Intercept\\[[0-9]+\\]$", rn)
  if (length(thr_rows) == 0) {
    stop("Could not find cumulative thresholds (Intercept[k]) in the brms fit; ",
         "is this a brms::cumulative() model?", call. = FALSE)
  }
  ord <- order(as.integer(sub("^Intercept\\[([0-9]+)\\]$", "\\1", rn[thr_rows])))
  as.numeric(fx[thr_rows[ord], "Estimate"])
}

#' Per-stratum predictions for a cumulative (clmm) fit
#'
#' Ordinal counterpart of \code{maihda_stratum_predictions_wemix()}: per-stratum
#' aggregates of the location prediction plus the stratum effect, on the latent
#' (link) scale or as the expected category score (response scale).
#'
#' @param object A \code{maihda_model} with engine \code{"ordinal"}.
#' @param summary_obj Its \code{maihda_summary} (for the stratum estimates).
#' @param scale "response" (expected category score) or "link" (latent).
#' @return A data frame as from \code{maihda_weighted_stratum_aggregate()}.
#' @keywords internal
maihda_stratum_predictions_ordinal <- function(object, summary_obj,
                                               scale = c("response", "link")) {
  scale <- match.arg(scale)
  data <- object$data
  if (!"stratum" %in% names(data)) {
    stop("'stratum' variable not found in fitted model data.")
  }

  cutpoints <- maihda_clmm_cutpoints(object$model)
  link <- object$family$link
  prior_w <- maihda_prediction_weights(object)
  eta_fixed <- maihda_clmm_linpred(object, include_re = FALSE)

  stratum_est <- summary_obj$stratum_estimates
  if (is.null(stratum_est) || nrow(stratum_est) == 0) {
    stop("No stratum estimates available.")
  }

  key <- as.character(data$stratum)
  idx <- match(key, as.character(stratum_est$stratum))
  transform_eta <- function(eta) {
    if (scale == "response") maihda_ordinal_eta_to_score(eta, cutpoints, link) else eta
  }

  pred_df <- data.frame(
    stratum = key,
    predicted_row = transform_eta(eta_fixed + stratum_est$random_effect[idx]),
    lower_row = transform_eta(eta_fixed + stratum_est$lower_95[idx]),
    upper_row = transform_eta(eta_fixed + stratum_est$upper_95[idx]),
    fixed_row = transform_eta(eta_fixed),
    weight = prior_w,
    stringsAsFactors = FALSE
  )

  maihda_weighted_stratum_aggregate(
    pred_df, c("predicted_row", "lower_row", "upper_row", "fixed_row")
  )
}

# ---- calibrated proportional-odds test ---------------------------------------

#' Parametric-bootstrap proportional-odds test for a cumulative MAIHDA fit
#'
#' Tests the proportional-odds (parallel-lines) assumption of a fitted
#' cumulative \code{clmm} MAIHDA model by calibrating the omnibus
#' nominal-effects likelihood-ratio statistic against its own null distribution
#' under the fitted model.
#'
#' @details
#' The statistic is the ordinary omnibus nominal-effects LRT: the fixed-effect
#' part of the model is refitted twice with \code{ordinal::clm()} -- once with
#' all covariates proportional, once with every covariate entering as a
#' threshold-specific (nominal) effect -- and twice the log-likelihood
#' difference is taken. Because \code{clm()} has no random effect, that statistic
#' is computed on the \emph{marginal} fit.
#'
#' Referring it to a chi-squared distribution, as an ordinary
#' \code{ordinal::nominal_test()} would, is not valid here. A conditional
#' cumulative model with a normal random intercept does not in general remain an
#' ordinary proportional-odds model after the random intercept is marginalised
#' away: the implied marginal cumulative-logit slopes differ across thresholds
#' for any non-zero stratum variance. The chi-squared null is therefore false
#' under the correctly specified model, and its rejection rate grows without
#' bound in the sample size -- at a stratum VPC near 7 percent, a plausible
#' MAIHDA value, a correctly specified model is rejected about a quarter of the
#' time at \code{n = 96,000}. Stratum heterogeneity and genuine non-proportional
#' odds are not separable by the fixed-only statistic alone.
#'
#' The one exception is an exactly symmetric threshold configuration, where the
#' marginal slopes coincide and the fixed-only statistic is valid: for symmetric
#' \eqn{u} the map \eqn{\eta \mapsto} logit \eqn{E[\mathrm{plogis}(\eta - u)]}
#' is odd, so its derivative is even and thresholds placed symmetrically about the
#' location share a slope. Three categories cut at \eqn{-c} and \eqn{+c} is the
#' case that arises in practice. Asymmetric thresholds are markedly worse rather
#' than better, so this exception narrows the problem without softening it.
#'
#' This function removes that confounding by simulating the null distribution
#' \emph{under the fitted \code{clmm} itself}: each replicate redraws the stratum
#' random effects from \eqn{N(0, \tau^2)} at the fitted variance, forms the
#' conditional category probabilities from the fitted thresholds and location
#' predictor, redraws the ordinal response, and recomputes the same fixed-only
#' statistic. The reported p-value is
#' \eqn{(1 + \#\{T_b \ge T_{obs}\}) / (1 + B)}, so proportional-odds data
#' generated from the fitted model rejects at the nominal rate by construction.
#'
#' \code{ordinal} supplies no \code{simulate()} method for \code{clmm}, so the
#' simulation is built directly from the fitted thresholds, location
#' coefficients and random-effect variance.
#'
#' The test is opt-in because it is expensive: every replicate refits two
#' \code{clm()} models, so the cost is roughly \code{n_sim} times the cost of the
#' fixed-only refit and grows with the sample size. It is not run automatically
#' at fit time.
#'
#' @param object A \code{maihda_model} fitted with \code{engine = "ordinal"}.
#' @param n_sim Number of parametric-bootstrap replicates (default 199).
#' @param seed Optional integer seed, for a reproducible bootstrap.
#' @return An object of class \code{maihda_po_test}: a list with \code{lrt},
#'   \code{df}, \code{n_terms}, \code{p_value} (the bootstrap p-value),
#'   \code{p_chisq}, \code{n_sim} (replicates that produced a usable statistic),
#'   \code{n_failed}, and \code{null_lrt} (the simulated null statistics).
#'   \code{p_value} is the only p-value the print method shows. \code{p_chisq}
#'   is the uncalibrated chi-squared p-value that the removed automatic screen
#'   used; it is retained on the object for comparison but deliberately not
#'   printed, and it is not evidence against the fitted model.
#' @seealso \code{\link{fit_maihda}}, \code{\link{maihda_cumulative}}
#' @examples
#' \donttest{
#' strata <- make_strata(maihda_sim_data, vars = c("gender", "race"))
#' d <- strata$data
#' d$y <- factor(cut(d$health_outcome, 3), labels = 1:3, ordered = TRUE)
#' m <- fit_maihda(y ~ age + (1 | stratum), data = d, family = "ordinal")
#' maihda_proportional_odds_test(m, n_sim = 99, seed = 1)
#' }
#' @export
maihda_proportional_odds_test <- function(object, n_sim = 199, seed = NULL) {
  if (!inherits(object, "maihda_model") || !inherits(object$model, "clmm")) {
    stop("maihda_proportional_odds_test() needs a maihda_model fitted with ",
         "engine = 'ordinal' (a cumulative clmm fit).", call. = FALSE)
  }
  if (!requireNamespace("ordinal", quietly = TRUE)) {
    stop("Package 'ordinal' is required for maihda_proportional_odds_test().",
         call. = FALSE)
  }
  if (length(n_sim) != 1L || !is.numeric(n_sim) || !is.finite(n_sim) ||
      n_sim < 1 || n_sim != round(n_sim)) {
    stop("'n_sim' must be a single whole number of at least 1.", call. = FALSE)
  }
  n_sim <- as.integer(n_sim)
  if (!is.null(seed)) {
    set.seed(seed)
  }

  obs <- maihda_ordinal_po_stat(object$model)
  if (is.null(obs)) {
    stop("The proportional-odds statistic could not be computed for this fit ",
         "(a null, covariate-free model has no covariate slopes to test).",
         call. = FALSE)
  }

  # Ingredients of the fitted conditional model. maihda_ordinal_check_formula()
  # guarantees the clmm path is the canonical single (1 | stratum) structure, so
  # one random-effect SD and one grouping factor are all that is needed.
  cutpoints <- maihda_clmm_cutpoints(object$model)
  link <- object$family$link
  tau <- sqrt(max(maihda_clmm_variances(object)$stratum, 0))
  eta_fixed <- maihda_clmm_linpred(object, include_re = FALSE)

  # The bootstrap statistic must be recomputed on the SAME frame the observed one
  # used, with only the response redrawn, so reuse that frame verbatim.
  f <- stats::formula(object$model)
  rhs_terms <- attr(stats::terms(f), "term.labels")
  fixed_terms <- rhs_terms[!grepl("\\|", rhs_terms)]
  resp <- all.vars(f)[1]
  dat <- tryCatch(stats::model.frame(object$model),
                  error = function(e) object$model$model)

  grp <- factor(as.character(object$data$stratum))
  if (is.null(dat) || nrow(dat) != length(eta_fixed) ||
      length(grp) != length(eta_fixed)) {
    stop("Could not align the clmm model frame with the analytic data; the ",
         "proportional-odds bootstrap cannot be run on this fit.", call. = FALSE)
  }
  gi <- as.integer(grp)
  n_grp <- nlevels(grp)
  K <- length(cutpoints) + 1L
  y_levels <- levels(dat[[resp]])
  if (length(y_levels) != K) {
    y_levels <- as.character(seq_len(K))
  }

  null_lrt <- rep(NA_real_, n_sim)
  for (b in seq_len(n_sim)) {
    eta <- eta_fixed + stats::rnorm(n_grp, 0, tau)[gi]
    probs <- maihda_ordinal_category_probs(eta, cutpoints, link = link)
    # Inverse-CDF draw from each row's category distribution: the first category
    # whose cumulative probability reaches the row's uniform draw. The final
    # cumulative column is pinned to 1 because a floating-point cumsum can land a
    # hair below it, which would leave a row all-FALSE and silently return
    # category 1 -- the opposite end of the scale -- instead of category K.
    cum <- t(apply(probs, 1, cumsum))
    cum[, ncol(cum)] <- 1
    yb <- max.col(stats::runif(nrow(cum)) <= cum, ties.method = "first")
    dat[[resp]] <- factor(y_levels[yb], levels = y_levels, ordered = TRUE)
    stat <- maihda_po_lrt(dat, resp, fixed_terms)
    if (!is.null(stat)) {
      null_lrt[b] <- stat$lrt
    }
  }

  ok <- null_lrt[is.finite(null_lrt)]
  if (length(ok) == 0) {
    stop("Every parametric-bootstrap replicate failed to refit; the ",
         "proportional-odds test cannot be calibrated for this model.",
         call. = FALSE)
  }

  structure(
    list(lrt = obs$lrt, df = obs$df, n_terms = obs$n_terms,
         p_value = (1 + sum(ok >= obs$lrt)) / (1 + length(ok)),
         p_chisq = stats::pchisq(obs$lrt, df = obs$df, lower.tail = FALSE),
         n_sim = length(ok), n_failed = n_sim - length(ok), null_lrt = ok),
    class = "maihda_po_test"
  )
}

#' @export
print.maihda_po_test <- function(x, ...) {
  pal <- maihda_palette()
  cat(pal$bold("Proportional-odds test (parametric bootstrap under the fitted clmm)"),
      "\n\n", sep = "")
  cat(sprintf("  Nominal-effects LRT : %.3f on %d df over %d covariate(s)\n",
              x$lrt, as.integer(x$df), as.integer(x$n_terms)))
  cat(sprintf("  Bootstrap p-value   : %.4f  (%d replicates)\n",
              x$p_value, as.integer(x$n_sim)))
  if (x$n_failed > 0) {
    cat(sprintf("  Replicates dropped  : %d (refit failed)\n",
                as.integer(x$n_failed)))
  }
  # The uncalibrated chi-squared p-value ($p_chisq) is deliberately NOT printed.
  # It is an inference against a null that is false under the fitted conditional
  # model (see the details of maihda_proportional_odds_test), so showing it beside
  # the bootstrap p-value would give the reader two answers to one question and
  # invite them to choose. It stays on the returned object for comparison.
  invisible(x)
}
