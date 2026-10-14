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
# $optRes (the return of the optimiser that made the fit, nlminb() or ucminf(),
# read by maihda_clmm_converged()). ordinal exports VarCorr(), ranef() and
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

# ---- clmm formula normalization: offset() must precede the random-effect terms ----

# Split a formula right-hand side into its top-level `+` operands, in order.
# Only `+` is split: a `-` subtree (y ~ a - b) stays one operand, so nothing is
# silently promoted out of a subtraction.
maihda_rhs_plus_terms <- function(expr) {
  if (is.call(expr) && identical(expr[[1L]], as.name("+")) && length(expr) == 3L) {
    return(c(maihda_rhs_plus_terms(expr[[2L]]), list(expr[[3L]])))
  }
  list(expr)
}

# Rebuild a right-hand side from an ordered list of `+` operands.
maihda_rhs_from_plus_terms <- function(terms_list) {
  Reduce(function(a, b) call("+", a, b), terms_list)
}

# TRUE for a random-effect operand -- (1 | g), (1 || g), or any parenthesized
# expression whose top-level operator is a bar.
maihda_is_bar_term <- function(expr) {
  while (is.call(expr) && identical(expr[[1L]], as.name("(")) && length(expr) == 2L) {
    expr <- expr[[2L]]
  }
  is.call(expr) &&
    (identical(expr[[1L]], as.name("|")) || identical(expr[[1L]], as.name("||")))
}

# TRUE for an operand that CARRIES a formula offset: offset(<expr>) itself, a
# redundantly parenthesized one, or a `+` group containing one. stats::terms()
# finds an offset through `(` and `+` (and through `-`, `*`, `:`), but NOT through
# an ordinary call -- I(offset(x)) and log(offset(x)) are plain variables, not
# offsets -- so this mirrors the `(`/`+` part of that traversal, which is every
# spelling the package or a user realistically writes. A `-`, `*` or `:` subtree is
# deliberately NOT rewritten (moving it would change the fixed design); the
# alignment check below catches such a formula and errors instead of mis-fitting.
maihda_is_offset_term <- function(expr) {
  while (is.call(expr) && identical(expr[[1L]], as.name("(")) && length(expr) == 2L) {
    expr <- expr[[2L]]
  }
  if (!is.call(expr)) {
    return(FALSE)
  }
  if (identical(expr[[1L]], as.name("offset"))) {
    return(TRUE)
  }
  if (identical(expr[[1L]], as.name("+")) && length(expr) == 3L) {
    return(any(vapply(maihda_rhs_plus_terms(expr), maihda_is_offset_term,
                      logical(1))))
  }
  FALSE
}

# Verify the invariant ordinal::clmm() depends on: every offset position read off
# the BAR-FREE terms must point at the SAME expression in the variables list of the
# BARRED (subbars) formula, because clmm takes the index from the former and the
# frame from the latter. maihda_offset_before_bars() establishes this for every
# spelling it can rewrite; this is the backstop for the ones it cannot (an offset
# buried in a `-`, `*` or `:` subtree after the random effect), so those ERROR with
# an actionable message instead of being fitted against the wrong column. Returns
# TRUE when the question cannot be decided (a formula terms() cannot process):
# conservative, never block a fit that is not provably broken.
maihda_clmm_offset_aligned <- function(formula) {
  ok <- tryCatch({
    tf <- stats::terms(maihda_nobars(formula))
    oi <- attr(tf, "offset")
    if (is.null(oi) || length(oi) == 0L) {
      TRUE
    } else {
      vf <- as.list(attr(tf, "variables"))[-1L]
      vfull <- as.list(attr(stats::terms(reformulas::subbars(formula)),
                            "variables"))[-1L]
      max(oi) <= length(vfull) &&
        all(vapply(oi, function(i) identical(vf[[i]], vfull[[i]]), logical(1)))
    }
  }, error = function(e) TRUE)
  isTRUE(ok)
}

# Move every formula offset() term ahead of the random-effect terms.
#
# ordinal::clmm() builds its model frame from subbars(formula) -- the bars turned
# into `+` -- but then OVERWRITES that frame's terms attribute with the terms of
# nobars(formula), and reads the offset with model.offset(), i.e. by POSITION in
# the bar-free variables list. The two variable lists only line up while every
# random-effect term follows every offset term. With the offset written last --
# which is exactly what stats::update(fixed ~ ., . ~ . + (1 | stratum)) produces,
# and what a user writing y ~ x + (1 | stratum) + offset(off) supplies -- the
# index lands on the grouping variable instead, so clmm silently fits the STRATUM
# COLUMN as the offset (integer stratum ids => wrong thresholds, wrong slopes and
# a between-stratum variance inflated to var(stratum id)); when the mis-indexed
# column is a factor or character it errors instead ("'offset' must be numeric",
# "non-numeric argument to binary operator"). Relocating the offset restores the
# prefix alignment: no bar term then precedes it, so its position in the bar-free
# variables list is the position it occupies in the full one.
#
# Nothing else moves -- the relative order of every other term, the response, the
# intercept and the formula environment are preserved -- and a formula with no
# offset, no random effect, or its offsets already ahead of the bars is returned
# UNCHANGED, so a no-offset fit is bit-identical. clmm-only: lme4 and brms read the
# offset from their own frame and are indifferent to the term order, and the wemix
# engine refuses offsets, which WeMix::mix() does not fit.
#
# A trailing subtraction is walked through rather than given up on. `- 1` is what
# stats::update() leaves behind for every no-intercept formula, so a user writing
# y ~ 0 + x + a + offset(off) + (1 | g) reaches the engine as
# y ~ x + a + (1 | g) + offset(off) - 1, whose OUTERMOST operator is `-`: the
# plus-splitter sees one operand, no relocation happens, and the backstop below
# refuses a fit that is perfectly relocatable. Only the MINUEND is reordered and
# only when no subtracted operand carries an offset -- an offset inside the
# subtraction (y ~ x + (1 | g) - offset(off)) must not be lifted out, because that
# would change the fixed design, and it stays for the backstop to refuse.
maihda_offset_before_bars_rhs <- function(expr) {
  if (is.call(expr) && identical(expr[[1L]], as.name("-")) && length(expr) == 3L) {
    if (maihda_is_offset_term(expr[[3L]])) {
      return(expr)
    }
    return(call("-", maihda_offset_before_bars_rhs(expr[[2L]]), expr[[3L]]))
  }
  parts <- maihda_rhs_plus_terms(expr)
  if (length(parts) < 2L) {
    return(expr)
  }
  is_off <- vapply(parts, maihda_is_offset_term, logical(1))
  is_bar <- vapply(parts, maihda_is_bar_term, logical(1))
  if (!any(is_off) || !any(is_bar)) {
    return(expr)
  }
  first_bar <- which(is_bar)[1L]
  if (!any(is_off & seq_along(parts) > first_bar)) {
    return(expr)
  }
  # Stable relocation: the offsets keep their relative order and land immediately
  # before the first random-effect term; everything else keeps its order too.
  idx <- seq_along(parts)
  reordered <- c(idx[!is_off & idx < first_bar], idx[is_off],
                 idx[!is_off & idx >= first_bar])
  maihda_rhs_from_plus_terms(parts[reordered])
}

maihda_offset_before_bars <- function(formula) {
  if (!inherits(formula, "formula") || length(formula) < 2L) {
    return(formula)
  }
  rhs <- formula[[length(formula)]]
  new_rhs <- maihda_offset_before_bars_rhs(rhs)
  # identical() on the language objects, so a formula nothing moved in is returned
  # as it came -- the bit-identical guarantee the no-offset and already-ordered
  # paths rely on.
  if (identical(new_rhs, rhs)) {
    return(formula)
  }
  out <- formula
  out[[length(out)]] <- new_rhs
  environment(out) <- environment(formula)
  out
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
#' ordinal fit dump the whole data set. Any formula \code{offset()} term is moved
#' ahead of the \code{(1 | stratum)} term first, because clmm reads the offset by
#' position in the bar-free variables list (see \code{maihda_offset_before_bars}).
#'
#' A Laplace fit that ends with every stratum effect at exactly zero, or with a
#' between-stratum SD below 1\% of the latent residual SD, is compared with the
#' boundary fit: \code{clmm()}'s general routine
#' (\code{clmm.control(useMatrix = TRUE)}) started from the first fit's thresholds
#' and coefficients and an SD of zero. The boundary fit is the one returned unless
#' its log-likelihood is below the first fit's by more than 1e-9 of its size.
#'
#' The model is fitted on its fixed-effects design with every column centred and
#' scaled, the boundary fit included, and returned as the \code{clmm()} fit of the
#' formula at that maximum. A covariate far from zero or on a large or tiny scale
#' can leave \code{clmm()}'s optimiser at its starting values or short of the maximum,
#' and \code{clmm()} takes the gradient and Hessian by differences with a fixed step,
#' too coarse for the coefficient of such a covariate. The returned fit carries the
#' gradient and Hessian of the rescaled fit, mapped exactly, so its standard errors
#' are those of that fit. Under \code{threshold = "symmetric2"} the thresholds cannot
#' take up a shift, and the columns are scaled only. A fit given a \code{start} is
#' made as the formula stands, since the start is in the formula's own terms.
#'
#' @param formula The resolved model formula (with \code{(1 | stratum)}).
#' @param data The data (after strata creation / response preparation).
#' @param family The cumulative family marker (link "logit" or "probit").
#' @param dot_vals Named list of evaluated \code{...} arguments forwarded to
#'   \code{ordinal::clmm()} (e.g. \code{nAGQ}, \code{control}).
#' @return A list with \code{model} (the \code{clmm} fit), \code{data} (the
#'   analytic data frame actually fitted) and \code{coding} (the fit's factor levels
#'   and contrast matrices, taken while the options that shaped the design are in
#'   force; \code{NULL} if they could not be read).
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

  # Put any formula offset() term ahead of the (1 | stratum) term before handing
  # the formula to clmm: clmm reads the offset by POSITION in the bar-free
  # variables list while indexing the frame it built from the barred one, so a
  # trailing offset makes it fit the grouping column as the offset (see
  # maihda_offset_before_bars). A no-offset formula is returned unchanged.
  formula <- maihda_offset_before_bars(formula)
  # Backstop: a spelling the relocation could not rewrite (an offset buried in a
  # `-`, `*` or `:` subtree after the random effect) would still be read off the
  # wrong column, silently. Refuse it with an actionable message instead.
  if (!maihda_clmm_offset_aligned(formula)) {
    stop("The offset() term in this formula cannot be moved ahead of the ",
         "random-effect term automatically, and ordinal::clmm() would read the ",
         "wrong column as the offset. Write the offset as its own term BEFORE ",
         "the random effect, e.g. y ~ x + offset(log(n)) + (1 | stratum).",
         call. = FALSE)
  }

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
  # Forward the extra arguments under the formals clmm() binds them to, so the coding
  # record below reads the contrasts argument clmm() used even when it was spelled
  # partially (see maihda_resolve_dot_names()).
  dot_vals <- maihda_resolve_dot_names(dot_vals, ordinal::clmm, names(args))
  fit_call <- as.call(c(list(quote(ordinal::clmm)), args, dot_vals))
  # The fit is made on the fixed-effects design with its columns centred and scaled
  # and returned in the formula's own terms (see maihda_clmm_rescaled_fit()); a fit
  # that cannot be made that way is made as the formula stands.
  model <- maihda_clmm_rescaled_fit(args, dot_vals, fit_env, fit_call)
  if (is.null(model)) {
    model <- maihda_clmm_fit_and_refit(fit_call, args, dot_vals, fit_env)
  }

  # Record the factor coding of the fitted design NOW, while the options and the
  # contrasts argument that shaped it are still in force: predictions rebuild the
  # design by hand and must reuse exactly this coding (see
  # maihda_engine_fixed_coding()). A failure here must not break a fit clmm accepted.
  coding <- tryCatch(
    maihda_engine_fixed_coding(model, maihda_nobars(formula), data,
                               contrasts = dot_vals[["contrasts"]]),
    error = function(e) NULL)

  list(model = model, data = data, coding = coding)
}

# A fit of clmm()'s single-intercept routine is made again by the general one when its
# between-stratum SD is below this share of the latent residual SD. The value decides
# only whether the second routine is tried, not what counts as a boundary fit.
.maihda_clmm_refit_below <- 0.01

# SD of the stratum random intercept of a clmm fit; NA when it cannot be read.
maihda_clmm_stratum_sd <- function(model) {
  v <- tryCatch(as.numeric(ordinal::VarCorr(model)[["stratum"]][1, 1]),
                error = function(e) NA_real_)
  if (length(v) == 1L && is.finite(v) && v >= 0) sqrt(v) else NA_real_
}

# TRUE when clmm()'s routine for a single random intercept returned every stratum
# effect at exactly zero. That routine tests the effects' convergence before it
# computes their Hessian, so effects that pass at their starting value of zero come
# back with the Hessian as it was handed in. The Laplace approximation is then the
# likelihood of the model without the random effect, which does not depend on the
# between-stratum SD, and a quadrature rule is left unadapted: the SD of such a fit is
# where the optimiser stood, not an estimate. At clmm()'s own starting point (the
# pooled fixed-effect estimates, SD 1) every stratum's score is zero whenever the
# strata have the same outcome proportions or the fixed part spans the strata, and the
# fit then returns that SD of 1.
maihda_clmm_effects_at_start <- function(model) {
  if (!inherits(model, "clmm") || !isFALSE(model$control$useMatrix)) {
    return(FALSE)
  }
  effects <- model$ranef
  is.numeric(effects) && length(effects) > 0L && isTRUE(all(effects == 0))
}

# What the fit record says of a fit maihda_clmm_effects_at_start() recognises and
# clmm()'s general routine did not take over.
maihda_clmm_effects_at_start_message <- function(model) {
  n_quad <- model$nAGQ
  quadrature <- is.numeric(n_quad) && length(n_quad) == 1L && !n_quad %in% c(0, 1)
  paste0(
    "ordinal::clmm() returned every stratum effect at exactly zero beside a ",
    "between-stratum SD of ", format(signif(maihda_clmm_stratum_sd(model), 3)),
    ". Effects of exactly zero go with an SD of zero, so this SD was not estimated: ",
    "do not use it, or the VPC and PCV built on it.",
    if (quadrature) {
      paste0(" With nAGQ = 1 (the Laplace approximation) such a fit is made again by ",
             "clmm()'s general routine.")
    })
}

# TRUE when a clmm fit is to be made again by clmm()'s general routine
# (clmm.control(useMatrix = TRUE)): a Laplace fit of the single-intercept routine whose
# stratum effects are all exactly zero (maihda_clmm_effects_at_start()), or whose
# between-stratum SD is below .maihda_clmm_refit_below of the latent residual SD. The
# single-intercept routine works on the log of the SD, so on its way to a zero variance
# it stops where the optimiser's relative tolerance is met, at an SD that depends on
# the data. The general routine bounds the SD at zero and computes the effects' Hessian
# where it returns them, so it can return the boundary fit. A quadrature fit stays
# where it is: only the single-intercept routine offers quadrature.
maihda_clmm_refit_wanted <- function(model) {
  if (!inherits(model, "clmm") || !isFALSE(model$control$useMatrix)) {
    return(FALSE)
  }
  n_quad <- model$nAGQ
  if (!is.numeric(n_quad) || length(n_quad) != 1L || !n_quad %in% c(0, 1)) {
    return(FALSE)
  }
  if (maihda_clmm_effects_at_start(model)) {
    return(TRUE)
  }
  sd_stratum <- maihda_clmm_stratum_sd(model)
  sd_latent <- if (identical(model$link, "probit")) 1 else pi / sqrt(3)
  is.finite(sd_stratum) && sd_stratum / sd_latent < .maihda_clmm_refit_below
}

# Where the general routine starts when it fits the model of `model` again: the first
# fit's thresholds and coefficients, and an SD of zero. The likelihood is level in the
# SD at zero, so the routine stays on that bound and returns the boundary fit -- the
# first fit's fixed part with no between-stratum variance -- after a few iterations;
# maihda_clmm_refit_no_worse() then says whether that fit or the first one is kept.
# From its own start (the estimates of the model without the random effect, SD 1) the
# routine is a poor judge of a small variance: it can use up nlminb()'s iterations
# before it arrives, most readily with a covariate on a scale far from 1, and it can
# come to rest at zero where a small positive SD has the higher likelihood.
maihda_clmm_general_start <- function(model) {
  list(as.numeric(model$coefficients), 0)
}

# The call that fits the model of a first clmm() call by the general routine. `args`
# and `dot_vals` are what the first call was built from. Loose arguments of clmm() are
# control settings and take precedence over the entries of `control`, so the routine
# and its optimiser are set that way and everything else the caller passed stays as it
# was. Left out: the caller's starting value, which the two routines read in different
# forms, and the quadrature setting, which the general routine refuses unless it is 1.
# The start is named, not written out: the caller binds .maihda_clmm_start where the
# call is evaluated (see maihda_clmm_general_start()), which keeps the stored call
# short.
maihda_clmm_general_call <- function(args, dot_vals) {
  tags <- names(dot_vals)
  if (is.null(tags)) {
    tags <- rep("", length(dot_vals))
  }
  # `method` sits ahead of clmm.control()'s dots, so a shortened spelling binds to it.
  left_out <- tags %in% c("start", "nAGQ", "useMatrix") |
    (nzchar(tags) & startsWith("method", tags))
  as.call(c(list(quote(ordinal::clmm)), args, dot_vals[!left_out],
            list(start = quote(.maihda_clmm_start), useMatrix = TRUE,
                 method = "nlminb")))
}

# TRUE when the between-stratum variance of a clmm fit is below the 1e-8 that
# maihda_fit_diagnostics() flags as a singular fit.
maihda_clmm_at_boundary <- function(model) {
  sd_stratum <- maihda_clmm_stratum_sd(model)
  is.finite(sd_stratum) && sd_stratum^2 < 1e-8
}

# clmm()'s optimiser verdict for a fit: TRUE, FALSE, or NA where it is not known --
# there is no return code to read, or ucminf() could go no further and the fit cannot
# say where that left it. The return is read by the rules of the optimiser that made
# it: ucminf()'s by maihda_ucminf_converged(), any other as nlminb()'s, where 0 is
# convergence. nlminb() can also report "singular convergence" (PORT code 7) with the
# between-stratum variance at zero; that is the boundary fit, which the `singular` flag
# reports. And a boundary fit that maihda_clmm_general_refit() kept on the strength of
# the first fit's own convergence carries that verdict with it as a mark.
maihda_clmm_converged <- function(model) {
  if (isTRUE(attr(model, "maihda_converged"))) {
    return(TRUE)
  }
  returned <- tryCatch(model$optRes, error = function(e) NULL)
  code <- tryCatch(returned$convergence, error = function(e) NULL)
  if (!is.numeric(code) || length(code) != 1L || is.na(code)) {
    return(NA)
  }
  if (maihda_is_ucminf_return(returned)) {
    return(maihda_ucminf_converged(returned, model))
  }
  if (code == 0) {
    return(TRUE)
  }
  maihda_clmm_at_boundary(model) &&
    identical(maihda_nlminb_port_code(returned$message), 7L)
}

# TRUE when `returned`, the optimiser's return a clmm fit keeps as $optRes, is
# ucminf()'s: it has the class "ucminf", which older versions of ucminf do not set, or
# the largest gradient and the count of evaluations that ucminf() reports in $info.
# The fit's control settings do not say which optimiser ran: clmm()'s general routine
# runs nlminb() when it is asked for ucminf().
maihda_is_ucminf_return <- function(returned) {
  inherits(returned, "ucminf") ||
    (is.list(returned) &&
       all(c("maxgradient", "neval") %in% names(returned[["info"]])))
}

# A clmm fit counts as at a stationary point of its likelihood when it is less than
# this many standard errors from one (maihda_clmm_stationary_distance()).
.maihda_clmm_stationary_within <- 0.1

# ucminf()'s verdict for the fit `model` it returned as `returned`. It stops for one of
# four reasons: its gradient tolerance is met (code 1), a step is shorter than its step
# tolerance (2), its evaluations are used up (3), or its line search finds no lower
# point (4). 1 is convergence; 3 is not, and nor are the negative codes of a search
# that did not start. 2 and 4 say only that the search could go no further. That is
# where it ends at a maximum, once the gradient tolerance asks for more than the
# finite-difference gradient of the likelihood can show, and also where it ends short
# of one: with the starting values when its first line search fails, or some steps on
# when the likelihood is rough. So for 2 and 4 the fit itself is asked, and it is
# converged when it lies within .maihda_clmm_stationary_within standard errors of a
# stationary point. Where the fit cannot say, a search that made no step is not
# converged -- it returned its starting values -- and any other is unknown, NA.
maihda_ucminf_converged <- function(returned, model) {
  code <- returned$convergence
  if (code == 1) {
    return(TRUE)
  }
  if (!code %in% c(2, 4)) {
    return(FALSE)
  }
  distance <- maihda_clmm_stationary_distance(model)
  if (!is.na(distance)) {
    return(distance < .maihda_clmm_stationary_within)
  }
  if (maihda_ucminf_made_no_step(returned)) FALSE else NA
}

# TRUE when ucminf() made no step and so returned its starting values: its count of
# evaluations stands at the 1 made there and grows with every step.
maihda_ucminf_made_no_step <- function(returned) {
  count <- tryCatch(as.numeric(returned[["info"]][["neval"]]),
                    error = function(e) NA_real_)
  length(count) == 1L && is.finite(count) && count <= 1
}

# The step of the central differences by which clmm() takes the gradient and the
# Hessian it stores with a fit, the same in every parameter (the default of ordinal's
# deriv12()).
.maihda_clmm_derivative_step <- 1e-4

# How coarse that step is for each parameter of a fit with the Hessian `hessian`, whose
# inverse is `inverse`:
#   step^2 * H_jj^1.5 * sqrt((H^-1)_jj),
# the squared step in the parameter's conditional standard errors times the ratio of
# its marginal standard error to its conditional one. A central difference along a
# parameter errs by a sixth of the squared step times the third derivative there, and
# a few thousandths of this is what that error adds to the distance of
# maihda_clmm_stationary_distance(). For a parameter tied to no other it is the squared
# ratio of the step to its standard error; a covariate far from zero ties its
# coefficient to the thresholds, and a step along the coefficient then moves every
# observation's predictor by the step times the covariate.
maihda_clmm_step_coarseness <- function(hessian, inverse) {
  .maihda_clmm_derivative_step^2 * diag(hessian)^1.5 * sqrt(diag(inverse))
}

# How far a clmm fit is from a stationary point of its likelihood, in standard errors:
# sqrt(g' H^-1 g) for the gradient g and the Hessian H of the negative log-likelihood
# that clmm() stores with the fit. It is the length of the Newton step in the metric of
# H, and half its square is the log-likelihood the fit is short of that point by. NA
# when the fit cannot say: it holds no gradient and Hessian that go together, the
# Hessian is not positive definite, or the step by which clmm() took them is too
# coarse for a parameter (maihda_clmm_step_coarseness() above 1): a coefficient with a
# standard error below the step, as with a covariate in the thousands, or one tied to
# the thresholds by a covariate far from zero, as with a calendar year. What such a
# gradient and Hessian give is the error of the differences, not the fit's distance:
# they have put fits at their maximum up to 5 standard errors from it. A fit made on
# rescaled columns (maihda_clmm_carry_back()) is asked by the derivatives clmm() took
# there; the distance is the same in either parametrisation, the coarseness is not.
maihda_clmm_stationary_distance <- function(model) {
  taken <- tryCatch(attr(model, "maihda_rescaled"), error = function(e) NULL)
  if (is.list(taken)) {
    model <- taken
  }
  gradient <- tryCatch(as.numeric(model$gradient), error = function(e) NULL)
  hessian <- tryCatch(model$Hessian, error = function(e) NULL)
  if (!is.numeric(gradient) || length(gradient) == 0L || !is.matrix(hessian) ||
      !is.numeric(hessian) || !identical(dim(hessian), rep(length(gradient), 2L)) ||
      !all(is.finite(gradient)) || !all(is.finite(hessian))) {
    return(NA_real_)
  }
  curvature <- tryCatch(eigen(hessian, symmetric = TRUE, only.values = TRUE)$values,
                        error = function(e) NA_real_)
  if (!all(is.finite(curvature)) || min(curvature) <= 0) {
    return(NA_real_)
  }
  inverse <- tryCatch(solve(hessian), error = function(e) NULL)
  if (is.null(inverse) || !all(is.finite(inverse)) || min(diag(inverse)) <= 0) {
    return(NA_real_)
  }
  if (max(maihda_clmm_step_coarseness(hessian, inverse)) > 1) {
    return(NA_real_)
  }
  squared <- sum(gradient * drop(inverse %*% gradient))
  if (is.finite(squared) && squared >= 0) sqrt(squared) else NA_real_
}

# What the fit record says of a fit its optimiser did not vouch for: the optimiser's
# own message. Of a ucminf() search that could go no further it says where that left
# the fit: at the starting values when no step was made, else how far from a
# stationary point.
maihda_clmm_optimiser_message <- function(model) {
  returned <- tryCatch(model$optRes, error = function(e) NULL)
  said <- tryCatch(as.character(returned$message), error = function(e) character(0))
  said <- said[nzchar(said)]
  if (!maihda_is_ucminf_return(returned) ||
      !isTRUE(returned$convergence %in% c(2, 4))) {
    return(said)
  }
  own <- if (length(said) > 0L) paste0(" (", sub("[.]$", "", said[1]), ")")
  if (maihda_ucminf_made_no_step(returned)) {
    return(paste0("ucminf() made no step from the starting values", own,
                  ": the estimates returned are those values."))
  }
  distance <- maihda_clmm_stationary_distance(model)
  if (isTRUE(distance >= .maihda_clmm_stationary_within)) {
    return(paste0("ucminf() stopped ", format(signif(distance, 2)),
                  " standard errors from a stationary point by the fit's own gradient ",
                  "and Hessian", own, "."))
  }
  said
}

# The PORT return code nlminb() closes its message with, as in
# "singular convergence (7)"; NA for any other text.
maihda_nlminb_port_code <- function(message) {
  if (!is.character(message) || length(message) != 1L || is.na(message)) {
    return(NA_integer_)
  }
  found <- regmatches(message, regexec("\\(([0-9]+)\\)[[:space:]]*$", message))[[1]]
  if (length(found) == 2L) as.integer(found[2]) else NA_integer_
}

# TRUE when a second fit can take the first one's place: a clmm fit with a
# between-stratum SD and a log-likelihood to read, which converged or sits on the
# boundary. nlminb() cannot always say of the boundary fit that it converged: started
# on the fit it is to return, it reports "false convergence" with a covariate on a
# scale far from 1. Whether such a fit is kept is for maihda_clmm_refit_no_worse() to
# say, and what its fit record says of convergence for maihda_clmm_general_refit().
maihda_clmm_refit_usable <- function(refit) {
  if (!inherits(refit, "clmm") || !is.finite(maihda_clmm_stratum_sd(refit)) ||
      !isTRUE(is.finite(as.numeric(refit$logLik)))) {
    return(FALSE)
  }
  isTRUE(maihda_clmm_converged(refit)) || maihda_clmm_at_boundary(refit)
}

# Two log-likelihoods closer than this share of their size are not told apart: ten
# times the relative tolerance nlminb() stops at by default.
.maihda_clmm_loglik_tol <- 1e-9

# TRUE when the second fit's log-likelihood is not below the first one's. The second
# fit is the boundary fit, so a first fit with the higher log-likelihood has a
# between-stratum variance that is not zero, and stays. A first fit left with every
# stratum effect at zero carries the log-likelihood of the model without the random
# effect, which is the boundary fit's own.
maihda_clmm_refit_no_worse <- function(refit, model) {
  ll_new <- as.numeric(refit$logLik)
  ll_old <- as.numeric(model$logLik)
  if (length(ll_old) != 1L || !is.finite(ll_old)) {
    return(TRUE)
  }
  isTRUE(ll_new >= ll_old - .maihda_clmm_loglik_tol * max(1, abs(ll_old)))
}

# The fit of clmm()'s general routine for the model of `model`, or `model` itself when
# that fit fails, is no fit to take the first one's place (maihda_clmm_refit_usable())
# or has the lower log-likelihood. `signalled` holds the warnings and messages of the
# first fit, which the second would repeat; what the second fit has to say beyond
# them is passed on only if that fit is the one kept.
maihda_clmm_general_refit <- function(model, args, dot_vals, fit_env, signalled) {
  assign(".maihda_clmm_start", maihda_clmm_general_start(model), envir = fit_env)
  refit_call <- maihda_clmm_general_call(args, dot_vals)
  fresh <- list()
  refit <- tryCatch(
    withCallingHandlers(
      eval(refit_call, fit_env),
      warning = function(w) {
        if (!conditionMessage(w) %in% signalled) {
          fresh[[length(fresh) + 1L]] <<- w
        }
        invokeRestart("muffleWarning")
      },
      message = function(m) {
        if (!conditionMessage(m) %in% signalled) {
          fresh[[length(fresh) + 1L]] <<- m
        }
        invokeRestart("muffleMessage")
      }),
    error = function(e) NULL)
  if (!maihda_clmm_refit_usable(refit) || !maihda_clmm_refit_no_worse(refit, model)) {
    return(model)
  }
  # A boundary fit nlminb() did not vouch for takes the first fit's convergence, when
  # that fit converged: it started from the same fixed part and loses nothing. The fit
  # record reads the mark; without it, it gives the optimiser's own return.
  if (!isTRUE(maihda_clmm_converged(refit)) && isTRUE(maihda_clmm_converged(model))) {
    attr(refit, "maihda_converged") <- TRUE
  }
  for (cond in fresh) {
    if (inherits(cond, "warning")) warning(cond) else message(cond)
  }
  refit
}

# The fit of `fit_call`, made again by clmm()'s general routine where its routine for a
# single random intercept leaves it at or near a between-stratum SD of zero (see
# maihda_clmm_refit_wanted()). What the first fit signals is noted on its way out, so
# that the second does not say it again; `signalled` holds what was said before it.
maihda_clmm_fit_and_refit <- function(fit_call, args, dot_vals, fit_env,
                                      signalled = character(0)) {
  model <- withCallingHandlers(
    eval(fit_call, fit_env),
    warning = function(w) signalled <<- c(signalled, conditionMessage(w)),
    message = function(m) signalled <<- c(signalled, conditionMessage(m)))
  if (maihda_clmm_refit_wanted(model)) {
    model <- maihda_clmm_general_refit(model, args, dot_vals, fit_env, signalled)
  }
  model
}

# The fit of the model of `args` made on its fixed-effects design with every column
# centred and scaled (maihda_clmm_rescaled_design()), the boundary fit included, and
# returned in the formula's own terms (maihda_clmm_carry_back()). NULL when the caller
# gave a start, the design has no column to rescale, or the fit could not be made or
# carried back; the model is then fitted as the formula stands. What the fits signal
# is held back until the fit to be returned is known, and is then passed on once.
maihda_clmm_rescaled_fit <- function(args, dot_vals, fit_env, fit_call) {
  heard <- list()
  outcome <- tryCatch(
    withCallingHandlers({
      design <- maihda_clmm_rescaled_design(args, dot_vals, fit_env)
      if (is.null(design)) {
        NULL
      } else {
        said <- vapply(heard, conditionMessage, character(1))
        model <- maihda_clmm_fit_and_refit(design$call, design$args, design$dot_vals,
                                           design$env, said)
        maihda_clmm_carry_back(model, design, args, dot_vals, fit_env, fit_call)
      }
    },
    warning = function(w) {
      heard[[length(heard) + 1L]] <<- w
      invokeRestart("muffleWarning")
    },
    message = function(m) {
      heard[[length(heard) + 1L]] <<- m
      invokeRestart("muffleMessage")
    }),
    error = function(e) NULL)
  if (is.null(outcome)) {
    return(NULL)
  }
  for (cond in heard) {
    if (inherits(cond, "warning")) warning(cond) else message(cond)
  }
  outcome
}

# The model of `args` on its fixed-effects design with the columns centred and scaled:
# the call that fits it and what maihda_clmm_carry_back() needs to take that fit back.
# The design is clmm()'s own (control method "model.frame"): its columns, the intercept
# left out, become .maihda_x1, .maihda_x2, ... in a copy of the data, beside the
# formula's offset and its random-effect terms. A column is centred on its mean, which
# the thresholds take up, and divided by its SD; a constant column is left as it is,
# and under threshold = "symmetric2", whose thresholds cannot take up a shift, no column
# is centred. NULL when the caller gave a start or there is nothing to rescale.
maihda_clmm_rescaled_design <- function(args, dot_vals, fit_env) {
  tags <- names(dot_vals)
  if (is.null(tags)) {
    tags <- rep("", length(dot_vals))
  }
  if ("start" %in% tags) {
    return(NULL)
  }
  # `method` sits ahead of clmm.control()'s dots, so a shortened spelling binds to it.
  frame_call <- as.call(c(list(quote(ordinal::clmm)), args,
                          dot_vals[!(nzchar(tags) & startsWith("method", tags))],
                          list(method = "model.frame")))
  frames <- eval(frame_call, fit_env)
  data <- get("data", envir = fit_env)
  design <- frames$X
  if (!is.matrix(design) || nrow(design) != nrow(data)) {
    return(NULL)
  }
  design <- design[, colnames(design) != "(Intercept)", drop = FALSE]
  spread <- apply(design, 2L, stats::sd)
  varies <- is.finite(spread) & spread > 0
  if (!any(varies)) {
    return(NULL)
  }
  threshold <- dot_vals[["threshold"]]
  threshold <- if (is.null(threshold)) "flexible" else
    match.arg(threshold, c("flexible", "symmetric", "symmetric2", "equidistant"))
  centre <- if (identical(threshold, "symmetric2")) rep(0, ncol(design)) else
    colMeans(design)
  centre[!varies] <- 0
  spread[!varies] <- 1
  generated <- paste0(".maihda_x", seq_len(ncol(design)))
  if (any(c(generated, ".maihda_offset") %in% names(data))) {
    return(NULL)
  }
  for (j in seq_len(ncol(design))) {
    data[[generated[j]]] <- (design[, j] - centre[j]) / spread[j]
  }
  terms <- lapply(generated, as.name)
  if (is.numeric(frames$off) && any(frames$off != 0)) {
    data$.maihda_offset <- frames$off
    terms <- c(terms, list(quote(offset(.maihda_offset))))
  }
  terms <- c(terms, lapply(reformulas::findbars(args$formula), function(bar) call("(", bar)))
  formula <- stats::as.formula(call("~", args$formula[[2L]],
                                    maihda_rhs_from_plus_terms(terms)))
  env <- new.env(parent = parent.env(fit_env))
  env$data <- data
  environment(formula) <- env
  args$formula <- formula
  # The rescaled design has no factor left for a `contrasts` argument to code.
  dot_vals <- dot_vals[tags != "contrasts"]
  list(call = as.call(c(list(quote(ordinal::clmm)), args, dot_vals)), args = args,
       dot_vals = dot_vals, env = env, columns = colnames(design), generated = generated,
       centre = unname(centre), spread = unname(spread))
}

# The map from the fixed-effect parameters of a fit on rescaled columns to those of the
# formula's own columns: (alpha, beta) = map %*% (alpha_s, beta_s). Column j entered as
# (x_j - centre_j) / spread_j, so beta_j = beta_s,j / spread_j, and every cut point moves
# by sum(centre * beta); `tJac` (cut points = tJac %*% alpha) says how the free threshold
# parameters take up that common shift. Pure in its inputs.
maihda_clmm_rescaled_map <- function(tJac, centre, spread) {
  n_alpha <- ncol(tJac)
  n_beta <- length(spread)
  map <- diag(n_alpha + n_beta)
  if (n_beta == 0L) {
    return(map)
  }
  shift <- rep(0, n_alpha)
  if (any(centre != 0)) {
    shift <- qr.solve(tJac, rep(1, nrow(tJac)))
    if (max(abs(drop(tJac %*% shift) - 1)) > 1e-8) {
      stop("These thresholds cannot take up a shift of the cut points.", call. = FALSE)
    }
  }
  beta <- n_alpha + seq_len(n_beta)
  map[seq_len(n_alpha), beta] <- outer(shift, centre / spread)
  map[beta, beta] <- diag(1 / spread, n_beta)
  map
}

# The call that rebuilds the fit of the formula at a point it is handed, without a
# search: clmm()'s start bound as .maihda_clmm_start, nlminb() allowed no iteration, and
# no Hessian (maihda_clmm_carry_back() carries one over). `general`: clmm()'s general
# routine, which takes no quadrature setting.
maihda_clmm_rebuild_call <- function(args, dot_vals, general) {
  tags <- names(dot_vals)
  if (is.null(tags)) {
    tags <- rep("", length(dot_vals))
  }
  left_out <- tags %in% c("start", "useMatrix", "Hess", "iter.max", if (general) "nAGQ") |
    (nzchar(tags) & startsWith("method", tags))
  args$Hess <- FALSE
  as.call(c(list(quote(ordinal::clmm)), args, dot_vals[!left_out],
            list(start = quote(.maihda_clmm_start), method = "nlminb", iter.max = 0L),
            if (general) list(useMatrix = TRUE)))
}

# The clmm() fit of the formula at the maximum of `model`, the fit on the rescaled
# columns of `design`. clmm() builds it at the carried-back point without a search, so
# its frames, coefficients, stratum effects and fitted values are those of the formula;
# the record of the search is the rescaled fit's -- the optimiser's return (its `par`
# carried back), the control settings, the start, the counts of evaluations and
# iterations -- and so are the gradient and Hessian, carried over exactly
# (g = T^-T g_s, H = T^-T H_s T^-1 for the map T of the free parameters), which give its
# standard errors. The derivatives as clmm() took them, on the rescaled columns, stay
# with the fit as attr(, "maihda_rescaled") for maihda_clmm_stationary_distance(). An
# error when the two fits do not line up.
maihda_clmm_carry_back <- function(model, design, args, dot_vals, fit_env, fit_call) {
  n_fixed <- length(model$coefficients)
  kept <- match(names(model$beta), design$generated)
  if (anyNA(kept) || !is.matrix(model$Hessian) || !is.numeric(model$gradient)) {
    stop("The fit on rescaled columns cannot be carried back.", call. = FALSE)
  }
  map <- maihda_clmm_rescaled_map(model$tJac, design$centre[kept], design$spread[kept])
  fixed <- seq_len(n_fixed)
  carried <- function(par) {
    par <- as.numeric(par)
    par[fixed] <- drop(map %*% par[fixed])
    par
  }
  par <- carried(model$optRes$par)
  general <- isTRUE(model$control$useMatrix)
  assign(".maihda_clmm_start", if (general) list(par[fixed], par[-fixed]) else par,
         envir = fit_env)
  rebuilt <- suppressWarnings(suppressMessages(
    eval(maihda_clmm_rebuild_call(args, dot_vals, general), fit_env)))
  free <- length(model$gradient)
  if (!inherits(rebuilt, "clmm") || length(rebuilt$coefficients) != n_fixed ||
      !identical(unname(names(rebuilt$beta)), design$columns[kept]) ||
      length(rebuilt$gradient) != free || free < n_fixed) {
    stop("The fit on rescaled columns cannot be carried back.", call. = FALSE)
  }
  # The free parameters are the fixed part and any variance parameter off its bound,
  # which the map leaves as they are.
  carry <- diag(free)
  carry[fixed, fixed] <- map
  inverse <- solve(carry)
  gradient <- drop(crossprod(inverse, as.numeric(model$gradient)))
  hessian <- crossprod(inverse, model$Hessian %*% inverse)
  names(gradient) <- names(rebuilt$gradient)
  dimnames(hessian) <- list(names(gradient), names(gradient))
  rebuilt$gradient <- gradient
  rebuilt$Hessian <- hessian
  returned <- model$optRes
  returned$par <- par
  rebuilt$optRes <- returned
  rebuilt$control <- model$control
  start <- model$start
  rebuilt$start <- if (is.list(start)) list(carried(start[[1L]]), start[[2L]]) else
    carried(start)
  rebuilt$Niter <- model$Niter
  rebuilt$neval <- model$neval
  rebuilt$info$niter <- model$info$niter
  rebuilt$info$max.grad <- formatC(max(abs(gradient)), digits = 2, format = "e")
  if (isTRUE(attr(model, "maihda_converged"))) {
    attr(rebuilt, "maihda_converged") <- TRUE
  }
  attr(rebuilt, "maihda_rescaled") <- list(
    columns = design$columns, centre = design$centre, spread = design$spread,
    gradient = model$gradient, Hessian = model$Hessian)
  # The call the fit stands for, as clmm() records a call it is given: the general
  # routine's for the boundary fit maihda_clmm_general_refit() made, whose call starts
  # from .maihda_clmm_start, and the caller's otherwise.
  refitted <- identical(model$call[["start"]], quote(.maihda_clmm_start))
  rebuilt$call <- match.call(ordinal::clmm,
                             if (refitted) maihda_clmm_general_call(args, dot_vals) else
                               fit_call)
  rebuilt
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
#' constructed with the fit's factor coding -- the levels it fitted and the contrast
#' matrix it applied to each factor, whatever \code{options(contrasts = )} says now
#' -- AND its transformation basis (so a
#' data-dependent term such as \code{scale(x)} uses the fit's centre and scale rather
#' than recomputing them from \code{newdata}) and multiplied by the
#' location coefficients \code{beta} (a clmm has \emph{no} intercept column --
#' it is absorbed by the thresholds -- so \code{beta}'s names select the right
#' columns), any formula offset term is evaluated on \code{newdata} and added
#' (so an offset-only null model still predicts its offset), and
#' \code{include_re} adds each row's stratum conditional mode (an
#' unseen stratum contributes 0 -- the zero-effect fallback that
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
  # recomputed from the prediction batch (see maihda_fitted_predict_terms()), and the
  # factors coded exactly as the fit coded them -- its levels and contrast matrices --
  # so a sum, custom or option-set contrast is not silently re-coded by default
  # treatment coding (see maihda_engine_fixed_coding()).
  tt <- maihda_fitted_predict_terms(object$formula, object$data)$terms
  design <- maihda_fixed_design(tt, newdata, maihda_object_fixed_coding(object))
  mf <- design$frame

  if (is.null(beta) || length(beta) == 0) {
    # Null (thresholds-only) model: the location fixed part is identically 0
    # (the offset, if any, is added below).
    eta <- rep(0, nrow(newdata))
  } else {
    X <- design$X
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
#' The statistic is the omnibus nominal-effects likelihood-ratio statistic: the
#' fixed-effect part of the model is refitted twice with \code{ordinal::clm()} --
#' once with all covariates proportional, once with every covariate entering as a
#' threshold-specific (nominal) effect -- and twice the log-likelihood difference
#' is taken.
#'
#' Its usual chi-squared reference (\code{ordinal::nominal_test()}) is not valid
#' here. After the stratum random intercept is integrated out, the marginal
#' cumulative-logit slopes generally differ across thresholds once the stratum
#' variance is non-zero. The chi-squared reference also treats the observations
#' as independent, whereas observations that share a stratum share its random
#' effect. A correctly specified model is therefore rejected too often, and more
#' often as the sample grows.
#'
#' The null distribution is instead simulated under the fitted \code{clmm}: each
#' replicate redraws the stratum random effects from \eqn{N(0, \tau^2)} at the
#' fitted variance, forms the conditional category probabilities from the fitted
#' thresholds and location predictor, redraws the response, and recomputes the
#' statistic. The p-value is
#' \eqn{(1 + \#\{T_b \ge T_{obs}\}) / (1 + B)}{(1 + #{T_b >= T_obs}) / (1 + B)}
#' over the \eqn{B} replicates that refitted. \code{n_sim} sets only its Monte
#' Carlo resolution: the smallest attainable value is \eqn{1 / (B + 1)}.
#'
#' This is a parametric-bootstrap approximation, not an exact test. The replicates
#' are drawn at the fitted stratum variance, which is estimated from only as many
#' units as there are strata, so the test is least reliable when the strata are
#' few, whatever \code{n_sim} is.
#'
#' The test is expensive -- every replicate refits two \code{clm()} models -- and
#' is not run at fit time.
#'
#' @param object A \code{maihda_model} fitted with \code{engine = "ordinal"}.
#' @param n_sim Number of parametric-bootstrap replicates (default 199).
#' @param seed Optional integer seed, for a reproducible bootstrap.
#' @return An object of class \code{maihda_po_test}: a list with \code{lrt},
#'   \code{df}, \code{n_terms}, \code{p_value} (the bootstrap p-value),
#'   \code{p_chisq}, \code{n_sim} (replicates that produced a usable statistic),
#'   \code{n_failed}, and \code{null_lrt} (the simulated null statistics).
#'   \code{p_value} is the one to report and the only one printed. \code{p_chisq}
#'   is the chi-squared p-value of the same statistic; it is not calibrated for a
#'   mixed model (see Details) and is returned for comparison only.
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

  # The observed statistic and every bootstrap replicate must be computed on the
  # SAME design with only the response redrawn, so build it once and reuse it.
  # maihda_po_refit_design() supplies the fitted fixed-effects columns under
  # generated names, which is what makes a transformed term such as log(x) or
  # factor(g) refittable at all -- see its comment for why reformulating the
  # original term labels against the clmm model frame cannot work.
  design <- maihda_po_refit_design(object$model)
  if (is.null(design)) {
    stop("A null, covariate-free cumulative model has no covariate slopes to ",
         "test for proportional odds.", call. = FALSE)
  }
  obs <- maihda_po_lrt(design$data, design$resp, design$terms)
  if (is.null(obs)) {
    # Distinct from the covariate-free case above: the covariates are there, but
    # the two fixed-only cumulative refits did not yield a usable statistic.
    # Measured causes are degenerate designs -- an exactly collinear pair of
    # covariates, a constant covariate -- where the refits mostly CONVERGE yet
    # their likelihoods give a non-positive df or a negative LRT, so
    # maihda_po_lrt()'s own sanity check rejects the pair. (An earlier version of
    # this comment blamed non-monotone thresholds under a continuous nominal
    # effect; that was never demonstrated and is not what the probes showed.)
    stop("The proportional-odds statistic could not be computed for this fit: ",
         "the fixed-only cumulative refit it is built from failed to converge.",
         call. = FALSE)
  }
  obs$n_terms <- design$n_terms

  # Ingredients of the fitted conditional model. maihda_ordinal_check_formula()
  # guarantees the clmm path is the canonical single (1 | stratum) structure, so
  # one random-effect SD and one grouping factor are all that is needed.
  cutpoints <- maihda_clmm_cutpoints(object$model)
  link <- object$family$link
  tau <- sqrt(max(maihda_clmm_variances(object)$stratum, 0))
  eta_fixed <- maihda_clmm_linpred(object, include_re = FALSE)

  dat <- design$data
  resp <- design$resp
  fixed_terms <- design$terms

  grp <- factor(as.character(object$data$stratum))
  if (nrow(dat) != length(eta_fixed) ||
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
