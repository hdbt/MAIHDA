#' Create Strata from Multiple Variables
#'
#' This function creates strata (intersectional categories) from multiple
#' categorical variables in a dataset.
#'
#' @param data A data frame containing the variables to create strata from.
#' @param vars Character vector of variable names to use for creating strata.
#' @param sep Separator to use between variable values when creating stratum labels.
#'   Default is " \\u00d7 " (a mathematical multiplication sign).
#' @param min_n Minimum number of observations required for a stratum to be included.
#'   Strata with fewer observations will be coded as NA. Default is 1.
#' @param autobin Logical indicating whether to automatically bin numeric grouping
#'   variables with more than 10 unique values into 3 categories (tertiles).
#'   Default is TRUE. When this happens a \code{message()} is emitted, because the
#'   resulting strata are data-dependent (tertile cut-points depend on the sample)
#'   and a continuous variable placed in the grouping term is usually unintended.
#'   Set \code{autobin = FALSE} to disable, or bin the variable yourself for
#'   explicit, reproducible cut-points.
#' @param bin_rows Optional logical vector, one element per row of \code{data},
#'   selecting the rows used to compute the auto-bin cut-points (every row is still
#'   assigned a stratum). Mainly for internal use by \code{\link{fit_maihda}}, which
#'   passes the \code{subset} rows so that \code{fit_maihda(..., subset = keep)}
#'   bins on the same sample as fitting \code{data[keep, ]}. Default \code{NULL}
#'   uses all rows.
#'
#' @return An object of class \code{maihda_strata}: a list with elements
#'   \item{data}{The original data frame with an added 'stratum' column. The
#'     strata_info is also attached as an attribute for use by fit_maihda()}
#'   \item{strata_info}{A data frame with information about each stratum including
#'     counts and the combination of variable values}
#'   \item{vars}{The stratum-defining variable names, as supplied.}
#'   \item{sep}{The label separator used.}
#'   \item{min_n}{The minimum stratum size applied.}
#'   \item{autobin_info}{A named list with the \code{breaks} and \code{labels}
#'     of each auto-binned numeric variable (empty when nothing was binned).}
#'
#' @details
#' If any of the specified variables has a missing value (NA) for a given observation,
#' that observation will be assigned to the NA stratum (stratum = NA), rather than
#' creating a stratum that includes the missing value.
#'
#' The strata_info data frame is also attached as an attribute to the data, which
#' allows fit_maihda() to automatically capture stratum labels for use in plots
#' and summaries.
#'
#' When \code{autobin} discretises a numeric grouping variable \code{v}, the
#' adjusted-model and prediction machinery later add an internal factor column named
#' \code{.maihda_dim_<v>}; the \code{.maihda_dim_} prefix is therefore reserved.
#' \code{make_strata()} errors if \code{data} already holds the \code{.maihda_dim_<v>}
#' column for a variable it is about to auto-bin, so an existing user column is never
#' silently overwritten (rename it, or pass \code{autobin = FALSE}).
#'
#' A numeric grouping variable with 10 or fewer unique values (or any numeric when
#' \code{autobin = FALSE}) is kept as its raw values: distinct values still define
#' distinct strata, but the variable is \emph{not} binned. If such a variable is then
#' used to build an adjusted MAIHDA (\code{\link{maihda}}), it enters the adjusted model
#' as a linear fixed effect rather than categorical main effects, and \code{maihda()}
#' warns (when it has three or more distinct values). Wrap category codes in
#' \code{factor()} before creating strata if you want them treated categorically.
#'
#' @examples
#' # Create strata from gender and race variables
#' result <- make_strata(maihda_sim_data, vars = c("gender", "race"))
#' print(result$strata_info)
#'
#' @export
#' @importFrom dplyr mutate group_by summarise n ungroup
#' @importFrom tidyr unite
#' @importFrom rlang .data
#' @importFrom stats quantile na.omit
make_strata <- function(data, vars, sep = " \u00d7 ", min_n = 1, autobin = TRUE,
                        bin_rows = NULL) {
  # Input validation
  if (!is.data.frame(data)) {
    stop("'data' must be a data frame")
  }
  # `bin_rows` (internal): a logical mask over the rows of `data` selecting the
  # sample used to COMPUTE the auto-bin cut-points; assignment via cut() still
  # covers every row. fit_maihda() passes the `subset` rows so that
  # fit_maihda(..., subset = keep) auto-bins on the same rows as fitting
  # data[keep, ] -- otherwise excluded rows would move the tertile cut-points.
  if (!is.null(bin_rows)) {
    if (!is.logical(bin_rows) || length(bin_rows) != nrow(data)) {
      stop("'bin_rows' must be a logical vector with one element per row of 'data'.",
           call. = FALSE)
    }
    if (!any(bin_rows, na.rm = TRUE)) {
      bin_rows <- NULL
    }
  }

  if (nrow(data) == 0) {
    stop("'data' has no rows; strata cannot be created from an empty data ",
         "frame.", call. = FALSE)
  }

  if (!is.character(vars) || length(vars) == 0) {
    stop("'vars' must be a character vector with at least one variable name")
  }

  # A dimension may define the strata only once: c("a", "a") would otherwise
  # build degenerate cells with labels like "a_val x a_val" and duplicate
  # metadata dimensions, silently mis-specifying the intersection.
  if (anyDuplicated(vars)) {
    stop("'vars' contains duplicated variable name(s): ",
         paste(unique(vars[duplicated(vars)]), collapse = ", "),
         ". Each dimension can define the strata only once.", call. = FALSE)
  }

  missing_vars <- setdiff(vars, names(data))
  if (length(missing_vars) > 0) {
    stop("Variables not found in data: ", paste(missing_vars, collapse = ", "))
  }

  if (!is.numeric(min_n) || length(min_n) != 1 ||
      is.na(min_n) || !is.finite(min_n) ||
      min_n < 1 || min_n != floor(min_n)) {
    stop("'min_n' must be a single positive whole number.", call. = FALSE)
  }
  min_n <- as.integer(min_n)

  # Create a copy of the data to avoid modifying the original. Numeric
  # auto-binning is applied only to the temporary strata-building columns so
  # original variables remain valid for fixed-effect model terms.
  result_data <- data
  strata_data <- data[, vars, drop = FALSE]

  # Auto-bin numeric variables with >10 unique values into 3 categories
  autobin_info <- list()
  if (autobin) {
    for (v in vars) {
      val <- strata_data[[v]]
      # Cut-points come from `bin_val` (the bin_rows sample, or all rows when
      # bin_rows is NULL); cut() below still bins every row. Whether the dimension
      # is auto-binned at all is likewise decided on that sample, so it matches
      # fitting the already-subset data.
      bin_val <- if (is.null(bin_rows)) val else val[bin_rows]
      if (is.numeric(val) && length(unique(stats::na.omit(bin_val))) > 10) {
        # This numeric dimension will be discretised, so the adjusted-model and
        # prediction machinery later add the reserved '.maihda_dim_<v>' factor column.
        # Reject a pre-existing user column of that name now -- before any
        # augmentation -- rather than silently overwriting it downstream.
        maihda_guard_reserved_dim_col(v, data)
        q <- stats::quantile(bin_val, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE)
        labels <- c(paste0(v, "_Low"), paste0(v, "_Mid"), paste0(v, "_High"))
        tertiles_ok <- length(unique(q)) == 4
        if (tertiles_ok) {
          breaks <- as.numeric(q)
        } else {
          # Tied quantiles (e.g. skewed/zero-inflated data): tertiles are not
          # defined, so fall back to equal-width bins. These are NOT tertiles and
          # can be highly imbalanced; warn rather than silently mislabel.
          rx <- range(bin_val, na.rm = TRUE)
          dx <- diff(rx)
          breaks <- seq(rx[1] - dx/1000, rx[2] + dx/1000, length.out = 4)
          warning("make_strata(): numeric variable '", v, "' has tied tertile ",
                  "cut-points, so equal-width bins were used instead of tertiles. ",
                  "The resulting groups may be very imbalanced; consider binning '",
                  v, "' yourself or setting autobin = FALSE.", call. = FALSE)
        }
        strata_data[[v]] <- cut(val, breaks = breaks, include.lowest = TRUE,
                                labels = labels)
        autobin_info[[v]] <- list(breaks = breaks, labels = labels)
        # Inform the user a continuous grouping variable was discretised: the
        # resulting strata are data-dependent. Pass autobin = FALSE to disable,
        # or pre-bin the variable yourself for explicit control.
        kind <- if (tertiles_ok) "tertiles" else "equal-width bins (tied tertiles)"
        counts <- table(strata_data[[v]])
        message("make_strata(): auto-binned numeric variable '", v, "' into ", kind,
                " (", paste(sprintf("%s=%d", names(counts), as.integer(counts)),
                            collapse = ", "), "). Set autobin = FALSE to disable.")
      }
    }
  }

  # Identify rows with any missing values in the specified variables
  has_missing <- apply(strata_data, 1, function(x) any(is.na(x)))

  # Build strata from the actual variable columns, not from pasted display
  # labels. This avoids collapsing distinct combinations whose values contain
  # the display separator.
  complete_strata_data <- strata_data[!has_missing, , drop = FALSE]
  unique_strata <- unique(complete_strata_data)
  combo_ids <- maihda_match_strata_rows(complete_strata_data, unique_strata, vars)
  stratum_counts <- tabulate(combo_ids, nbins = nrow(unique_strata))

  # Filter strata based on minimum count
  valid_idx <- which(stratum_counts >= min_n)

  # Create numeric stratum ID
  result_data$stratum <- NA_integer_

  # Assign stratum IDs only to rows without missing values that meet minimum count
  if (length(valid_idx) > 0) {
    result_data$stratum[!has_missing] <- match(combo_ids, valid_idx)
  }

  valid_strata <- unique_strata[valid_idx, , drop = FALSE]
  # Labels are built per column via as.character() (maihda_paste_label_rows), not
  # apply(): apply()'s as.matrix() coercion format()-pads numeric columns in a
  # mixed-type frame, producing labels like "m x  1" next to "m x 10".
  labels <- if (nrow(valid_strata) > 0) {
    maihda_paste_label_rows(valid_strata, sep)
  } else {
    character()
  }
  duplicated_labels <- duplicated(labels) | duplicated(labels, fromLast = TRUE)
  if (any(duplicated_labels)) {
    dup_rows <- valid_strata[duplicated_labels, , drop = FALSE]
    for (v in vars) {
      dup_rows[[v]] <- paste0(v, "=", as.character(dup_rows[[v]]))
    }
    labels[duplicated_labels] <- maihda_paste_label_rows(dup_rows, sep)
  }

  # Create stratum information table
  strata_info <- data.frame(
    stratum = seq_along(valid_idx),
    label = labels,
    n = as.integer(stratum_counts[valid_idx]),
    stringsAsFactors = FALSE
  )

  # Add the stratum-defining values without parsing the display label.
  for (var in vars) {
    strata_info[[var]] <- valid_strata[[var]]
  }

  # Attach strata_info as an attribute to the data for easy access
  attr(result_data, "strata_info") <- strata_info
  attr(result_data, "strata_vars") <- vars
  attr(result_data, "strata_sep") <- sep
  attr(result_data, "strata_autobin_info") <- autobin_info

  # Return results
  structure(
    list(
      data = result_data,
      strata_info = strata_info,
      vars = vars,
      sep = sep,
      min_n = min_n,
      autobin_info = autobin_info
    ),
    class = "maihda_strata"
  )
}

#' Print method for maihda_strata objects
#'
#' @param x A maihda_strata object
#' @param ... Additional arguments (not used)
#' @return No return value, called for side effects.
#' @export
print.maihda_strata <- function(x, ...) {
  pal <- maihda_palette()
  cat(pal$bold("MAIHDA Strata Object"), "\n", sep = "")
  cat("====================\n\n")
  cat("Variables used:", paste(x$vars, collapse = ", "), "\n")
  cat("Number of strata:", nrow(x$strata_info), "\n")
  cat("Total observations:", nrow(x$data), "\n")
  cat("Observations with valid strata:", sum(!is.na(x$data$stratum)), "\n\n")
  cat(pal$bold("Stratum summary:"), "\n", sep = "")
  print(x$strata_info)
  invisible(x)
}
