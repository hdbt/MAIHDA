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
#' @param autobin Logical indicating whether to automatically bin numeric grouping variables
#'   with more than 10 unique values into 3 categories (tertiles). Default is TRUE.
#'
#' @return A list with two elements:
#'   \item{data}{The original data frame with an added 'stratum' column. The
#'     strata_info is also attached as an attribute for use by fit_maihda()}
#'   \item{strata_info}{A data frame with information about each stratum including
#'     counts and the combination of variable values}
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
make_strata <- function(data, vars, sep = " \u00d7 ", min_n = 1, autobin = TRUE) {
  # Input validation
  if (!is.data.frame(data)) {
    stop("'data' must be a data frame")
  }

  if (!is.character(vars) || length(vars) == 0) {
    stop("'vars' must be a character vector with at least one variable name")
  }

  missing_vars <- setdiff(vars, names(data))
  if (length(missing_vars) > 0) {
    stop("Variables not found in data: ", paste(missing_vars, collapse = ", "))
  }

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
      if (is.numeric(val) && length(unique(stats::na.omit(val))) > 10) {
        q <- stats::quantile(val, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE)
        labels <- c(paste0(v, "_Low"), paste0(v, "_Mid"), paste0(v, "_High"))
        if (length(unique(q)) == 4) {
          breaks <- as.numeric(q)
        } else {
          rx <- range(val, na.rm = TRUE)
          dx <- diff(rx)
          breaks <- seq(rx[1] - dx/1000, rx[2] + dx/1000, length.out = 4)
        }
        strata_data[[v]] <- cut(val, breaks = breaks, include.lowest = TRUE,
                                labels = labels)
        autobin_info[[v]] <- list(breaks = breaks, labels = labels)
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
  labels <- if (nrow(valid_strata) > 0) {
    apply(valid_strata, 1, function(x) paste(x, collapse = sep))
  } else {
    character()
  }
  duplicated_labels <- duplicated(labels) | duplicated(labels, fromLast = TRUE)
  if (any(duplicated_labels)) {
    labels[duplicated_labels] <- apply(
      valid_strata[duplicated_labels, , drop = FALSE],
      1,
      function(x) paste(paste0(vars, "=", x), collapse = sep)
    )
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
  cat("MAIHDA Strata Object\n")
  cat("====================\n\n")
  cat("Variables used:", paste(x$vars, collapse = ", "), "\n")
  cat("Number of strata:", nrow(x$strata_info), "\n")
  cat("Total observations:", nrow(x$data), "\n")
  cat("Observations with valid strata:", sum(!is.na(x$data$stratum)), "\n\n")
  cat("Stratum summary:\n")
  print(x$strata_info)
  invisible(x)
}
