#' Create Strata from Multiple Variables
#'
#' This function creates strata (intersectional categories) from multiple
#' categorical variables in a dataset.
#'
#' @param data A data frame containing the variables to create strata from.
#' @param vars Character vector of variable names to use for creating strata.
#' @param sep Separator to use between variable values when creating stratum labels.
#'   Default is "_".
#' @param min_n Minimum number of observations required for a stratum to be included.
#'   Strata with fewer observations will be coded as NA. Default is 1.
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
#' \dontrun{
#' # Create strata from gender and race variables
#' data <- data.frame(
#'   gender = c("M", "F", "M", "F"),
#'   race = c("White", "Black", "White", "Black"),
#'   outcome = c(1, 2, 3, 4)
#' )
#' result <- make_strata(data, vars = c("gender", "race"))
#' print(result$strata_info)
#' }
#'
#' @export
#' @importFrom dplyr mutate group_by summarise n ungroup
#' @importFrom tidyr unite
#' @importFrom rlang .data
make_strata <- function(data, vars, sep = "_", min_n = 1) {
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
  
  # Create a copy of the data to avoid modifying the original
  result_data <- data
  
  # Identify rows with any missing values in the specified variables
  has_missing <- apply(result_data[, vars, drop = FALSE], 1, function(x) any(is.na(x)))
  
  # Create stratum variable by combining the specified variables
  # Only for rows without missing values
  result_data$stratum_label <- NA_character_
  result_data$stratum_label[!has_missing] <- apply(
    result_data[!has_missing, vars, drop = FALSE], 1, 
    function(x) paste(x, collapse = sep)
  )
  
  # Count observations per stratum (excluding rows with missing values)
  stratum_counts <- table(result_data$stratum_label[!has_missing])
  
  # Filter strata based on minimum count
  valid_strata <- names(stratum_counts[stratum_counts >= min_n])
  
  # Create numeric stratum ID
  result_data$stratum <- NA_integer_
  
  # Assign stratum IDs only to rows without missing values that meet minimum count
  if (length(valid_strata) > 0) {
    result_data$stratum[!has_missing] <- as.integer(
      factor(result_data$stratum_label[!has_missing], levels = valid_strata)
    )
  }
  
  # Create stratum information table
  strata_info <- data.frame(
    stratum = seq_along(valid_strata),
    label = valid_strata,
    n = as.integer(stratum_counts[valid_strata])
  )
  
  # Add the original variable values to strata_info
  if (nrow(strata_info) > 0) {
    for (var in vars) {
      strata_info[[var]] <- sapply(strsplit(strata_info$label, sep, fixed = TRUE),
                                   function(x) x[which(vars == var)])
    }
  }
  
  # Remove temporary label column from result_data
  result_data$stratum_label <- NULL
  
  # Attach strata_info as an attribute to the data for easy access
  attr(result_data, "strata_info") <- strata_info
  
  # Return results
  structure(
    list(
      data = result_data,
      strata_info = strata_info,
      vars = vars,
      sep = sep,
      min_n = min_n
    ),
    class = "maihda_strata"
  )
}

#' Print method for maihda_strata objects
#'
#' @param x A maihda_strata object
#' @param ... Additional arguments (not used)
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
