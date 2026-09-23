# Shared validation; this file has no package dependencies or side effects.

parse_replicates <- function(value, label, minimum = 1L) {
  number <- suppressWarnings(as.numeric(value))
  if (length(number) != 1L || !is.finite(number) ||
      number != floor(number) || number < minimum ||
      number > .Machine$integer.max) {
    stop(label, " must be a whole number between ", minimum, " and ",
         .Machine$integer.max, ".", call. = FALSE)
  }
  as.integer(number)
}

validate_unique_ids <- function(data, label) {
  if (!"ID" %in% names(data)) {
    stop(label, " has no ID column.", call. = FALSE)
  }
  ids <- as.character(data$ID)
  if (anyNA(ids) || any(!nzchar(trimws(ids)))) {
    stop(label, " contains missing or blank IDs.", call. = FALSE)
  }
  if (any(ids != trimws(ids))) {
    stop(label, " contains IDs with leading or trailing spaces.", call. = FALSE)
  }
  duplicates <- unique(ids[duplicated(ids)])
  if (length(duplicates)) {
    stop(label, " contains duplicate IDs: ",
         paste(duplicates, collapse = ", "), call. = FALSE)
  }
  invisible(TRUE)
}

