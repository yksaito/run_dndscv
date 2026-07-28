library(tidyverse)

input_dir <- file.path("output", "dndscv")
output_file <- file.path(input_dir, "sel_cv_merged_pos_selection.csv")

input_files <- list.files(
  input_dir,
  pattern = "^sel_cv_.*\\.csv$",
  full.names = TRUE
)
input_files <- setdiff(input_files, output_file)

if (length(input_files) == 0L) {
  stop("Error: No sel_cv_*.csv files were found in ", input_dir)
}

required_cols <- c(
  "gene_name", "wmis_cv", "wnon_cv", "wspl_cv",
  "qtrunc_cv", "qallsubs_cv"
)

merged <- map_dfr(input_files, function(input_file) {
  cohort <- basename(input_file) %>%
    str_remove("^sel_cv_") %>%
    str_remove("\\.csv$")

  result <- read_csv(input_file, show_col_types = FALSE)

  missing_cols <- setdiff(required_cols, names(result))
  if (length(missing_cols) > 0L) {
    stop(
      "Error: Missing column(s) in ", input_file, ": ",
      paste(missing_cols, collapse = ", ")
    )
  }

  if ("qglobal_cv" %in% names(result)) {
    q_driver_column <- "qglobal_cv"

    missing_indel_cols <- setdiff(c("wind_cv", "qind_cv"), names(result))
    if (length(missing_indel_cols) > 0L) {
      stop(
        "Error: qglobal_cv is present but indel column(s) are missing in ",
        input_file, ": ", paste(missing_indel_cols, collapse = ", ")
      )
    }
  } else {
    unexpected_indel_cols <- intersect(c("wind_cv", "qind_cv"), names(result))
    if (length(unexpected_indel_cols) > 0L) {
      stop(
        "Error: qglobal_cv is missing but indel column(s) are present in ",
        input_file, ": ", paste(unexpected_indel_cols, collapse = ", ")
      )
    }

    q_driver_column <- "qallsubs_cv"
  }

  q_driver <- suppressWarnings(as.numeric(result[[q_driver_column]]))
  invalid_q_driver <-
    !is.na(result[[q_driver_column]]) & is.na(q_driver)
  if (any(invalid_q_driver)) {
    stop(
      "Error: Non-numeric ", q_driver_column,
      " value(s) were found in ", input_file
    )
  }

  positive_direction <-
    result$wmis_cv > 1 |
    result$wnon_cv > 1 |
    result$wspl_cv > 1

  if ("wind_cv" %in% names(result)) {
    positive_direction <- positive_direction | result$wind_cv > 1
  }
  positive_direction[is.na(positive_direction)] <- FALSE

  significant_q <- result$qtrunc_cv < 0.1 | q_driver < 0.1
  if ("qind_cv" %in% names(result)) {
    significant_q <- significant_q | result$qind_cv < 0.1
  }
  significant_q[is.na(significant_q)] <- FALSE

  result %>%
    mutate(
      analysis_type = if_else(
        str_ends(cohort, "_MSIHexcluded"),
        "MSIHexcluded",
        "normal"
      ),
      cohort = cohort,
      q_driver = .env$q_driver,
      q_driver_source = .env$q_driver_column,
      .before = 1
    ) %>%
    filter(.env$positive_direction & .env$significant_q)
}) %>%
  arrange(cohort, q_driver, gene_name)

write_excel_csv(merged, output_file)

cat("Input files:", length(input_files), "\n")
cat("Selected rows:", nrow(merged), "\n")
cat("Saved result to:", output_file, "\n")
