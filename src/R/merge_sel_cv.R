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
  "wmis_cv", "wnon_cv", "wspl_cv", "wind_cv",
  "qtrunc_cv", "qind_cv", "qglobal_cv"
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

  result %>%
    mutate(cohort = cohort, .before = 1) %>%
    filter(
      (wmis_cv > 1 | wnon_cv > 1 | wspl_cv > 1 | wind_cv > 1) &
        (qtrunc_cv < 0.1 | qind_cv < 0.1 | qglobal_cv < 0.1)
    )
}) %>%
  arrange(cohort, qglobal_cv, gene_name)

write_excel_csv(merged, output_file)

cat("Input files:", length(input_files), "\n")
cat("Selected rows:", nrow(merged), "\n")
cat("Saved result to:", output_file, "\n")
