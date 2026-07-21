library(tidyverse)

# -----------------------------------------------------------------------------
# Settings
# -----------------------------------------------------------------------------
input_file <- "input/merged_sample_sheet4_primary_sample_only.csv"
output_dir <- "output/manifest"
minimum_cases <- 30

# -----------------------------------------------------------------------------
# Load sample sheet
# -----------------------------------------------------------------------------
samplesheet <- read_csv(input_file, show_col_types = FALSE) %>%
  filter(
    !is.na(parabricks),
    !is.na(tumor_bam),
    !is.na(tumor_sample_name),
    tumor_sample_name != ""
  )

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Make group names safe for use as file names
sanitize_filename <- function(x) {
  x %>%
    as.character() %>%
    str_trim() %>%
    str_replace_all("[/\\\\:*?\"<>|]", "_") %>%
    str_replace_all("\\s+", "_") %>%
    str_replace_all("_+", "_") %>%
    str_remove_all("^_|_$")
}

# Write one manifest per category
write_manifests <- function(data,
                            group_col,
                            prefix,
                            min_cases = NULL) {

  group_sym <- sym(group_col)

  group_counts <- data %>%
    filter(
      !is.na(!!group_sym),
      !!group_sym != ""
    ) %>%
    distinct(!!group_sym, tumor_sample_name) %>%
    count(!!group_sym, name = "n")

  if (!is.null(min_cases)) {
    group_counts <- group_counts %>%
      filter(n >= min_cases)
  }

  if (nrow(group_counts) == 0) {
    message("No eligible categories found for ", group_col)
    return(invisible(NULL))
  }

  walk(group_counts[[group_col]], function(group_value) {
    safe_group_value <- sanitize_filename(group_value)

    output_file <- file.path(
      output_dir,
      paste0(prefix, "_", safe_group_value, "_manifest.txt")
    )

    manifest <- data %>%
      filter(!!group_sym == group_value) %>%
      distinct(tumor_sample_name, CODE) %>%
      arrange(tumor_sample_name)

    write_csv(manifest, output_file)

    message(
      "Written: ", output_file,
      " (", nrow(manifest), " samples)"
    )
  })

  invisible(group_counts)
}

# -----------------------------------------------------------------------------
# CODE manifests
# Exclude reduction_team and patient_return
# -----------------------------------------------------------------------------
code_samplesheet <- samplesheet %>%
  filter(
    !is.na(cancer_description_en),
    !cancer_description_en %in% c("reduction_team", "patient_return")
  )

write_manifests(
  data = code_samplesheet,
  group_col = "CODE",
  prefix = "CODE"
)

# -----------------------------------------------------------------------------
# level_1 to level_5 manifests
# Only categories with at least 30 unique samples
# -----------------------------------------------------------------------------
for (level_col in paste0("level_", 1:5)) {
  write_manifests(
    data = samplesheet,
    group_col = level_col,
    prefix = str_replace(level_col, "level_", "level"),
    min_cases = minimum_cases
  )
}
