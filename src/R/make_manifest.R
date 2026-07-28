library(tidyverse)

script_version <- "1.2"
message("make_manifest.R version ", script_version)

# -----------------------------------------------------------------------------
# Settings
# -----------------------------------------------------------------------------
input_file <- "input/merged_sample_sheet4_primary_sample_only.csv"
output_dir <- "output/manifest"
minimum_cases <- 30

# -----------------------------------------------------------------------------
# Load sample sheet
# -----------------------------------------------------------------------------
samplesheet <- read_csv(input_file, show_col_types = FALSE)

required_cols <- c(
  "tumor_sample_name", "CODE", "parabricks", "tumor_bam",
  "cancer_description_en", "MSI_status", paste0("level_", 1:5)
)
missing_cols <- setdiff(required_cols, names(samplesheet))
if (length(missing_cols) > 0L) {
  stop(
    "Error: Missing required column(s) in ", input_file, ": ",
    paste(missing_cols, collapse = ", ")
  )
}

samplesheet <- samplesheet %>%
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
                            min_cases = NULL,
                            analysis_suffix = "") {

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
      paste0(
        prefix, "_", safe_group_value, analysis_suffix, "_manifest.txt"
      )
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

# 全categoryについて、明示的にMSI-Hと判定された症例だけを除外した
# 解析用manifestも作成する。MSI_statusが欠損または空欄の症例は残す。
samplesheet_msi_h_excluded <- samplesheet %>%
  filter(is.na(MSI_status) | str_trim(MSI_status) != "MSI-H")

# 再実行時に、現在のsample sheetでは対象外または30例未満になったgroupの
# 古いmanifestを誤って後段へ投入しないよう、全MSI-H除外manifestを作り直す。
old_msi_h_excluded_manifests <- list.files(
  output_dir,
  pattern = "^(CODE|level[1-5])_.*_MSIHexcluded_manifest\\.txt$",
  full.names = TRUE
)
if (length(old_msi_h_excluded_manifests) > 0L) {
  unlink(old_msi_h_excluded_manifests)
}

write_manifests(
  data = samplesheet_msi_h_excluded %>%
    filter(
      !is.na(cancer_description_en),
      !cancer_description_en %in% c("reduction_team", "patient_return")
    ),
  group_col = "CODE",
  prefix = "CODE",
  analysis_suffix = "_MSIHexcluded"
)

# -----------------------------------------------------------------------------
# level_1 to level_5 manifests
# Only categories with at least 30 unique samples
# -----------------------------------------------------------------------------
for (level_col in paste0("level_", 1:5)) {
  level_prefix <- str_replace(level_col, "level_", "level")

  write_manifests(
    data = samplesheet,
    group_col = level_col,
    prefix = level_prefix,
    min_cases = minimum_cases
  )

  write_manifests(
    data = samplesheet_msi_h_excluded,
    group_col = level_col,
    prefix = level_prefix,
    min_cases = minimum_cases,
    analysis_suffix = "_MSIHexcluded"
  )
}
