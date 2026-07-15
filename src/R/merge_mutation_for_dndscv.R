# 1 manifestを処理し、dNdScv入力用変異ファイルを作成するスクリプト
library(data.table)

options(readr.num_threads = 1)
Sys.setenv(VROOM_THREADS = 1)
setDTthreads(1)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 1L) {
  stop(
    "Usage: Rscript src/R/merge_mutation_for_dndscv.R <manifest_file>\n",
    "Example: Rscript src/R/merge_mutation_for_dndscv.R output/manifest/CODE_AB_manifest.txt"
  )
}

sheet_file <- normalizePath(args[[1]], mustWork = TRUE)
cat("Processing manifest:", sheet_file, "\n")

mutation_base_dir <- "/home/nh1sy/analysis/wgs/3_mutation/pon_filtering/output/final"
out_base_dir <- file.path("output", "merged_mutation")
dir.create(out_base_dir, recursive = TRUE, showWarnings = FALSE)

samplesheet <- fread(
  sheet_file,
  colClasses = "character",
  na.strings = c("", "NA"),
  showProgress = FALSE
)

required_manifest_cols <- c("tumor_sample_name", "CODE")
missing_manifest_cols <- setdiff(required_manifest_cols, names(samplesheet))
if (length(missing_manifest_cols) > 0L) {
  stop(
    "Error: Missing column(s) in manifest: ",
    paste(missing_manifest_cols, collapse = ", ")
  )
}

samplesheet <- unique(
  samplesheet[
    !is.na(tumor_sample_name) & tumor_sample_name != "" &
      !is.na(CODE) & CODE != "",
    .(tumor_sample_name, CODE)
  ]
)

if (nrow(samplesheet) == 0L) {
  stop("Error: No valid samples were found in ", sheet_file)
}

mutation_paths <- file.path(
  mutation_base_dir,
  samplesheet$CODE,
  paste0(samplesheet$tumor_sample_name, ".filtered.hg38_multianno.txt.gz")
)

missing_paths <- mutation_paths[!file.exists(mutation_paths)]
if (length(missing_paths) > 0L) {
  stop(
    "Error: The following mutation file(s) were not found:\n",
    paste(missing_paths, collapse = "\n")
  )
}

manifest_id <- sub("_manifest\\.txt$", "", basename(sheet_file))
out_file <- file.path(
  out_base_dir,
  paste0(manifest_id, "_merged_mutation_for_dndscv.txt")
)
tmp_file <- paste0(out_file, ".tmp_", Sys.getpid())

if (file.exists(tmp_file)) {
  unlink(tmp_file, force = TRUE)
}

on.exit({
  if (file.exists(tmp_file)) {
    unlink(tmp_file, force = TRUE)
  }
}, add = TRUE)

required_mutation_cols <- c("Chr", "Start", "Ref", "Alt")
header_written <- FALSE
total_rows <- 0L

# 全症例をメモリ上で結合せず、1症例ずつ読み込んで追記する
for (i in seq_along(mutation_paths)) {
  path <- mutation_paths[[i]]
  cat("Reading:", path, "\n")

  dt <- fread(
    cmd = paste("gzip -cd --", shQuote(path)),
    select = required_mutation_cols,
    colClasses = "character",
    na.strings = c(".", "NA", ""),
    showProgress = FALSE
  )

  missing_mutation_cols <- setdiff(required_mutation_cols, names(dt))
  if (length(missing_mutation_cols) > 0L) {
    stop(
      "Error: Missing column(s) in ", path, ": ",
      paste(missing_mutation_cols, collapse = ", ")
    )
  }

  setnames(
    dt,
    old = required_mutation_cols,
    new = c("chr", "pos", "ref", "mut")
  )
  dt[, sampleID := samplesheet$tumor_sample_name[[i]]]
  setcolorder(dt, c("sampleID", "chr", "pos", "ref", "mut"))

  fwrite(
    dt,
    tmp_file,
    sep = "\t",
    na = "",
    quote = FALSE,
    append = header_written,
    col.names = !header_written
  )

  header_written <- TRUE
  total_rows <- total_rows + nrow(dt)
}

if (!header_written) {
  stop("Error: No readable mutation data were found for ", sheet_file)
}

if (file.exists(out_file)) {
  unlink(out_file, force = TRUE)
}

if (!file.rename(tmp_file, out_file)) {
  stop("Error: Failed to move temporary output file to ", out_file)
}

cat("Saved result to:", out_file, "\n")
cat("Total rows:", total_rows, "\n")
