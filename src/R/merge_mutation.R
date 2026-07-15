# 1 manifestを処理し、全変異情報を染色体別に統合するスクリプト
library(data.table)

options(readr.num_threads = 1)
Sys.setenv(VROOM_THREADS = 1)
setDTthreads(1)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 1L) {
  stop(
    "Usage: Rscript src/R/merge_mutation.R <manifest_file>\n",
    "Example: Rscript src/R/merge_mutation.R output/manifest/CODE_AB_manifest.txt"
  )
}

sheet_file <- normalizePath(args[[1]], mustWork = TRUE)
cat("Processing manifest:", sheet_file, "\n")

# 入出力設定
mutation_base_dir <- "/home/nh1sy/analysis/wgs/3_mutation/pon_filtering/output/final"
target_chrs <- paste0("chr", c(1:22, "X", "Y"))

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

# 不正行と完全重複行を除外
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

# CODE_AB_manifest.txt -> CODE_AB
manifest_id <- sub("_manifest\\.txt$", "", basename(sheet_file))

# 出力例:
# output/merged_mutation/CODE_AB_by_chr/CODE_AB_merged_mutation_chrY.txt
out_dir <- file.path(
  "output",
  "merged_mutation",
  paste0(manifest_id, "_by_chr")
)

# 途中で失敗した場合に不完全な本番出力を残さないため、一時ディレクトリへ出力
parent_out_dir <- dirname(out_dir)
dir.create(parent_out_dir, recursive = TRUE, showWarnings = FALSE)

tmp_out_dir <- paste0(out_dir, ".tmp_", Sys.getpid())
if (dir.exists(tmp_out_dir)) {
  unlink(tmp_out_dir, recursive = TRUE, force = TRUE)
}
dir.create(tmp_out_dir, recursive = TRUE, showWarnings = FALSE)

on.exit({
  if (dir.exists(tmp_out_dir)) {
    unlink(tmp_out_dir, recursive = TRUE, force = TRUE)
  }
}, add = TRUE)

out_prefix <- paste0(manifest_id, "_merged_mutation")
out_files <- setNames(
  file.path(tmp_out_dir, paste0(out_prefix, "_", target_chrs, ".txt")),
  target_chrs
)

total_rows_by_chr <- setNames(rep(0L, length(target_chrs)), target_chrs)
header_written <- FALSE

# 各症例ファイルを1回だけ読み込み、該当染色体の出力へ追記する
for (path in mutation_paths) {
  cat("Reading:", path, "\n")

  dt <- fread(
    cmd = paste("gzip -cd --", shQuote(path)),
    colClasses = "character",
    na.strings = c(".", "NA", ""),
    showProgress = FALSE
  )

  if (!"Chr" %in% names(dt)) {
    stop("Error: Column 'Chr' was not found in ", path)
  }

  # 変異が0件の染色体にもヘッダー付きファイルを作成
  if (!header_written) {
    empty_dt <- dt[0]
    for (chr_name in target_chrs) {
      fwrite(
        empty_dt,
        out_files[[chr_name]],
        sep = "\t",
        na = "",
        quote = FALSE
      )
    }
    header_written <- TRUE
  }

  # 1 / chr1 / X / chrX のいずれにも対応
  chr_body <- toupper(sub("^chr", "", dt$Chr, ignore.case = TRUE))
  dt[, Chr_norm := paste0("chr", chr_body)]

  # chr1-22, X, Yだけを対象とし、存在する染色体だけ書き込む
  dt_target <- dt[Chr_norm %chin% target_chrs]

  if (nrow(dt_target) > 0L) {
    chr_tables <- split(dt_target, by = "Chr_norm", keep.by = FALSE)

    for (chr_name in names(chr_tables)) {
      fwrite(
        chr_tables[[chr_name]],
        out_files[[chr_name]],
        sep = "\t",
        na = "",
        quote = FALSE,
        append = TRUE,
        col.names = FALSE
      )
      total_rows_by_chr[[chr_name]] <-
        total_rows_by_chr[[chr_name]] + nrow(chr_tables[[chr_name]])
    }
  }
}

if (!header_written) {
  stop("Error: No readable mutation data were found for ", sheet_file)
}

# 一時出力が完成してから本番ディレクトリへ置換
if (dir.exists(out_dir)) {
  unlink(out_dir, recursive = TRUE, force = TRUE)
}

if (!file.rename(tmp_out_dir, out_dir)) {
  stop("Error: Failed to move temporary output directory to ", out_dir)
}

cat("\nSaved chromosome-split results:\n")
for (chr_name in target_chrs) {
  final_file <- file.path(out_dir, basename(out_files[[chr_name]]))
  cat("  ", final_file, " (", total_rows_by_chr[[chr_name]], " rows)\n", sep = "")
}
