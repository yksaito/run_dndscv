# q_driver < 0.1の遺伝子に該当するANNOVAR変異をcohort・解析種別ごとに統合する
suppressPackageStartupMessages(library(data.table))

script_version <- "1.2"
cat("merge_mutation_drivers.R version", script_version, "\n")

options(readr.num_threads = 1)
Sys.setenv(VROOM_THREADS = 1)
setDTthreads(1)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 1L) {
  stop(
    "Usage: Rscript src/R/merge_mutation_drivers.R <manifest_file>\n",
    "Example: Rscript src/R/merge_mutation_drivers.R ",
    "output/manifest/CODE_AB_manifest.txt"
  )
}

sheet_file <- normalizePath(args[[1]], mustWork = TRUE)
manifest_id <- sub("_manifest\\.txt$", "", basename(sheet_file))
is_msi_h_excluded <- grepl("_MSIHexcluded$", manifest_id)
cohort <- sub("_MSIHexcluded$", "", manifest_id)
analysis_type <- if (is_msi_h_excluded) "MSIHexcluded" else "normal"

cat("Processing manifest:", sheet_file, "\n")
cat("Cohort:", cohort, "\n")
cat("Analysis type:", analysis_type, "\n")

mutation_base_dir <-
  "/home/nh1sy/analysis/wgs/3_mutation/pon_filtering/output/final_phased_local_dedup"
dndscv_file <- file.path(
  "output", "dndscv", paste0("sel_cv_", manifest_id, ".csv")
)
out_dir <- file.path("output", "merged_mutation_drivers")
target_chrs <- paste0("chr", c(1:22, "X", "Y"))

if (!file.exists(dndscv_file) || file.info(dndscv_file)$size == 0L) {
  stop("Error: dNdScv result was not found or is empty: ", dndscv_file)
}

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

dndscv_result <- fread(
  dndscv_file,
  na.strings = c("", "NA"),
  showProgress = FALSE
)

required_dndscv_cols <- c("gene_name", "qallsubs_cv")
missing_dndscv_cols <- setdiff(required_dndscv_cols, names(dndscv_result))
if (length(missing_dndscv_cols) > 0L) {
  stop(
    "Error: Missing column(s) in ", dndscv_file, ": ",
    paste(missing_dndscv_cols, collapse = ", ")
  )
}

if ("qglobal_cv" %in% names(dndscv_result)) {
  missing_indel_cols <- setdiff(c("wind_cv", "qind_cv"), names(dndscv_result))
  if (length(missing_indel_cols) > 0L) {
    stop(
      "Error: qglobal_cv is present but indel column(s) are missing in ",
      dndscv_file, ": ", paste(missing_indel_cols, collapse = ", ")
    )
  }
  q_driver_column <- "qglobal_cv"
} else {
  unexpected_indel_cols <- intersect(c("wind_cv", "qind_cv"), names(dndscv_result))
  if (length(unexpected_indel_cols) > 0L) {
    stop(
      "Error: qglobal_cv is missing but indel column(s) are present in ",
      dndscv_file, ": ", paste(unexpected_indel_cols, collapse = ", ")
    )
  }
  q_driver_column <- "qallsubs_cv"
}

q_driver_numeric <- suppressWarnings(
  as.numeric(dndscv_result[[q_driver_column]])
)
invalid_q_driver <-
  !is.na(dndscv_result[[q_driver_column]]) & is.na(q_driver_numeric)
if (any(invalid_q_driver)) {
  stop(
    "Error: Non-numeric ", q_driver_column,
    " value(s) were found in ", dndscv_file
  )
}
dndscv_result[, q_driver := q_driver_numeric]
dndscv_result[, q_driver_source := q_driver_column]

normalize_dndscv_gene <- function(x) {
  x <- trimws(as.character(x))
  x[x %chin% c("CDKN2A.p14arf", "CDKN2A.p16INK4a")] <- "CDKN2A"
  x
}

significant_genes <- dndscv_result[
  !is.na(gene_name) & gene_name != "" &
    !is.na(q_driver) & q_driver < 0.1,
  .(
    dndscv_gene = as.character(gene_name),
    ensGene_gene = normalize_dndscv_gene(gene_name),
    q_driver,
    q_driver_source
  )
]
significant_genes <- unique(significant_genes)
setorder(significant_genes, q_driver, dndscv_gene)

driver_genes <- unique(significant_genes$ensGene_gene)

cat("Significant dNdScv rows:", nrow(significant_genes), "\n")
cat("Driver q-value source:", q_driver_column, "\n")
cat("Unique ensGene matching names:", length(driver_genes), "\n")

extract_gene_tokens <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""

  tokens <- unlist(strsplit(x, "[;,]", perl = TRUE), use.names = FALSE)
  tokens <- trimws(tokens)
  tokens <- sub("\\(dist=.*\\)$", "", tokens)
  unique(tokens[tokens != "" & tokens != "."])
}

match_driver_genes <- function(x, drivers) {
  if (length(drivers) == 0L) {
    return(rep("", length(x)))
  }

  vapply(
    strsplit(replace(as.character(x), is.na(x), ""), "[;,]", perl = TRUE),
    function(tokens) {
      tokens <- trimws(tokens)
      tokens <- sub("\\(dist=.*\\)$", "", tokens)
      matched <- sort(unique(tokens[tokens %chin% drivers]))
      paste(matched, collapse = ";")
    },
    character(1)
  )
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

suffix <- if (is_msi_h_excluded) "_MSIHexcluded" else ""
mutation_out <- file.path(
  out_dir,
  paste0(cohort, "_merged_mutation_drivers", suffix, ".txt")
)
gene_out <- file.path(
  out_dir,
  paste0(cohort, "_significant_genes_q_driver_lt_0_1", suffix, ".tsv")
)
missing_gene_out <- file.path(
  out_dir,
  paste0(cohort, "_significant_genes_not_in_ensGene", suffix, ".tsv")
)

mutation_tmp <- paste0(mutation_out, ".tmp_", Sys.getpid())
gene_tmp <- paste0(gene_out, ".tmp_", Sys.getpid())
missing_gene_tmp <- paste0(missing_gene_out, ".tmp_", Sys.getpid())
tmp_files <- c(mutation_tmp, gene_tmp, missing_gene_tmp)

unlink(tmp_files[file.exists(tmp_files)], force = TRUE)
on.exit(unlink(tmp_files[file.exists(tmp_files)], force = TRUE), add = TRUE)

seen_ensgenes <- character()
expected_annovar_cols <- NULL
header_written <- FALSE
total_input_rows <- 0L
total_driver_rows <- 0L

for (i in seq_along(mutation_paths)) {
  path <- mutation_paths[[i]]
  sample_id <- samplesheet$tumor_sample_name[[i]]
  cat("Reading:", path, "\n")

  dt <- fread(
    cmd = paste("gzip -cd --", shQuote(path)),
    colClasses = "character",
    na.strings = c(".", "NA", ""),
    showProgress = FALSE
  )

  required_annovar_cols <- c("Chr", "Gene.ensGene")
  missing_annovar_cols <- setdiff(required_annovar_cols, names(dt))
  if (length(missing_annovar_cols) > 0L) {
    stop(
      "Error: Missing column(s) in ", path, ": ",
      paste(missing_annovar_cols, collapse = ", ")
    )
  }

  if (is.null(expected_annovar_cols)) {
    expected_annovar_cols <- copy(names(dt))
  } else if (!identical(names(dt), expected_annovar_cols)) {
    stop("Error: ANNOVAR column structure differs in ", path)
  }

  total_input_rows <- total_input_rows + nrow(dt)

  chr_body <- toupper(sub("^chr", "", dt$Chr, ignore.case = TRUE))
  canonical <- paste0("chr", chr_body) %chin% target_chrs

  if (any(canonical)) {
    seen_ensgenes <- union(
      seen_ensgenes,
      extract_gene_tokens(dt[["Gene.ensGene"]][canonical])
    )
  }

  matched_genes <- match_driver_genes(dt[["Gene.ensGene"]], driver_genes)
  keep <- canonical & matched_genes != ""

  dt[, sampleID := sample_id]
  dt[, driver_gene := matched_genes]
  setcolorder(dt, c("sampleID", "driver_gene", expected_annovar_cols))

  if (!header_written) {
    fwrite(
      dt[0],
      mutation_tmp,
      sep = "\t",
      na = "",
      quote = FALSE
    )
    header_written <- TRUE
  }

  if (any(keep)) {
    fwrite(
      dt[keep],
      mutation_tmp,
      sep = "\t",
      na = "",
      quote = FALSE,
      append = TRUE,
      col.names = FALSE
    )
    total_driver_rows <- total_driver_rows + sum(keep)
  }
}

if (!header_written) {
  stop("Error: No readable mutation data were found for ", sheet_file)
}

significant_genes[, present_in_Gene.ensGene := ensGene_gene %chin% seen_ensgenes]
not_in_ensgene <- significant_genes[present_in_Gene.ensGene == FALSE]

fwrite(
  significant_genes,
  gene_tmp,
  sep = "\t",
  na = "",
  quote = FALSE
)
fwrite(
  not_in_ensgene,
  missing_gene_tmp,
  sep = "\t",
  na = "",
  quote = FALSE
)

final_files <- c(mutation_out, gene_out, missing_gene_out)
for (path in final_files[file.exists(final_files)]) {
  unlink(path, force = TRUE)
}

if (!file.rename(mutation_tmp, mutation_out)) {
  stop("Error: Failed to move temporary output file to ", mutation_out)
}
if (!file.rename(gene_tmp, gene_out)) {
  stop("Error: Failed to move temporary output file to ", gene_out)
}
if (!file.rename(missing_gene_tmp, missing_gene_out)) {
  stop("Error: Failed to move temporary output file to ", missing_gene_out)
}

cat("\nSaved result to:", mutation_out, "\n")
cat("Saved significant genes to:", gene_out, "\n")
cat("Saved genes absent from ensGene to:", missing_gene_out, "\n")
cat("Input mutation rows:", total_input_rows, "\n")
cat("Driver mutation rows:", total_driver_rows, "\n")
cat("Significant genes absent from ensGene:", nrow(not_in_ensgene), "\n")
