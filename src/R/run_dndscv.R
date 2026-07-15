# dNdScvを1つの統合変異ファイルに対して実行するスクリプト
suppressPackageStartupMessages({
  library(tidyverse)
  library(dndscv)
  library(openxlsx)
})

options(readr.num_threads = 1)
Sys.setenv(VROOM_THREADS = 1)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1L || length(args) > 2L) {
  stop(
    "Usage: Rscript src/R/run_dndscv.R <merged_mutation_for_dndscv.txt> [refdb]\n",
    "Example: Rscript src/R/run_dndscv.R ",
    "output/merged_mutation/CODE_AB_merged_mutation_for_dndscv.txt hg38"
  )
}

input_file <- normalizePath(args[[1]], mustWork = TRUE)
refdb <- if (length(args) == 2L) args[[2]] else "hg38"

# 入力ファイル名から解析名を取得する
# 例: CODE_AB_merged_mutation_for_dndscv.txt -> CODE_AB
input_basename <- basename(input_file)
analysis_id <- str_remove(
  input_basename,
  "_merged_mutation_for_dndscv\\.txt$"
)

# 想定したsuffixを持たない場合も、拡張子を除いた名前で実行可能にする
if (identical(analysis_id, input_basename)) {
  analysis_id <- tools::file_path_sans_ext(input_basename)
}

if (is.na(analysis_id) || analysis_id == "") {
  stop("Error: Could not determine analysis ID from input file: ", input_file)
}

cat("========================================\n")
cat("Starting dNdScv analysis\n")
cat("Input file :", input_file, "\n")
cat("Analysis ID:", analysis_id, "\n")
cat("Reference  :", refdb, "\n")
cat("========================================\n")

mutations <- read_tsv(
  input_file,
  col_types = cols(
    sampleID = col_character(),
    chr = col_character(),
    pos = col_double(),
    ref = col_character(),
    mut = col_character()
  ),
  progress = FALSE
) %>%
  as.data.frame()

required_cols <- c("sampleID", "chr", "pos", "ref", "mut")
missing_cols <- setdiff(required_cols, colnames(mutations))

if (length(missing_cols) > 0L) {
  stop(
    "Error: Missing required column(s): ",
    paste(missing_cols, collapse = ", ")
  )
}

# chr1形式と1形式のどちらでも処理できるように統一する
mutations <- mutations %>%
  transmute(
    sampleID = sampleID,
    chr = str_remove(chr, "^chr"),
    pos = pos,
    ref = ref,
    mut = mut
  ) %>%
  filter(
    !is.na(sampleID), sampleID != "",
    !is.na(chr), chr != "",
    !is.na(pos),
    !is.na(ref), ref != "",
    !is.na(mut), mut != ""
  ) %>%
  filter(chr %in% c(as.character(1:22), "X", "Y"))

if (nrow(mutations) == 0L) {
  stop("Error: No valid mutations remained after input validation.")
}

cat("Number of mutations:", nrow(mutations), "\n")
cat("Number of samples  :", n_distinct(mutations$sampleID), "\n")

# dNdScvの実行
dndsout <- dndscv(
  mutations,
  refdb = refdb
)

out_dir <- file.path("output", "dndscv")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

sel_cv_file <- file.path(
  out_dir,
  paste0("sel_cv_", analysis_id, ".csv")
)

all_output_file <- file.path(
  out_dir,
  paste0("dndsCV_", analysis_id, ".xlsx")
)

if (is.null(dndsout$sel_cv)) {
  stop("Error: dNdScv output does not contain 'sel_cv'.")
}

write_excel_csv(dndsout$sel_cv, sel_cv_file)
openxlsx::write.xlsx(dndsout, all_output_file, overwrite = TRUE)

cat("Saved sel_cv output to:", sel_cv_file, "\n")
cat("Saved all outputs to :", all_output_file, "\n")
cat("dNdScv analysis completed successfully.\n")

# 警告1：Mutations observed in contiguous sites...
# 訳：同一サンプル内で隣接するサイトに変異が観察されました。
# 意味: ゲノム上で隣り合った2つの塩基が両方とも変異している（例：chr1:1000 A>C, chr1:1001 T>G）箇所があるという警告です。
# 解釈: 喫煙に伴う特徴的な変異シグネチャや特定の変異原曝露においては、ジヌクレオチド変異（DNV: Dinucleotide Variant）や複雑な置換（Complex substitution）が頻発することがあります。これらを個別の1塩基置換（SNV）として扱うと、同義/非同義の判定やアミノ酸変化の予測が狂い、dN/dS比の計算にバイアスがかかる可能性があります。
# 対策: そのままでも計算は完遂しますが、厳密なドライバー遺伝子探索を行う場合は、バリアントコールの段階でこれらをDNV/MNVとして適切に結合しておくことが推奨されます。

# 警告2：Same mutations observed in different sampleIDs...
# 訳：異なるサンプルIDで全く同じ変異が観察されました。
# 意味: 複数の患者（サンプル）間で、染色体・座標・置換パターンが全く同じ変異が存在するという警告です。
# 解釈: KRAS や EGFR などの有名なドライバー遺伝子のホットスポット変異であれば、複数の患者で共通して検出されるのは自然なことです。しかし、もしこの共通変異が数千〜数万レベルで大量に存在する場合は注意が必要です。
# 対策: ホットスポット由来であれば無視して問題ありません。ただし、WGSデータのマージ処理におけるサンプルの重複、PCR増幅アーティファクト、あるいはアライメントエラーの可能性もゼロではないため、共通している変異が既知のドライバー遺伝子に集中しているかどうかを軽く確認（サニティチェック）しておくことをお勧めします。