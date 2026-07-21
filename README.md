# dndscv

## prep
```
module use /usr/local/package/modulefiles/
module load R/4.4.0
R

> library(devtools); install_github("im3sanger/dndscv")
# Rsamtoolsのみうまくインストールできず手動でインストール
```

## run

### 1, Make manifest
```
bash run_preparation.sh
```
dndscvを行うsample groupをまとめたmanifest fileを作成する。


### 2, Merge mutation
Mutation fileをがん種ごとにまとめる。
```
# dndscv用のchr, pos, ref, altだけのファイルを作成する。
bash submit_merge_mutation.sh dndscv CODE

# 確認用に全てのcolumnを保持したファイルを作成する。重すぎるためにby_chrというディレクトリに染色体ごとに出力。
bash submit_merge_mutation.sh merge
```

### 3, dndscv
```
bash submit_run_dndscv.sh CODE
```

### 4, Merge positively selected genes
`output/dndscv/sel_cv_*.csv`を統合し、positive directionかつ`qtrunc_cv`、`qind_cv`、`qglobal_cv`のいずれかが0.1未満の遺伝子を抽出する。
```
Rscript src/R/merge_sel_cv.R
```

出力：
```
output/dndscv/sel_cv_merged_pos_selection.csv
```

### 5, Merge mutations in significant driver genes
cohortごとの`sel_cv_<cohort>.csv`から`qglobal_cv < 0.1`の遺伝子を選び、`Gene.ensGene`が一致するANNOVAR変異を全列保持して統合する。
`CDKN2A.p14arf`と`CDKN2A.p16INK4a`は`CDKN2A`として照合する。
```bash
# log/merge_drv/は事前に作成する
bash submit_merge_mutation_drivers.sh CODE
```

出力：
```text
output/merged_mutation_drivers/<cohort>_merged_mutation_drivers.txt
output/merged_mutation_drivers/<cohort>_significant_genes_qglobal_lt_0_1.tsv
output/merged_mutation_drivers/<cohort>_significant_genes_not_in_ensGene.tsv
```
