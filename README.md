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
bash submit_merge_mutation.sh dndscv

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
