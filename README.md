# dndscv version 1.2

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
`input/merged_sample_sheet4_primary_sample_only.csv`から、dndscvを行うsample groupをまとめたmanifest fileを作成する。

`CODE`および`level1`～`level5`の各groupについて、以下の2解析を作成する。

```text
CODE_AB_manifest.txt
CODE_AB_MSIHexcluded_manifest.txt
level2_Skin_manifest.txt
level2_Skin_MSIHexcluded_manifest.txt
```

通常manifestは`MSI_status`によらず全対象症例を含む。`_MSIHexcluded` manifestは`MSI_status`が明示的に`MSI-H`の症例だけを除外し、欠損・空欄・その他のstatusは残す。level groupでは通常版、MSI-H除外版のそれぞれに30例以上の基準を適用する。


### 2, Merge mutation
Mutation fileをがん種ごとにまとめる。
```
# dndscv用のchr, pos, ref, altだけのファイルを作成する。
bash submit_merge_mutation.sh dndscv CODE

# 確認用に全てのcolumnを保持したファイルを作成する。重すぎるためにby_chrというディレクトリに染色体ごとに出力。
bash submit_merge_mutation.sh merge
```

指定した全categoryについて、通常cohortと`_MSIHexcluded` cohortを別々のarray jobへ投入する。既存の完全な出力は各arrayのinput listから除外する。不完全出力はerrorとする。

### 3, dndscv
```
bash submit_run_dndscv.sh CODE
```

指定した全categoryの通常cohortと`_MSIHexcluded` cohortを別々のarray jobで解析する。完了済みcohortは各解析で独立にskipする。

### 4, Merge positively selected genes
`output/dndscv/sel_cv_*.csv`を統合し、positive directionかつ`qtrunc_cv`、`qind_cv`、`q_driver`のいずれかが0.1未満の遺伝子を抽出する。`q_driver`は`qglobal_cv`が存在するcohortでは同列、存在しないcohortでは`qallsubs_cv`を使用し、`q_driver_source`に採用列を記録する。`analysis_type`列で通常版とMSI-H除外版を識別する。
```
Rscript src/R/merge_sel_cv.R
```

出力：
```
output/dndscv/sel_cv_merged_pos_selection.csv
```

### 5, Merge mutations in significant driver genes
各cohortの`sel_cv_<cohort>.csv`から`q_driver < 0.1`の遺伝子を選び、`Gene.ensGene`が一致するANNOVAR変異を全列保持して、cohort・解析種別ごとに統合する。通常版とMSI-H除外版は別々のarray taskとして実行する。
`CDKN2A.p14arf`と`CDKN2A.p16INK4a`は`CDKN2A`として照合する。
```bash
# log/merge_drv/は事前に作成する
bash submit_merge_mutation_drivers.sh CODE
```

1 cohortあたりの出力：
```text
output/merged_mutation_drivers/<cohort>_merged_mutation_drivers.txt
output/merged_mutation_drivers/<cohort>_significant_genes_q_driver_lt_0_1.tsv
output/merged_mutation_drivers/<cohort>_significant_genes_not_in_ensGene.tsv

output/merged_mutation_drivers/<cohort>_merged_mutation_drivers_MSIHexcluded.txt
output/merged_mutation_drivers/<cohort>_significant_genes_q_driver_lt_0_1_MSIHexcluded.tsv
output/merged_mutation_drivers/<cohort>_significant_genes_not_in_ensGene_MSIHexcluded.tsv
```

初回実行では、通常版とMSI-H除外版を合わせて対象cohort数×2 taskとなる。各taskの3出力がすべて非空なら完了済みとしてskipする。3出力の一部だけが存在する、またはいずれかが0 byteなら不完全出力としてqsub前に停止する。
