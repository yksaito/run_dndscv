# dNdScv解析パイプライン詳細

## 1. 目的

本パイプラインは、population filter後のANNOVAR注釈済み変異ファイルをがん種・疾患階層ごとに統合し、dNdScvによる遺伝子単位のpositive selection解析を実施するためのものである。

主な処理は次の4段階である。

1. サンプルシートから解析cohortごとのmanifestを作成する。
2. manifestに含まれる症例の変異をdNdScv入力形式へ統合する。
3. cohortごとにdNdScvを実行する。
4. 各cohortの`sel_cv`を統合し、positive selection候補遺伝子を抽出する。

通常は`CODE`単位、すなわちがん種ごとにarray jobを投入する。

---

## 2. 実行環境

### R

SHIROKANE上では以下を使用する。

```bash
module use /usr/local/package/modulefiles/
module load R/4.4.0
```

### 主なR package

- `tidyverse`
- `data.table`
- `dndscv`
- `openxlsx`

`dndscv`が未導入の場合の例：

```r
library(devtools)
install_github("im3sanger/dndscv")
```

R packageが別のR 4.4.xでbuildされていたという警告や、`translate`関数のimport置換に関する警告は、dNdScv本体が実行を継続している場合には通常、停止原因ではない。一方、`ベクトルを割り当てることができません`はメモリ不足を示す。

---

## 3. ディレクトリ構成

```text
dndscv/
├── input/
│   └── merged_sample_sheet3.csv
├── output/
│   ├── manifest/
│   ├── manifest_lists/
│   ├── merged_mutation/
│   ├── dndscv_lists/
│   └── dndscv/
├── log/
│   ├── prep/
│   ├── merge_dnd/
│   ├── merge_mut/
│   └── dndscv/
├── src/
│   ├── R/
│   │   ├── make_manifest.R
│   │   ├── merge_mutation_for_dndscv.R
│   │   ├── merge_mutation.R
│   │   ├── run_dndscv.R
│   │   └── merge_sel_cv.R
│   ├── merge_mutation.sh
│   └── run_dndscv.sh
├── run_preparation.sh
├── submit_merge_mutation.sh
└── submit_run_dndscv.sh
```

このパイプラインは、ANNOVARおよびpopulation filterが完了したファイルを入力とする。ANNOVAR stageとpopulation filter stageのログは、それぞれ上流パイプラインの`log/annovar/`、`log/filter_annovar/`を使用する。本パイプライン内ではこれらのstageを再実行しない。

---

## 4. 入出力の概要

| Stage | 主な入力 | 実施内容 | 主な出力 |
|---|---|---|---|
| 1. Manifest | `input/merged_sample_sheet3.csv` | cohortごとの症例一覧を作成 | `output/manifest/*_manifest.txt` |
| 2. dNdScv用統合 | population filter後の`*.filtered.hg38_multianno.txt.gz`とmanifest | 症例IDを付加し、必要な5列へ変換してcohort単位で統合 | `output/merged_mutation/*_merged_mutation_for_dndscv.txt` |
| 2b. 確認用統合 | 同上 | 全ANNOVAR列を保持したまま染色体別に統合 | `output/merged_mutation/*_by_chr/*.txt` |
| 3. dNdScv | Stage 2の5列ファイル | `dndscv()`を実行 | `output/dndscv/sel_cv_*.csv`、`dndsCV_*.xlsx` |
| 4. 候補統合 | `output/dndscv/sel_cv_*.csv` | cohort列を追加し、指定条件でpositive selection候補を抽出 | `output/dndscv/sel_cv_merged_pos_selection.csv` |

---

## 5. Stage 1: Manifestの作成

### 実行

プロジェクトルートで実行する。

```bash
bash run_preparation.sh
```

`run_preparation.sh`はR 4.4.0をloadし、`src/R/make_manifest.R`を実行する。SGE resourceは`s_vmem=8G`、ログ出力先は`log/prep/`である。`log/prep/`は事前に作成しておく。

### 入力

```text
input/merged_sample_sheet3.csv
```

使用する主な列：

| 列 | 用途 |
|---|---|
| `tumor_sample_name` | 症例IDおよび変異ファイル名 |
| `CODE` | がん種および変異ファイルのサブディレクトリ |
| `parabricks` | 解析対象症例の確認 |
| `tumor_bam` | 解析対象症例の確認 |
| `cancer_description_en` | CODE manifestから除外する特殊群の判定 |
| `level_1`～`level_5` | 疾患階層別manifestの作成 |

最初に、`parabricks`、`tumor_bam`、`tumor_sample_name`が欠損している行と、`tumor_sample_name`が空の行を除外する。

### Manifest作成条件

#### CODE

- `CODE`ごとに作成する。
- 最低症例数による除外は行わない。
- `cancer_description_en`が欠損している行を除外する。
- `cancer_description_en`が`reduction_team`または`patient_return`の行を除外する。

#### level 1～5

- `level_1`から`level_5`の各値ごとに作成する。
- uniqueな`tumor_sample_name`が30症例以上のgroupだけを出力する。

ファイル名に使用できない記号と空白は`_`へ置換される。

### 出力

```text
output/manifest/CODE_<group>_manifest.txt
output/manifest/level1_<group>_manifest.txt
...
output/manifest/level5_<group>_manifest.txt
```

拡張子は`.txt`だが、内容はヘッダー付きcomma-separated形式である。

```csv
tumor_sample_name,CODE
sample_001,AM
sample_002,AM
```

同じ`tumor_sample_name`と`CODE`の完全重複は除外し、`tumor_sample_name`順に並べる。

---

## 6. Stage 2: dNdScv用変異ファイルの統合

### 通常の実行

がん種ごとに処理する場合：

```bash
bash submit_merge_mutation.sh dndscv CODE
```

複数階層を指定する場合：

```bash
bash submit_merge_mutation.sh dndscv CODE level1 level2
```

categoryを省略した場合は、`CODE level1 level2 level3 level4 level5`をすべて対象とする。

同時実行task数は環境変数`MAX_CONCURRENT`で変更できる。現行実装のdefaultは20である。

```bash
MAX_CONCURRENT=5 bash submit_merge_mutation.sh dndscv CODE
```

### Array jobの構成

- `submit_merge_mutation.sh`が対象manifestを列挙する。
- manifest listを`output/manifest_lists/`へ出力する。
- manifest 1ファイルをSGE array jobの1 taskへ割り当てる。
- worker名は`src/merge_mutation.sh`である。
- R script名は`src/R/merge_mutation_for_dndscv.R`である。
- job名は`merge_dnd`、ログ出力先は`log/merge_dnd/`である。

manifest listの例：

```text
output/manifest_lists/dndscv_CODE_YYYYMMDD_HHMMSS.txt
```

### 入力変異ファイル

各manifest行について、次のファイルを参照する。

```text
/home/nh1sy/analysis/wgs/3_mutation/pon_filtering/output/annovar_filtered/<CODE>/<tumor_sample_name>.filtered.hg38_multianno.txt.gz
```

必要なANNOVAR列は以下である。

| ANNOVAR列 | dNdScv入力列 | 内容 |
|---|---|---|
| manifestの`tumor_sample_name` | `sampleID` | 症例ID |
| `Chr` | `chr` | 染色体 |
| `Start` | `pos` | ANNOVAR形式の変異開始座標 |
| `Ref` | `ref` | reference allele |
| `Alt` | `mut` | alternate allele |

`Otherinfo4/5/7/8`ではなく、ANNOVARが正規化した`Chr/Start/Ref/Alt`を使用する。dNdScvのindel入力はANNOVAR形式を想定しているため、VCFの`POS/REF/ALT`表現へ戻さない。これにより、挿入の`ref = "-"`、欠失の`mut = "-"`など、ANNOVAR形式のindel表現をそのまま渡す。

MNVまたは複合置換がANNOVARファイルで1 eventとして記録されている場合も、`Ref`と`Alt`を分割せず1 eventのまま出力する。

### 実施内容

1. manifestから`tumor_sample_name`と`CODE`を読み込む。
2. 不正な行と完全重複行を除外する。
3. 対応する全症例の変異ファイルが存在することを確認する。
4. 各症例ファイルから`Chr/Start/Ref/Alt`のみを読む。
5. `sampleID/chr/pos/ref/mut`へ変換する。
6. 全症例を一度にメモリへ保持せず、症例ごとに一時ファイルへ追記する。
7. 全症例の処理完了後に最終ファイルへrenameする。

manifest内のファイルが1つでも存在しない場合は停止する。欠損症例を自動skipしない。

### 出力

```text
output/merged_mutation/<manifest_id>_merged_mutation_for_dndscv.txt
```

例：

```text
output/merged_mutation/CODE_AM_merged_mutation_for_dndscv.txt
```

出力はヘッダー付きtab-separated形式である。

```text
sampleID	chr	pos	ref	mut
sample_001	1	123456	A	G
sample_001	2	234567	-	AT
sample_002	3	345678	CT	-
```

---

## 7. Stage 2b: 全ANNOVAR列を保持した確認用統合

この処理はdNdScv本体の入力作成には不要である。注釈列を含む元データをcohort単位で確認したい場合だけ実行する。

### 実行

```bash
bash submit_merge_mutation.sh merge CODE
```

### 実施内容

- worker名は`src/merge_mutation.sh`のまま使用する。
- R scriptは`src/R/merge_mutation.R`を使用する。
- job名は`merge_mut`、ログ出力先は`log/merge_mut/`である。
- ANNOVARの全列を保持する。
- `Chr`が`1`、`chr1`、`X`、`chrX`などのいずれでも判定できるよう、一時的な正規化列を作る。
- `chr1`～`chr22`、`chrX`、`chrY`だけを対象とする。
- 出力サイズを分割するため、染色体ごとに出力する。
- 変異が0件の染色体にもヘッダー付きファイルを作成する。

### 出力

```text
output/merged_mutation/<manifest_id>_by_chr/
├── <manifest_id>_merged_mutation_chr1.txt
├── ...
├── <manifest_id>_merged_mutation_chr22.txt
├── <manifest_id>_merged_mutation_chrX.txt
└── <manifest_id>_merged_mutation_chrY.txt
```

---

## 8. Stage 3: dNdScvの実行

### 通常の実行

がん種ごとに実行する。

```bash
bash submit_run_dndscv.sh CODE
```

categoryを省略した場合は、`CODE level1 level2 level3 level4 level5`をすべて対象とする。

同時実行task数のdefaultは2である。

```bash
MAX_CONCURRENT=1 bash submit_run_dndscv.sh CODE
```

reference databaseのdefaultは`hg38`である。別のdNdScv referenceを使う場合は`REFDB`で指定する。

```bash
REFDB=RefCDS_human_GRCh38.p12_dNdScv.0.1.0.rda \
  bash submit_run_dndscv.sh CODE
```

### Array jobの構成

- `submit_run_dndscv.sh`が指定categoryの統合変異ファイルを列挙する。
- input listを`output/dndscv_lists/`へ出力する。
- cohortの統合変異ファイル1つをSGE array jobの1 taskへ割り当てる。
- worker名は`src/run_dndscv.sh`である。
- R script名は`src/R/run_dndscv.R`である。
- job名は`run_dndscv`、ログ出力先は`log/dndscv/`である。
- 1 task内のBLAS等のthread数を1へ制限する。

input listの例：

```text
output/dndscv_lists/dndscv_CODE_YYYYMMDD_HHMMSS.txt
```

### 入力

```text
output/merged_mutation/<cohort>_merged_mutation_for_dndscv.txt
```

必須列は次の5列である。

| 列 | Rでの型 | 内容 |
|---|---|---|
| `sampleID` | character | 症例ID |
| `chr` | character | 染色体 |
| `pos` | double | 座標 |
| `ref` | character | reference allele |
| `mut` | character | mutant allele |

### 実施内容

1. 入力ファイル名から解析IDを作成する。例：`CODE_AM_merged_mutation_for_dndscv.txt`から`CODE_AM`。
2. `chr`先頭の`chr`を除去する。
3. 必須値が欠損または空の行を除外する。
4. 染色体`1`～`22`、`X`、`Y`だけを残す。
5. 次の呼び出しでdNdScvを実行する。

```r
dndsout <- dndscv(
  mutations,
  refdb = refdb
)
```

`onesided`などは明示的に変更せず、導入されているdndscv packageのdefaultを使用する。

### 出力

#### `sel_cv` CSV

```text
output/dndscv/sel_cv_<cohort>.csv
```

例：

```text
output/dndscv/sel_cv_CODE_AM.csv
```

`dndsout$sel_cv`をCSVとして保存する。positive selection候補の統合にはこのファイルを使用する。

#### 全出力Excel

```text
output/dndscv/dndsCV_<cohort>.xlsx
```

`dndscv()`が返すlist全体をsheet別に保存する。`sel_cv`以外の結果やannotated mutationsを確認するときに使用する。

### メモリに関する注意

現行の`submit_run_dndscv.sh`は`qsub`時に`s_vmem`を明示していない。そのため、clusterのdefault memoryが小さい環境では、dNdScvのenvironment load中または解析中に次のようなエラーで停止することがある。

```text
エラー: サイズ 110 Kb のベクトルを割り当てることができません
```

このメッセージは110 Kbのデータ自体が大きいという意味ではなく、taskがmemory上限へ到達し、新しいallocationができなかったことを示す。大規模cohortでは、`submit_run_dndscv.sh`の`qsub` resourceに、実測に応じた十分な`s_vmem`、例えば32Gを指定して再投入する。

```bash
-l s_vmem=32G \
```

同時実行数`MAX_CONCURRENT`を下げることはnode全体の混雑緩和には有効だが、1 taskあたりのmemory上限そのものは増加させない。

---

## 9. Stage 4: Positive selection候補の統合

### 実行

Stage 3の全array taskが完了した後、プロジェクトルートで実行する。

```bash
Rscript src/R/merge_sel_cv.R
```

### 入力

```text
output/dndscv/sel_cv_*.csv
```

最終出力である`sel_cv_merged_pos_selection.csv`自身は入力から除外する。

各ファイルには少なくとも次の列が必要である。

```text
wmis_cv, wnon_cv, wspl_cv, wind_cv,
qtrunc_cv, qind_cv, qglobal_cv
```

必要列が1つでも欠けているファイルがあれば停止する。

### cohort列

ファイル名から`sel_cv_`と`.csv`を除いた文字列を`cohort`として先頭列へ追加する。

```text
sel_cv_CODE_AM.csv  ->  cohort = CODE_AM
```

### 抽出条件

次の両方を満たす行を残す。

#### 1. いずれかの変異classでpositive direction

```r
wmis_cv > 1 |
wnon_cv > 1 |
wspl_cv > 1 |
wind_cv > 1
```

#### 2. いずれかの統計検定でFDR 10%未満

```r
qtrunc_cv < 0.1 |
qind_cv < 0.1 |
qglobal_cv < 0.1
```

実際のfilterは次のとおりである。

```r
(wmis_cv > 1 | wnon_cv > 1 | wspl_cv > 1 | wind_cv > 1) &
  (qtrunc_cv < 0.1 | qind_cv < 0.1 | qglobal_cv < 0.1)
```

これは`qglobal_cv < 0.1`だけに限定するより広い候補集合である。truncating mutationまたはindelの検定だけで有意な遺伝子も残し、後の解釈・確認に利用できるようにしている。最終報告でglobal testだけを採用する場合は、この統合結果からさらに`qglobal_cv < 0.1`で絞り込む。

結果は`cohort`、`qglobal_cv`、`gene_name`の順に並べる。

### 出力

```text
output/dndscv/sel_cv_merged_pos_selection.csv
```

出力には、追加した`cohort`列と、元の`sel_cv`の全列を保持する。

---

## 10. 推奨実行順序

通常のCODE単位解析は以下の順で実行する。

```bash
# 事前に必要なlog directoryを用意する
mkdir -p log/prep log/merge_dnd log/dndscv

# 1. Manifest
bash run_preparation.sh

# 2. dNdScv用変異統合
bash submit_merge_mutation.sh dndscv CODE

# merge_dnd array jobの完了を確認する

# 3. dNdScv
bash submit_run_dndscv.sh CODE

# run_dndscv array jobの完了を確認する

# 4. Positive selection候補の統合
Rscript src/R/merge_sel_cv.R
```

確認用に全ANNOVAR列も統合する場合だけ、別途以下を実行する。

```bash
mkdir -p log/merge_mut
bash submit_merge_mutation.sh merge CODE
```

`submit_merge_mutation.sh`と`submit_run_dndscv.sh`は、listおよび当該jobのlog directoryを現行実装内で作成する。一方、`run_preparation.sh`は`log/prep/`を作成しないため、初回実行前に用意する。

---

## 11. 主な停止条件と確認事項

| 状況 | 動作・確認事項 |
|---|---|
| manifest directoryがない | merge submissionを停止する |
| 指定categoryに一致するmanifestがない | merge submissionを停止する |
| manifestの必須列がない | 当該taskを停止する |
| manifest内の変異ファイルが1つでもない | 当該taskを停止する。自動skipしない |
| ANNOVARの`Chr/Start/Ref/Alt`がない | dNdScv用merge taskを停止する |
| dNdScv入力の必須5列がない | dNdScv taskを停止する |
| 入力validation後に変異が0件 | dNdScv taskを停止する |
| `dndsout$sel_cv`がない | 出力処理を停止する |
| `sel_cv`に統合用必須列がない | Stage 4を停止する |
| memory allocation error | 1 taskあたりの`s_vmem`を増やして再投入する |

dNdScvから隣接変異に関する警告が出た場合は、DNV/MNVがSNVへ分割されていないかを確認する。異なるsampleID間で同一変異が観察されたという警告は、既知のhotspotでは起こり得るが、件数が極端に多い場合は症例重複や上流のartifactを確認する。

---

## 12. 最終成果物

解析の主要成果物は次の3種類である。

1. cohort別のgene-level selection結果

   ```text
   output/dndscv/sel_cv_<cohort>.csv
   ```

2. cohort別のdNdScv全出力

   ```text
   output/dndscv/dndsCV_<cohort>.xlsx
   ```

3. 全cohortのpositive selection候補統合結果

   ```text
   output/dndscv/sel_cv_merged_pos_selection.csv
   ```

最終候補をglobal testで限定する場合は、3のファイルに対して`qglobal_cv < 0.1`を追加適用する。

---

## 13. 参考資料

- dndscv package: [im3sanger/dndscv](https://github.com/im3sanger/dndscv)。入力は`sampleID, chr, pos, ref, mut`の5列で、indelはANNOVAR形式で指定するというpackage側の仕様を、本パイプラインの`Chr/Start/Ref/Alt`選択の根拠としている。
- Martincorena I, et al. *Universal Patterns of Selection in Cancer and Somatic Tissues*. Cell. 2017;171:1029-1041.e21. [doi:10.1016/j.cell.2017.09.042](https://doi.org/10.1016/j.cell.2017.09.042)
