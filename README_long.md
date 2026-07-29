# dNdScv解析パイプライン詳細（version 1.2）

## 1. 目的

本パイプラインは、population filter後のANNOVAR注釈済み変異ファイルをがん種・疾患階層ごとに統合し、dNdScvによる遺伝子単位のpositive selection解析を実施するためのものである。

主な処理は次の5段階である。

1. サンプルシートから解析cohortごとのmanifestを作成する。
2. manifestに含まれる症例の変異をdNdScv入力形式へ統合する。
3. cohortごとにdNdScvを実行する。
4. 各cohortの`sel_cv`を統合し、positive selection候補遺伝子を抽出する。
5. `q_driver < 0.1`の遺伝子に該当する変異を、cohort・解析種別ごとに統合する。

通常は`CODE`単位、すなわちがん種ごとにarray jobを投入する。version 1.2では、`CODE`および`level1`～`level5`の全categoryについて通常解析とMSI-H除外解析を作成し、変異統合とdNdScvは両解析を別々のarray jobへ投入する。

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
│   └── merged_sample_sheet4_primary_sample_only.csv
├── output/
│   ├── manifest/
│   ├── manifest_lists/
│   ├── merged_mutation/
│   ├── merged_mutation_drivers/
│   ├── dndscv_lists/
│   └── dndscv/
├── log/
│   ├── manifest/
│   ├── merge_mutation/
│   ├── merge_drv/
│   └── dndscv/
├── src/
│   ├── R/
│   │   ├── make_manifest.R
│   │   ├── merge_mutation_for_dndscv.R
│   │   ├── merge_mutation.R
│   │   ├── merge_mutation_drivers.R
│   │   ├── run_dndscv.R
│   │   └── merge_sel_cv.R
│   ├── merge_mutation.sh
│   └── run_dndscv.sh
├── run_preparation.sh
├── submit_merge_mutation.sh
├── submit_merge_mutation_drivers.sh
└── submit_run_dndscv.sh
```

このパイプラインは、ANNOVARおよびpopulation filterが完了したファイルを入力とする。ANNOVAR stageとpopulation filter stageのログは、それぞれ上流パイプラインの`log/annovar/`、`log/filter_annovar/`を使用する。本パイプライン内ではこれらのstageを再実行しない。

---

## 4. 入出力の概要

| Stage | 主な入力 | 実施内容 | 主な出力 |
|---|---|---|---|
| 1. Manifest | `input/merged_sample_sheet4_primary_sample_only.csv` | cohortごとの症例一覧を作成。全categoryで通常解析とMSI-H除外解析の2種類を作成 | `output/manifest/*_manifest.txt` |
| 2. dNdScv用統合 | population filter後の`*.filtered.hg38_multianno.txt.gz`とmanifest | 症例IDを付加し、必要な5列へ変換してcohort単位で統合 | `output/merged_mutation/*_merged_mutation_for_dndscv.txt` |
| 2b. 確認用統合 | 同上 | 全ANNOVAR列を保持したまま染色体別に統合 | `output/merged_mutation/*_by_chr/*.txt` |
| 3. dNdScv | Stage 2の5列ファイル | `dndscv()`を実行 | `output/dndscv/sel_cv_*.csv`、`dndsCV_*.xlsx` |
| 4. 候補統合 | `output/dndscv/sel_cv_*.csv` | cohort列を追加し、指定条件でpositive selection候補を抽出 | `output/dndscv/sel_cv_merged_pos_selection.csv` |
| 5. Driver変異統合 | cohort manifest、`sel_cv_<cohort>.csv`、ANNOVAR変異 | `q_driver < 0.1`の遺伝子に該当する変異をcohort・解析種別ごとに統合 | cohortごとの通常版3ファイル、MSI-H除外版3ファイル |

---

## 5. Stage 1: Manifestの作成

### 実行

プロジェクトルートで実行する。

```bash
bash run_preparation.sh
```

`run_preparation.sh`はR 4.4.0をloadし、`src/R/make_manifest.R`を実行する。SGE resourceは`s_vmem=8G`、SGEとして投入する場合のログ出力先は`log/manifest/`である。`log/manifest/`は事前に作成しておく。

### 入力

```text
input/merged_sample_sheet4_primary_sample_only.csv
```

使用する主な列：

| 列 | 用途 |
|---|---|
| `tumor_sample_name` | 症例IDおよび変異ファイル名 |
| `CODE` | がん種および変異ファイルのサブディレクトリ |
| `parabricks` | 解析対象症例の確認 |
| `tumor_bam` | 解析対象症例の確認 |
| `cancer_description_en` | CODE manifestから除外する特殊群の判定 |
| `MSI_status` | 全categoryのMSI-H除外解析で、`MSI-H`症例を除外するために使用 |
| `level_1`～`level_5` | 疾患階層別manifestの作成 |

最初に、`parabricks`、`tumor_bam`、`tumor_sample_name`が欠損している行と、`tumor_sample_name`が空の行を除外する。

### Manifest作成条件

#### CODE

- `CODE`ごとに作成する。
- 最低症例数による除外は行わない。
- `cancer_description_en`が欠損している行を除外する。
- `cancer_description_en`が`reduction_team`または`patient_return`の行を除外する。
- 同じCODEについて、MSI statusによらない通常manifestと、末尾が`_MSIHexcluded`のmanifestを作成する。
- `_MSIHexcluded`では、前後空白を除いた`MSI_status`が明示的に`MSI-H`である症例だけを除外する。
- `MSI_status`が欠損、空欄、または`MSI-H`以外の症例は`_MSIHexcluded`にも残す。
- 再実行時には全categoryの既存`*_MSIHexcluded_manifest.txt`を一度削除して作り直し、現在のsample sheetと一致しない古いMSI-H除外manifestを残さない。

#### level 1～5

- `level_1`から`level_5`の各値ごとに作成する。
- 通常manifestと、明示的なMSI-H症例だけを除いた`_MSIHexcluded` manifestを作成する。
- uniqueな`tumor_sample_name`が30症例以上のgroupだけを出力する。この30例基準は通常版とMSI-H除外版で独立に判定する。
- 例えば通常版が35例でも、MSI-H除外後が29例なら通常manifestのみを出力する。

ファイル名に使用できない記号と空白は`_`へ置換される。

### 出力

```text
output/manifest/CODE_<group>_manifest.txt
output/manifest/CODE_<group>_MSIHexcluded_manifest.txt
output/manifest/level1_<group>_manifest.txt
output/manifest/level1_<group>_MSIHexcluded_manifest.txt
...
output/manifest/level5_<group>_manifest.txt
output/manifest/level5_<group>_MSIHexcluded_manifest.txt
```

拡張子は`.txt`だが、内容はヘッダー付きcomma-separated形式である。

```csv
tumor_sample_name,CODE
sample_001,AM
sample_002,AM
```

同じ`tumor_sample_name`と`CODE`の完全重複は除外し、`tumor_sample_name`順に並べる。

例として、`CODE_AB_manifest.txt`はMSI statusによらないAB全対象症例、`CODE_AB_MSIHexcluded_manifest.txt`はそのうち明示的なMSI-H症例を除外した症例集合である。level manifestも同じ2種類を作成する。MSI-H除外後に出力基準を満たさないgroupについては、空manifestを作成しない。

---

## 6. Stage 2: dNdScv用変異ファイルの統合

### 通常の実行

がん種ごとに処理する場合：

```bash
bash submit_merge_mutation.sh dndscv CODE
```

同じコマンドで各がん種の通常解析とMSI-H除外解析を処理するが、両者は別々のmanifest listと別々のSGE array jobとして投入される。

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
- 通常版とMSI-H除外版のmanifest listを`output/manifest_lists/`へ別々に出力する。
- manifest 1ファイルをSGE array jobの1 taskへ割り当てる。
- worker名は`src/merge_mutation.sh`である。
- R script名は`src/R/merge_mutation_for_dndscv.R`である。
- job名は通常版`prep_dnds_norm`、MSI-H除外版`prep_dnds_msi`、ログ出力先は`log/merge_mutation/`である。

`dndscv` modeでは`<cohort>_merged_mutation_for_dndscv.txt`が非空なら完了済みとしてskipする。`merge` modeではchr1～chr22、chrX、chrYの24ファイルがすべて非空の場合だけ完了済みとしてskipする。出力の一部だけが存在する場合や0 byte出力がある場合は、不完全出力としてqsub前に停止する。通常版とMSI-H除外版の両方を検証してから、いずれかのqsubを実行する。

manifest listの例：

```text
output/manifest_lists/dndscv_CODE_normal_YYYYMMDD_HHMMSS.txt
output/manifest_lists/dndscv_CODE_MSIHexcluded_YYYYMMDD_HHMMSS.txt
```

### 入力変異ファイル

各manifest行について、次のファイルを参照する。

```text
/home/nh1sy/analysis/wgs/3_mutation/pon_filtering/output/final_phased_local_dedup/<CODE>/<tumor_sample_name>.filtered.hg38_multianno.txt.gz
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
output/merged_mutation/CODE_AM_MSIHexcluded_merged_mutation_for_dndscv.txt
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
- job名は通常版`merge_mut_norm`、MSI-H除外版`merge_mut_msi`、ログ出力先は`log/merge_mutation/`である。
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

Stage 2と同様に、通常cohortと`_MSIHexcluded` cohortの両方が対象になる。両者は独立したanalysis IDを持ち、別々のinput list・別々のarray job・別々の完了済みskip判定を持つ。

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
- 通常版とMSI-H除外版のinput listを`output/dndscv_lists/`へ別々に出力する。
- cohortの統合変異ファイル1つをSGE array jobの1 taskへ割り当てる。
- worker名は`src/run_dndscv.sh`である。
- R script名は`src/R/run_dndscv.R`である。
- job名は通常版`run_dnds_norm`、MSI-H除外版`run_dnds_msi`、ログ出力先は`log/dndscv/`である。
- 1 task内のBLAS等のthread数を1へ制限する。

### 完了済みcohortのskip

各cohortについて、次の2ファイルが両方とも非空で存在する場合だけ完了済みと判定し、array input listから除外する。

```text
output/dndscv/sel_cv_<cohort>.csv
output/dndscv/dndsCV_<cohort>.xlsx
```

片方だけ存在する場合、またはいずれかが0 byteの場合は不完全出力としてsubmissionを停止する。全cohortが完了済みならqsubせず正常終了する。作成されるinput listには未完了cohortだけが含まれる。

input listの例：

```text
output/dndscv_lists/dndscv_CODE_normal_YYYYMMDD_HHMMSS.txt
output/dndscv_lists/dndscv_CODE_MSIHexcluded_YYYYMMDD_HHMMSS.txt
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
wmis_cv, wnon_cv, wspl_cv, qtrunc_cv, qallsubs_cv
```

indel modelが実行されたcohortでは、さらに`wind_cv`、`qind_cv`、`qglobal_cv`が存在する。`qglobal_cv`が存在する場合は`q_driver = qglobal_cv`、存在しない場合は`q_driver = qallsubs_cv`とし、採用した列名を`q_driver_source`へ記録する。必須列が欠けているファイル、または`qglobal_cv`とindel関連列の構成が矛盾するファイルがあれば停止する。

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

indel modelが実行されていないcohortでは`wind_cv`が存在しないため、substitutionの3列だけで判定する。

#### 2. いずれかの統計検定でFDR 10%未満

```r
qtrunc_cv < 0.1 |
qind_cv < 0.1 |
q_driver < 0.1
```

実際のfilterは次のとおりである。

```r
(wmis_cv > 1 | wnon_cv > 1 | wspl_cv > 1 | wind_cv > 1) &
  (qtrunc_cv < 0.1 | qind_cv < 0.1 | q_driver < 0.1)
```

これは`q_driver < 0.1`だけに限定するより広い候補集合である。truncating mutationまたはindelの検定だけで有意な遺伝子も残し、後の解釈・確認に利用できるようにしている。最終報告でglobal testだけを採用する場合は、この統合結果からさらに`q_driver < 0.1`で絞り込む。

結果は`cohort`、`q_driver`、`gene_name`の順に並べる。

### 出力

```text
output/dndscv/sel_cv_merged_pos_selection.csv
```

出力には、追加した`analysis_type`、`cohort`、`q_driver`、`q_driver_source`列と、元の`sel_cv`の全列を保持する。`analysis_type`は`normal`または`MSIHexcluded`であり、通常cohortと`_MSIHexcluded` cohortを明示的に区別する。

---

## 9b. Stage 5: q_driverで有意な遺伝子の変異統合

### 実行

```bash
# log/merge_drv/は事前に作成する
bash submit_merge_mutation_drivers.sh CODE
```

指定categoryのmanifestをarray input listへ登録し、共通worker `src/merge_mutation.sh`から`src/R/merge_mutation_drivers.R`をmanifest引数付きで実行する。通常版とMSI-H除外版は別々のtaskであり、初回実行では対象cohort数×2 taskとなる。job名は`merge_drv`、ログ出力先は`log/merge_drv/`、defaultの同時実行上限は20 taskである。

categoryは`CODE`、`level1`～`level5`から指定できる。省略時は全categoryを対象とする。各解析種別のmanifestが実際に作成されているcohortだけが対象になる。

### 入力

- `output/manifest/<cohort>_manifest.txt`
- `output/manifest/<cohort>_MSIHexcluded_manifest.txt`
- 対応する`output/dndscv/sel_cv_<cohort>.csv`または`sel_cv_<cohort>_MSIHexcluded.csv`
- `/home/nh1sy/analysis/wgs/3_mutation/pon_filtering/output/final_phased_local_dedup/<CODE>/<sample>.filtered.hg38_multianno.txt.gz`

### 遺伝子の選択と照合

1. `sel_cv_<cohort>.csv`に`qglobal_cv`があれば同列、なければ`qallsubs_cv`を`q_driver`として、`q_driver < 0.1`の行を選ぶ。
2. `CDKN2A.p14arf`および`CDKN2A.p16INK4a`は、ANNOVARとの照合時に`CDKN2A`へ変換する。
3. ANNOVARの`Gene.ensGene`を`;`または`,`で分割し、遺伝子名を完全一致で照合する。
4. `chr1`～`chr22`、`chrX`、`chrY`だけを出力する。
5. ANNOVARの全列を保持し、先頭に`sampleID`と`driver_gene`を追加する。

manifest内の変異ファイル、cohortに対応する`sel_cv`、または必須列が欠けている場合は停止し、自動skipしない。

### 出力

```text
output/merged_mutation_drivers/<cohort>_merged_mutation_drivers.txt
output/merged_mutation_drivers/<cohort>_significant_genes_q_driver_lt_0_1.tsv
output/merged_mutation_drivers/<cohort>_significant_genes_not_in_ensGene.tsv

output/merged_mutation_drivers/<cohort>_merged_mutation_drivers_MSIHexcluded.txt
output/merged_mutation_drivers/<cohort>_significant_genes_q_driver_lt_0_1_MSIHexcluded.tsv
output/merged_mutation_drivers/<cohort>_significant_genes_not_in_ensGene_MSIHexcluded.tsv
```

各taskは3ファイルを1つの完了bundleとして作成する。`*_merged_mutation_drivers*.txt`は該当ANNOVAR変異、`*_significant_genes_q_driver_lt_0_1*.tsv`は採用した有意遺伝子、`*_significant_genes_not_in_ensGene*.tsv`は当該cohortの`Gene.ensGene`に一度も現れなかった有意遺伝子である。有意遺伝子または未照合遺伝子が0件でもヘッダー付きファイルを作成する。

各taskについて3ファイルがすべて非空なら完了済みとしてarray input listから除外する。3ファイルの一部だけが存在する、またはいずれかが0 byteなら不完全bundleとして、いずれのqsubも行う前に停止する。対応する`sel_cv`がない、または0 byteの場合も同様に停止する。

---

## 10. 推奨実行順序

通常のCODE単位解析は以下の順で実行する。

```bash
# 事前に必要なlog directoryを用意する
mkdir -p log/manifest log/merge_mutation log/dndscv log/merge_drv

# 1. Manifest
bash run_preparation.sh

# 2. dNdScv用変異統合
bash submit_merge_mutation.sh dndscv CODE

# prep_dnds_normとprep_dnds_msiの完了を確認する

# 3. dNdScv
bash submit_run_dndscv.sh CODE

# run_dnds_normとrun_dnds_msiの完了を確認する

# 4. Positive selection候補の統合
Rscript src/R/merge_sel_cv.R

# 5. q_driver < 0.1の遺伝子に該当する変異を統合
# log/merge_drv/は事前に作成する
bash submit_merge_mutation_drivers.sh CODE
```

確認用に全ANNOVAR列も統合する場合だけ、別途以下を実行する。

```bash
mkdir -p log/merge_mutation
bash submit_merge_mutation.sh merge CODE
```

`submit_merge_mutation_drivers.sh`はlog directoryを作成しないため、`log/merge_drv/`を初回実行前に用意する。

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
| Driver抽出対象の`sel_cv`がない、または0 byte | Driver arrayをqsubする前に停止する |
| Driver抽出の3出力が一部だけ存在、またはいずれかが0 byte | 不完全bundleとしてDriver arrayをqsubする前に停止する |
| memory allocation error | 1 taskあたりの`s_vmem`を増やして再投入する |

dNdScvから隣接変異に関する警告が出た場合は、DNV/MNVがSNVへ分割されていないかを確認する。異なるsampleID間で同一変異が観察されたという警告は、既知のhotspotでは起こり得るが、件数が極端に多い場合は症例重複や上流のartifactを確認する。

---

## 12. 最終成果物

解析の主要成果物は次の4種類である。

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

4. cohort・解析種別ごとの有意driver遺伝子、未照合遺伝子、該当ANNOVAR変異

   ```text
   output/merged_mutation_drivers/<cohort>_merged_mutation_drivers[_MSIHexcluded].txt
   output/merged_mutation_drivers/<cohort>_significant_genes_q_driver_lt_0_1[_MSIHexcluded].tsv
   output/merged_mutation_drivers/<cohort>_significant_genes_not_in_ensGene[_MSIHexcluded].tsv
   ```

   上記の`[_MSIHexcluded]`は説明上の任意suffixであり、実ファイル名では空文字または`_MSIHexcluded`となる。

最終候補をglobal testで限定する場合は、3のファイルに対して`q_driver < 0.1`を追加適用する。

---

## 13. 参考資料

- dndscv package: [im3sanger/dndscv](https://github.com/im3sanger/dndscv)。入力は`sampleID, chr, pos, ref, mut`の5列で、indelはANNOVAR形式で指定するというpackage側の仕様を、本パイプラインの`Chr/Start/Ref/Alt`選択の根拠としている。
- Martincorena I, et al. *Universal Patterns of Selection in Cancer and Somatic Tissues*. Cell. 2017;171:1029-1041.e21. [doi:10.1016/j.cell.2017.09.042](https://doi.org/10.1016/j.cell.2017.09.042)
