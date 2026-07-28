#!/bin/bash

set -euo pipefail

PIPELINE_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
TEST_ROOT=$(mktemp -d)
trap 'rm -rf -- "$TEST_ROOT"' EXIT

mkdir -p \
  "$TEST_ROOT/bin" \
  "$TEST_ROOT/output/manifest" \
  "$TEST_ROOT/output/merged_mutation" \
  "$TEST_ROOT/output/dndscv" \
  "$TEST_ROOT/log/merge_mutation" \
  "$TEST_ROOT/log/dndscv" \
  "$TEST_ROOT/log/merge_drv"

printf '%s\n' \
  '#!/bin/bash' \
  'printf '\''%s\n'\'' "$*" >> "$QSUB_LOG"' \
  > "$TEST_ROOT/bin/qsub"
chmod +x "$TEST_ROOT/bin/qsub"
export PATH="$TEST_ROOT/bin:$PATH"
export QSUB_LOG="$TEST_ROOT/qsub.log"

for manifest in \
  CODE_A_manifest.txt \
  CODE_B_manifest.txt \
  CODE_A_MSIHexcluded_manifest.txt \
  CODE_B_MSIHexcluded_manifest.txt; do
  printf 'tumor_sample_name,CODE\n' > "$TEST_ROOT/output/manifest/$manifest"
done

assert_two_jobs() {
  local first_name="$1"
  local second_name="$2"
  [[ $(wc -l < "$QSUB_LOG") -eq 2 ]]
  grep -q -- "-N $first_name" "$QSUB_LOG"
  grep -q -- "-N $second_name" "$QSUB_LOG"
}

# submit_merge_mutation.sh dndscv: one completed cohort in each analysis set.
printf 'done\n' > "$TEST_ROOT/output/merged_mutation/CODE_A_merged_mutation_for_dndscv.txt"
printf 'done\n' > "$TEST_ROOT/output/merged_mutation/CODE_B_MSIHexcluded_merged_mutation_for_dndscv.txt"
(
  cd "$TEST_ROOT"
  bash "$PIPELINE_DIR/submit_merge_mutation.sh" dndscv CODE
)
assert_two_jobs prep_dnds_norm prep_dnds_msi
normal_list=$(find "$TEST_ROOT/output/manifest_lists" -name 'dndscv_CODE_normal_*.txt')
msi_list=$(find "$TEST_ROOT/output/manifest_lists" -name 'dndscv_CODE_MSIHexcluded_*.txt')
grep -q 'CODE_B_manifest.txt' "$normal_list"
! grep -q 'CODE_A_manifest.txt' "$normal_list"
grep -q 'CODE_A_MSIHexcluded_manifest.txt' "$msi_list"
! grep -q 'CODE_B_MSIHexcluded_manifest.txt' "$msi_list"

# submit_merge_mutation.sh merge: complete means all 24 chromosome files.
: > "$QSUB_LOG"
for cohort in CODE_A CODE_B_MSIHexcluded; do
  out_dir="$TEST_ROOT/output/merged_mutation/${cohort}_by_chr"
  mkdir -p "$out_dir"
  for chr in {1..22} X Y; do
    printf 'header\n' > "$out_dir/${cohort}_merged_mutation_chr${chr}.txt"
  done
done
(
  cd "$TEST_ROOT"
  bash "$PIPELINE_DIR/submit_merge_mutation.sh" merge CODE
)
assert_two_jobs merge_mut_norm merge_mut_msi

# submit_run_dndscv.sh: independent output skip and separate arrays.
: > "$QSUB_LOG"
printf 'done\n' > "$TEST_ROOT/output/merged_mutation/CODE_B_merged_mutation_for_dndscv.txt"
printf 'done\n' > "$TEST_ROOT/output/merged_mutation/CODE_A_MSIHexcluded_merged_mutation_for_dndscv.txt"
printf 'done\n' > "$TEST_ROOT/output/dndscv/sel_cv_CODE_A.csv"
printf 'done\n' > "$TEST_ROOT/output/dndscv/dndsCV_CODE_A.xlsx"
printf 'done\n' > "$TEST_ROOT/output/dndscv/sel_cv_CODE_B_MSIHexcluded.csv"
printf 'done\n' > "$TEST_ROOT/output/dndscv/dndsCV_CODE_B_MSIHexcluded.xlsx"
(
  cd "$TEST_ROOT"
  bash "$PIPELINE_DIR/submit_run_dndscv.sh" CODE
)
assert_two_jobs run_dnds_norm run_dnds_msi
normal_list=$(find "$TEST_ROOT/output/dndscv_lists" -name 'dndscv_CODE_normal_*.txt')
msi_list=$(find "$TEST_ROOT/output/dndscv_lists" -name 'dndscv_CODE_MSIHexcluded_*.txt')
grep -q 'CODE_B_merged_mutation_for_dndscv.txt' "$normal_list"
grep -q 'CODE_A_MSIHexcluded_merged_mutation_for_dndscv.txt' "$msi_list"

# All dNdScv outputs complete: no qsub.
: > "$QSUB_LOG"
printf 'done\n' > "$TEST_ROOT/output/dndscv/sel_cv_CODE_B.csv"
printf 'done\n' > "$TEST_ROOT/output/dndscv/dndsCV_CODE_B.xlsx"
printf 'done\n' > "$TEST_ROOT/output/dndscv/sel_cv_CODE_A_MSIHexcluded.csv"
printf 'done\n' > "$TEST_ROOT/output/dndscv/dndsCV_CODE_A_MSIHexcluded.xlsx"
(
  cd "$TEST_ROOT"
  bash "$PIPELINE_DIR/submit_run_dndscv.sh" CODE
)
[[ ! -s "$QSUB_LOG" ]]

# Partial dNdScv output: fail before either qsub.
rm -f "$TEST_ROOT/output/dndscv/dndsCV_CODE_B.xlsx"
: > "$QSUB_LOG"
if (
  cd "$TEST_ROOT"
  bash "$PIPELINE_DIR/submit_run_dndscv.sh" CODE
); then
  echo "Expected partial-output test to fail" >&2
  exit 1
fi
[[ ! -s "$QSUB_LOG" ]]

# Driver merge wrapper: four manifests become four independent array tasks.
printf 'done\n' > "$TEST_ROOT/output/dndscv/dndsCV_CODE_B.xlsx"
for cohort in CODE_A CODE_B CODE_A_MSIHexcluded CODE_B_MSIHexcluded; do
  printf 'done\n' > "$TEST_ROOT/output/dndscv/sel_cv_${cohort}.csv"
done

: > "$QSUB_LOG"
(
  cd "$TEST_ROOT"
  bash "$PIPELINE_DIR/submit_merge_mutation_drivers.sh" CODE
)
[[ $(wc -l < "$QSUB_LOG") -eq 1 ]]
grep -q -- '-N merge_drv' "$QSUB_LOG"
grep -q -- '-t 1-4' "$QSUB_LOG"
driver_list=$(find "$TEST_ROOT/output/manifest_lists" -name 'drivers_CODE_*.txt')
[[ $(wc -l < "$driver_list") -eq 4 ]]

# One complete three-file bundle is skipped; the other three tasks remain.
mkdir -p "$TEST_ROOT/output/merged_mutation_drivers"
for output in \
  CODE_A_merged_mutation_drivers.txt \
  CODE_A_significant_genes_q_driver_lt_0_1.tsv \
  CODE_A_significant_genes_not_in_ensGene.tsv; do
  printf 'header\n' > "$TEST_ROOT/output/merged_mutation_drivers/$output"
done
: > "$QSUB_LOG"
(
  cd "$TEST_ROOT"
  bash "$PIPELINE_DIR/submit_merge_mutation_drivers.sh" CODE
)
grep -q -- '-t 1-3' "$QSUB_LOG"

# A partial bundle stops submission before qsub.
printf 'header\n' > \
  "$TEST_ROOT/output/merged_mutation_drivers/CODE_B_merged_mutation_drivers_MSIHexcluded.txt"
: > "$QSUB_LOG"
if (
  cd "$TEST_ROOT"
  bash "$PIPELINE_DIR/submit_merge_mutation_drivers.sh" CODE
); then
  echo "Expected partial driver-output bundle test to fail" >&2
  exit 1
fi
[[ ! -s "$QSUB_LOG" ]]

# All four bundles complete: no qsub.
rm -f \
  "$TEST_ROOT/output/merged_mutation_drivers/CODE_B_merged_mutation_drivers_MSIHexcluded.txt"
for cohort in CODE_B; do
  for output in \
    "${cohort}_merged_mutation_drivers.txt" \
    "${cohort}_significant_genes_q_driver_lt_0_1.tsv" \
    "${cohort}_significant_genes_not_in_ensGene.tsv"; do
    printf 'header\n' > "$TEST_ROOT/output/merged_mutation_drivers/$output"
  done
done
for cohort in CODE_A CODE_B; do
  for output in \
    "${cohort}_merged_mutation_drivers_MSIHexcluded.txt" \
    "${cohort}_significant_genes_q_driver_lt_0_1_MSIHexcluded.tsv" \
    "${cohort}_significant_genes_not_in_ensGene_MSIHexcluded.tsv"; do
    printf 'header\n' > "$TEST_ROOT/output/merged_mutation_drivers/$output"
  done
done
: > "$QSUB_LOG"
(
  cd "$TEST_ROOT"
  bash "$PIPELINE_DIR/submit_merge_mutation_drivers.sh" CODE
)
[[ ! -s "$QSUB_LOG" ]]

echo "All submission-wrapper tests passed."
