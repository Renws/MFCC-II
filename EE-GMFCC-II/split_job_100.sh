#!/usr/bin/env bash
# Split jobs in rungjf.raw into up to 100 Slurm submission scripts
# Assumption: each line in rungjf.raw is one Gaussian job command

set -euo pipefail

RAW="rungjf.raw"
TEMPLATE="run_gauss.sh"
MAX_SPLITS=100

TMP_PREFIX="rungjf"
OUT_PREFIX="run_gauss"

if [[ ! -s "$RAW" ]]; then
    echo "ERROR: $RAW does not exist or is empty"
    exit 1
fi

if [[ ! -f "$TEMPLATE" ]]; then
    echo "ERROR: $TEMPLATE does not exist"
    exit 1
fi

job_counts=$(wc -l < "$RAW" | tr -d ' ')

# If there are fewer than 100 jobs, only generate non-empty scripts
split_counts=$(( job_counts < MAX_SPLITS ? job_counts : MAX_SPLITS ))

echo "Total jobs: $job_counts"
echo "Generate Slurm scripts: $split_counts"

# Remove old files to avoid mixing with previous results
rm -f ${TMP_PREFIX}[0-9][0-9] ${OUT_PREFIX}[0-9][0-9].sh

# Distribute lines evenly into split_counts temporary files
awk -v n="$split_counts" -v total="$job_counts" -v prefix="$TMP_PREFIX" '
{
    idx = int((NR - 1) * n / total)
    file = sprintf("%s%02d", prefix, idx)
    print $0 >> file
}
' "$RAW"

# Generate Slurm submission scripts
for ((i=0; i<split_counts; i++)); do
    idx=$(printf "%02d" "$i")
    chunk="${TMP_PREFIX}${idx}"
    out="${OUT_PREFIX}${idx}.sh"

    cp "$TEMPLATE" "$out"

    # Change Slurm job name, e.g. gauss00, gauss01, ..., gauss99
    sed -i "s/^#SBATCH[[:space:]]\+-J[[:space:]]\+.*/#SBATCH -J gauss${idx}/" "$out"

    printf "\n" >> "$out"
    cat "$chunk" >> "$out"

    chmod +x "$out"
    rm -f "$chunk"
done

echo "Done."
echo "Example single submission: sbatch ${OUT_PREFIX}00.sh"
echo "Batch submission command:"
echo "for f in ${OUT_PREFIX}[0-9][0-9].sh; do sbatch \"\$f\"; done"
