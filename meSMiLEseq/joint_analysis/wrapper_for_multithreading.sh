#!/bin/bash

# Number of cores per script (adjust as needed, total = 6 for 01_.py and 6 * 4 for 04_.py, so 30 cores)

CORES_SCRIPT1=6  
CORES_SCRIPT2=6

# List of experiments
EXPERIMENTS=$(printf "SmSAG%02d\n" {3..3})

# Paths to your scripts
SCRIPT1="01_kmer_analysis.py"
SCRIPT2="04_pyProBound_analysis.py"

# Barcode groups (split BC1–BC12 into 3 groups of 4)
BARCODE_GROUPS=(
    "BC1 BC2 BC3 BC4"
    "BC5 BC6 BC7 BC8"
    "BC9 BC10 BC11 BC12"
)

# Define function for script 1: run for one experiment + one barcode group
run_kmer_split() {
    exp=$1
    bc_group="$2"
    echo "[KMER] Starting $exp with barcodes: $bc_group"
    python3 "$SCRIPT1" -sms "$exp" -bc $bc_group
    echo "[KMER] Finished $exp with barcodes: $bc_group"
}

# Define function for script 2 (still runs one job per experiment)
run_probound_split() {
    exp=$1
    bc_group="$2"
    echo "[PROBOUND] Starting $exp with barcodes: $bc_group"
    python3 -u "$SCRIPT2" -sms "$exp" -bc $bc_group
    echo "[PROBOUND] Finished $exp with barcodes: $bc_group"
}

# Export functions and variables
export -f run_kmer_split
export -f run_probound_split
export SCRIPT1 SCRIPT2

# from here new#######################
# Prepare input for script 1 (each line = "SmSAGxx barcode_chunk")
kmer_jobs=()
for exp in $EXPERIMENTS; do
    for bc_group in "${BARCODE_GROUPS[@]}"; do
        kmer_jobs+=("$exp:::${bc_group}")
    done
done

# Generate jobs for Script 2 (experiment + barcode group)
probound_jobs=()
for exp in $EXPERIMENTS; do
    for bc_group in "${BARCODE_GROUPS[@]}"; do
        probound_jobs+=("$exp:::${bc_group}")
    done
done

# Run Script 1: parallelized by experiment and barcode group
printf "%s\n" "${kmer_jobs[@]}" | \
    parallel --line-buffer -j "$CORES_SCRIPT1" --colsep ':::' run_kmer_split {1} {2} &

# Run Script2 in parallel
printf "%s\n" "${probound_jobs[@]}" | \
    parallel --line-buffer -j "$CORES_SCRIPT2" --colsep ':::' run_probound_split {1} {2}

# Wait for both background jobs
wait

echo "All experiments completed for both scripts."




