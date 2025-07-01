#!/bin/bash

# Number of cores per script (adjust as needed, total = 3 for 01_.py and 12 for 04_.py, so 15 cores)
CORES_PER_SCRIPT=3

# List of experiments
EXPERIMENTS=$(printf "SmSAG%02d\n" {1..2})

# Paths to your scripts
SCRIPT1="01_kmer_analysis.py"
SCRIPT2="04_pyProBound_analysis.py"

# Define function for script 1
run_kmer() {
    exp=$1
    echo "[KMER] Starting $exp"
    python3 "$SCRIPT1" -sms "$exp"
    echo "[KMER] Finished $exp"
}

# Define function for script 2
run_probound() {
    exp=$1
    echo "[PROBOUND] Starting $exp"
    python3 "$SCRIPT2" -sms "$exp"
    echo "[PROBOUND] Finished $exp"
}

# Export functions and variables
export -f run_kmer
export -f run_probound
export SCRIPT1 SCRIPT2

# Run both in parallel (background & foreground)
echo "$EXPERIMENTS" | parallel -j "$CORES_PER_SCRIPT" run_kmer {} &
echo "$EXPERIMENTS" | parallel -j "$CORES_PER_SCRIPT" run_probound {}

# Wait for all background jobs to finish
wait

echo "All experiments completed for both scripts."