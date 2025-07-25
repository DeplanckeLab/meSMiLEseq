#!/bin/bash

# Define your variables
experiment_name="exp10"
transcription_factor="ZNF891_FL"
kmer="6"
binding_size="12"

# Step 1: 
echo "Running 00_read_in_data.py..."
python3 00_read_in_data.py -sms "$experiment_name" -tf "$transcription_factor" 

# Check if the script completed successfully
if [ $? -eq 0 ]; then
    echo "00_read_in_data.py completed successfully."
else
    echo "00_read_in_data.py failed."
    exit 1
fi

# Step 2: 
echo "Running 01_kmer_analysis.py..."
python3 01_kmer_analysis.py -sms "$experiment_name" -tf "$transcription_factor" -k "$kmer"

# Check if the script completed successfully
if [ $? -eq 0 ]; then
    echo "01_kmer_analysis.py completed successfully."
else
    echo "01_kmer_analysis.py failed."
    exit 1
fi

# Step 3: 
echo "Running 02_fishers_exact_test.py..."
python3 02_fishers_exact_test.py -sms "$experiment_name" -tf "$transcription_factor" -k "$kmer"

# Check if the script completed successfully
if [ $? -eq 0 ]; then
    echo "02_fishers_exact_test.py completed successfully."
else
    echo "02_fishers_exact_test.py failed."
    exit 1
fi

# Step 4: 
echo "Running 03_calculate_ratios.py..."
python3 03_calculate_ratios.py -sms "$experiment_name" -tf "$transcription_factor" -k "$kmer"

# Check if the script completed successfully
if [ $? -eq 0 ]; then
    echo "03_calculate_ratios.py completed successfully."
else
    echo "03_calculate_ratios.py failed."
    exit 1
fi

# Step 5:
echo "Running 04_pyProBound_analysis.py..."
python3 04_pyProBound_analysis.py -sms "$experiment_name" -tf "$transcription_factor" -bs "$binding_size"

# Check if the script completed successfully
if [ $? -eq 0 ]; then
    echo "04_pyProBound_analysis.py completed successfully."
else
    echo "04_pyProBound_analysis.py failed."
    exit 1
fi

echo "Pipeline completed successfully!"
