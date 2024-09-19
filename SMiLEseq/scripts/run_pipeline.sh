#!/bin/bash

# Step 1: Run the first script
echo "Running 00_readin_data.py..."
python3 00_readin_data.py

# Check if the script completed successfully
if [ $? -eq 0 ]; then
    echo "00_readin_data.py completed successfully."
else
    echo "00_readin_data.py failed."
    exit 1
fi

# Step 2: Run the second script
echo "Running 01_fishers_exact_test.py..."
python3 01_fishers_exact_test.py

# Check if the script completed successfully
if [ $? -eq 0 ]; then
    echo "01_fishers_exact_test.py completed successfully."
else
    echo "01_fishers_exact_test.py failed."
    exit 1
fi

# Step 3: Run the third script
echo "Running 02_extract_significant_fastqs.py..."
python3 02_extract_significant_fastqs.py

# Check if the script completed successfully
if [ $? -eq 0 ]; then
    echo "02_extract_significant_fastqs.py completed successfully."
else
    echo "02_extract_significant_fastqs.py failed."
    exit 1
fi

# Step 4: Run the fourth script
echo "Running 03_ProBound_significant_TFs.py..."
python3 03_ProBound_significant_TFs.py

# Check if the script completed successfully
if [ $? -eq 0 ]; then
    echo "03_ProBound_significant_TFs.py completed successfully."
else
    echo "03_ProBound_significant_TFs.py failed."
    exit 1
fi

echo "Pipeline completed successfully!"
