#!/bin/bash

# Get the target folder from the first argument
target_dir="$1"

# Check the folder exists
if [ ! -d "$target_dir" ]; then
    echo "Error: Directory '$target_dir' does not exist."
    exit 1
fi

cd "$target_dir" || exit

# Rename all files matching pattern: *_6kmer_*.csv, *_7kmer_*.csv, etc.
for file in *[0-9]kmer*.csv; do
    # Use sed to remove 'k' in '6kmer' → '6mer'
    newname=$(echo "$file" | sed -E 's/([0-9])kmer/\1mer/')
    if [[ "$file" != "$newname" ]]; then
        mv "$file" "$newname"
        echo "Renamed: $file → $newname"
    fi
done