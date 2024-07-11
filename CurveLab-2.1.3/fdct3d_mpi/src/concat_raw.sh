#!/bin/bash

if [ "$#" -lt 1 ]; then
        echo "Need an input (number of raw files)"
    exit 1
fi

# Define the output file
output_file="outdata_3d_mpi.raw"
prefix="outdata_3d_mpi_"
suffix=".raw"
arg1=$1
last_ind=$((arg1 - 1))

# Create or empty the output file
> "$output_file"

# Loop through the list of input files and concatenate them
for num in $(seq 0 $last_ind); do
    file="${prefix}${num}${suffix}"
    if [ -f "$file" ]; then
        cat "$file" >> "$output_file"
    else
        echo "Warning: $file does not exist."
    fi
done

# Print a success message
echo "Files concatenated into $output_file"