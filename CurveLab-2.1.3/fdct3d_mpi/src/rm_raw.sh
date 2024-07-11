#!/bin/bash

if [ "$#" -lt 1 ]; then
        echo "Need an input (number of files to remove)"
    exit 1
fi

# Define the output file
prefix="outdata_3d_mpi_"
suffix=".raw"
arg1=$1
last_ind=$((arg1 - 1))


# Loop through the list of input files and concatenate them
for num in $(seq 0 $last_ind); do
    file="${prefix}${num}${suffix}"
    if [ -f "$file" ]; then
        rm "$file"
        echo "Removed $file"
    else
        echo "Warning: $file does not exist."
    fi
done

# Print a success message
echo "Output files 0 to $last_ind removed"
