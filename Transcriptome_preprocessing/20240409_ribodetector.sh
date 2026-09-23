#!/bin/bash
# Define the desired chunk size 
chunk_size=256

# Loop over each pair of forward and reverse read files
for forward_file in *_x_1.fq; do
    # Extract the base filename without extension
    base_name="${forward_file%_x_1.fq}"

    # Define the corresponding reverse read file
    reverse_file="${base_name}_x_2.fq"

    # Define the output filename
    output_file_1="${base_name}_nonrrna_1.fq"
    output_file_2="${base_name}_nonrrna_2.fq"

    # Run ribodetector_cpu for the current pair of files
    ribodetector_cpu -t 20 -l 70 -i "$forward_file" "$reverse_file" --e none -o "$output_file_1" "$output_file_2"  --chunk_size "$chunk_size"
done
