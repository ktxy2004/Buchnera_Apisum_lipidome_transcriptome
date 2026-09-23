#!/bin/bash

# Destination directory to store renamed quant.sf files
destination_dir=~/Desktop/rnaseq/rna_analysis_aphid/quant_sf

# Create the destination directory if it doesn't exist
mkdir -p "$destination_dir"

# Loop through each directory containing quant.sf files
for folder in ~/Desktop/rnaseq/rna_analysis_aphid/*/; do
    # Extract the folder name
    folder_name=$(basename "$folder")
    
    # Copy and rename each quant.sf file
    cp "$folder"/*.sf "$destination_dir"/"${folder_name%_quant}_quant.sf"
done
