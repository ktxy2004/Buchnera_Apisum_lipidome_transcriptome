#!/bin/bash
#export PATH="$HOME/Documents/bbmap:$PATH"  
# Path to reference genome
refApisum="/Users/kathrine/Desktop/rnaseq/rna_analysis_aphid/aphid_rna/20240408_trimmomatic_paired/GCF_005508785.2_pea_aphid_22Mar2018_4r6ur_v2_genomic.fna.gz"
refBuch="/Users/kathrine/Desktop/rnaseq/rna_analysis_aphid/aphid_rna/20240408_trimmomatic_paired/GCA_000009605.1_ASM960v1_genomic.fna.gz"

# Output directory
output_dir="out_bbsplit2"

# Memory allocation
memory="700g"

# Number of threads (change 2 to the desired value)
threads="6"

# Loop over each pair of forward and reverse read files
for forward_file in *_R1_paired.fastq.gz; do
    # Extract the base filename without extension
    base_name="${forward_file%_R1_paired.fastq.gz}"

    # Define the corresponding reverse read file
    reverse_file="${base_name}_R2_paired.fastq.gz"

    # Run bbsplit for the current pair of files
    bbsplit.sh -Xmx"$memory" threads="$threads" \
        in1="$forward_file" \
        in2="$reverse_file" \
        ref_x="$refApisum" \
        ref_y="$refBuch" \
        basename="${output_dir}/${base_name}_%_#.fq" \
        outu="${output_dir}/${base_name}_cleaned_#.fq"
done
