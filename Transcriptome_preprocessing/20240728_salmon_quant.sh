#!/bin/bash

#generated a decoy.text file extracted names of the genome targets

grep "^>" <(gunzip -c GCF_005508785.2_pea_aphid_22Mar2018_4r6ur_v2_genomic.fna.gz) | cut -d " " -f 1 > decoys.txt        
sed -i.bak -e 's/>//g' decoys.txt 


#####concatenate transcriptome and genome reference file (from NCBI) genome and rna.fna files

cat rna.fna.gz GCF_005508785.2_pea_aphid_22Mar2018_4r6ur_v2_genomic.fna.gz > ApisumGentrome.fa.gz


###salmon indexing 
salmon index -t ApisumGentrome.fa.gz -d decoys.txt -p 12 -i salmon_index --gencode


# Path to Salmon index
index="/Users/kathrine/Desktop/20240725_Apisum_reads/20240728_Apisum_nonrrna/salmon_index"

# Path to gene-to-transcript mapping file
gene_map_file="/Users/kathrine/Desktop/20240725_Apisum_reads/20240728_Apisum_nonrrna/genomic.gtf"

# Output directory
output_dir="20240728_Apisum_nonrrna_salmonOut"

# Loop over each pair of forward and reverse read files
for forward_file in *_nonrrna_1.fq; do
    # Extract the base filename without extension
    base_name="${forward_file%_nonrrna_1.fq}"

    # Define the corresponding reverse read file
    reverse_file="${base_name}_nonrrna_2.fq"

    # Run Salmon quantification for the current pair of files
    salmon quant -i "$index" -l A --geneMap "$gene_map_file" \
        -1 "$forward_file" \
        -2  "$reverse_file" \
        --validateMappings -o "${output_dir}/${base_name}_quant"
done
