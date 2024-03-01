#!/bin/bash

mkdir -p 14_PHYLO-ALIGN

# Loop through each file in 13_PHYLO-FASTA directory
for file in 13_PHYLO-FASTA/*; do
  # Extract the filename without the path and cut on '_'
  short_name=$(basename "$file" | cut -d "_" -f1,2)
  
  # Align
  clustalo -i "$file" -o "14_PHYLO-ALIGN/${short_name}-align.faa" --auto
done
