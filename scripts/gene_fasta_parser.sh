#!/bin/bash

# Make directory for FASTA-GENES, if it doesn't already exist
mkdir -p 10_FASTA-GENES

# Change to the 10_FASTA-GENES directory
cd 10_FASTA-GENES

# Loop through samples.csv and run gffParse.pl for each line
while IFS=, read -r taxa genome gff; do
  # Run parser for each genome
  gffParse.pl -c -i ../"$genome" -g ../"$gff" -b ../"$taxa" -d "$taxa" -p -f CDS
done < ../samples.csv

# Change back to the original directory
cd ..
