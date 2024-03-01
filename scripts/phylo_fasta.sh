#!/bin/bash

phylo="13_PHYLO-FASTA"
dir="12_BUSCO-ORTHO"
tsv="run_apicomplexa_odb10/full_table.tsv"

mkdir -p $phylo

# Read through the BUSCO IDs
while read -r id; do

  # Append each sequence record for each taxa to BUSCO ID fasta file.
  for taxa in Ht Pb Pc Pf Pk Pv Py Tg; do

    # Get the header
    header=$(cat $dir/BUSCO_$taxa-clean.faa/$tsv | grep $id | cut -f3)
    
    # Get name of gene without whitespace or other special characters
    name=$(cat $dir/BUSCO_$taxa-clean.faa/$tsv | grep $id | cut -f7 | \
      tr " " "_" | sed "s/[,()]//g")

    fasta="10_FASTA-GENES/$taxa-clean.faa"

    # Append the fasta record (header + sequence)
    cat $fasta | grep -A 1 -w "^>$header" >> "${phylo}/${id}_${name}.faa"
  
  done
done < $dir/ALL-SHARE-ORTHO.txt