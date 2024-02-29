#!/bin/bash

# Set directory
dir="10_FASTA-GENES"

# Create directory if it doesnt exist
mkdir -p 12_BUSCO-ORTHO

for taxa in Ht Pb Pc Pf Pk Pv Py Tg; do

  busco -i $dir/$taxa-clean.faa -m prot -c 64 \
	  --out_path 12_BUSCO-ORTHO -l apicomplexa_odb10

done
