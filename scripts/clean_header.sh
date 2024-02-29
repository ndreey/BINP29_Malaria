#!/bin/bash

# Set directory
dir="10_FASTA-GENES"

for taxa in Ht Pb Pc Pf Pk Pv Py Tg; do
  fix="_$taxa"

  # Rename to raw.
  mv ${dir}/${taxa}.faa ${dir}/raw-${taxa}.faa

  # Fix header to keep gene number and taxa, and remove asterisks.
  awk -v fix="$fix" '/^>/ {print $1 fix; next} {print}' ${dir}/raw-${taxa}.faa \
   | sed 's/\*//g' > ${dir}/${taxa}-clean.faa

done
