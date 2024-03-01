#!/bin/bash

mkdir -p 15_IQTREE/
cd 15_IQTREE/

# Loop through each file
for file in ../14_PHYLO-ALIGN/*; do
  # Extract the filename without the path and cut on '_'
  short_name=$(basename "$file" | cut -d "-" -f1)
  
  mkdir -p $short_name
  
  # Model
  iqtree2 -nt 32 -s $file -m MFP -B 1000 --prefix ${short_name}/${short_name}
done