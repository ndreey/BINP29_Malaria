#!/bin/bash

tsv="run_apicomplexa_odb10/full_table.tsv"
for taxa in Ht Pb Pc Pf Pk Pv Py; do

cat 12_BUSCO-ORTHO/BUSCO_$taxa-clean.faa/$tsv | grep -w "Complete" | \
 cut -f 1 >> 12_BUSCO-ORTHO/no-Tg_BUSCO-ID.txt

done

numb=$(cat 12_BUSCO-ORTHO/all_BUSCO-ID.txt | sort | uniq -c | grep -c -w "7")

echo "  Number of orthologous genes shared: $numb"
