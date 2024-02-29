#!/bin/bash

tsv="run_apicomplexa_odb10/full_table.tsv"
for taxa in Ht Pb Pc Pf Pk Pv Py Tg; do

cat 12_BUSCO-ORTHO/BUSCO_$taxa-clean.faa/$tsv | grep -w "Complete" | \
 cut -f 1 >> 12_BUSCO-ORTHO/all_BUSCO-ID.txt

done