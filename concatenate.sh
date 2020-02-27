#!/bin/bash

cd results/execution

#head -1 results_level1.csv > results_level$1.csv
#for f in results_level$1.csv; do (cat "${f}"; echo) >> results_all.csv; done

cat results_level*.csv >> results_all.csv

#sort -u results_all.csv >> results_all_sorted.csv
#rm -r *_$1.csv