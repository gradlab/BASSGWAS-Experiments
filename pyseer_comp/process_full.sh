#!/bin/bash
mkdir $1'_full_pyseer'
Rscript comp_pyseer2.R ../data/$2/data_bin.csv ../data/$2 ./$1'_full_pyseer' ../annotate/references.txt
