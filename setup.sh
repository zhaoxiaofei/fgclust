#!/usr/bin/env sh

#find $PWD/benchmark/input/pdbents/ -type f | python benchmark/src/pdbents-to-seqres.py > benchmark/output/pdbent-seqres.faa
cat benchmark/output/pdbent-seqres.faa    | bin/len-revname-sort.out > benchmark/output/pdbent-seqres_len-revname-sort.faa
cat benchmark/input/uniref100-2017-03.faa | bin/len-revname-sort.out > benchmark/output/uniref100-2017-03_len-revname-sort.faa
benchmark/bin/mmseqs createdb benchmark/output/pdbent-seqres_len-revname-sort.faa benchmark/output/pdbent-seqres_len-revname-sort.mmseqs-db
benchmark/bin/mmseqs createdb benchmark/output/uniref100-2017-03_len-revname-sort.faa benchmark/output/uniref100-2017-03_len-revname-sort.mmseqs-db

