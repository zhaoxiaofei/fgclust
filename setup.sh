#!/usr/bin/env sh

OUT=$1
mkdir -p ${OUT}/
find $PWD/benchmark/input/pdbents/ -type f | python benchmark/src/pdbents-to-seqres.py > ${OUT}/pdbent-seqres.faa
cat ${OUT}/pdbent-seqres.faa    | bin/len-revname-sort.out > ${OUT}/pdbent-seqres_len-revname-sort.faa
cat benchmark/input/uniref100-2017-03.faa | bin/len-revname-sort.out > ${OUT}/uniref100-2017-03_len-revname-sort.faa
benchmark/bin/mmseqs createdb ${OUT}/pdbent-seqres_len-revname-sort.faa ${OUT}/pdbent-seqres_len-revname-sort.mmseqs-db
benchmark/bin/mmseqs createdb ${OUT}/uniref100-2017-03_len-revname-sort.faa ${OUT}/uniref100-2017-03_len-revname-sort.mmseqs-db

