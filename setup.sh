#!/usr/bin/env sh

OUT=$1
mkdir -p ${OUT}/
#find $PWD/benchmark/input/pdbents/ -type f | python benchmark/src/pdbents-to-seqres.py > ${OUT}/pdbent-seqres.faa
#cat ${OUT}/pdbent-seqres.faa    | bin/len-revname-sort.out > ${OUT}/pdbent-seqres_len-revname-sort.faa
for uniref in "uniref100-2011-01 uniref100-2014-01 uniref100-2017-01" Pfam-A.seed; do
    cat "benchmark/input/${uniref}.faa" | bin/len-revname-sort.out > "${OUT}/${uniref}_len-revname-sort.faa"
    benchmark/bin/mmseqs createdb "${OUT}/${uniref}_len-revname-sort.faa" "${OUT}/${uniref}_len-revname-sort.mmseqs-db"
done

#benchmark/bin/mmseqs createdb ${OUT}/pdbent-seqres_len-revname-sort.faa ${OUT}/pdbent-seqres_len-revname-sort.mmseqs-db

