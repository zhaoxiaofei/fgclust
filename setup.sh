#!/usr/bin/env sh

mkdir -p benchmark/setup-out/
find $PWD/benchmark/input/pdbents/ -type f | python benchmark/src/pdbents-to-seqres.py > benchmark/setup-out/pdbent-seqres.faa
cat benchmark/setup-out/pdbent-seqres.faa    | bin/len-revname-sort.out > benchmark/setup-out/pdbent-seqres_len-revname-sort.faa
cat benchmark/input/uniref100-2017-03.faa | bin/len-revname-sort.out > benchmark/setup-out/uniref100-2017-03_len-revname-sort.faa
benchmark/bin/mmseqs createdb benchmark/setup-out/pdbent-seqres_len-revname-sort.faa benchmark/setup-out/pdbent-seqres_len-revname-sort.mmseqs-db
benchmark/bin/mmseqs createdb benchmark/setup-out/uniref100-2017-03_len-revname-sort.faa benchmark/setup-out/uniref100-2017-03_len-revname-sort.mmseqs-db

