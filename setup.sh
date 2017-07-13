#!/usr/bin/env sh

ROOTDIR=$(dirname `which $0`)

IN=$1 ; OUT=$2 ; mkdir -p ${OUT}/

cat "${IN}/Rfam.seed.fna" | "${ROOTDIR}/bin/len-revname-sort.out" > "${OUT}/Rfam.seed_len-revname-sort.fna"

## find "${IN}/pdbents/" -type f | python "${ROOTDIR}/benchmark/src/pdbents-to-seqres.py" > ${IN}/pdbent-seqres.faa

for faafile in pdbent-seqres Pfam-A.seed uniref100-2011-01 uniref100-2014-01 uniref100-2017-01; do
    cat "${IN}/${faafile}.faa" | "${ROOTDIR}/bin/len-revname-sort.out" > "${OUT}/${faafile}_len-revname-sort.faa"
    "${ROOTDIR}/benchmark/bin/mmseqs" createdb "${OUT}/${faafile}_len-revname-sort.faa" "${OUT}/${faafile}_len-revname-sort.mmseqs-db"
done

