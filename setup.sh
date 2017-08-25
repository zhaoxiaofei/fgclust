#!/usr/bin/env sh

ROOTDIR=$(dirname `which $0`)

IN=$1 ; OUT=$2 ; mkdir -p ${OUT}/

cat "${IN}/Rfam.seed.fna" | "${ROOTDIR}/bin/len-revname-sort.out" > "${OUT}/Rfam.seed_shuf.fna"

## find "${IN}/pdbents/" -type f | python "${ROOTDIR}/benchmark/src/pdbents-to-seqres.py" > ${IN}/pdbent-seqres.faa

for faafile in pdbent-seqres Pfam-A.seed uniref100-02 uniref100-12 uniref100-2011-01 uniref100-2014-01 uniref100-2017-01; do
#for faafile in uniref100-2; do
    date; cat "${IN}/${faafile}.faa" | "${ROOTDIR}/bin/len-revname-sort.out" --israndom 1 > "${OUT}/${faafile}_shuf.faa"
    date; "${ROOTDIR}/benchmark/bin/mmseqs" createdb "${OUT}/${faafile}_shuf.faa" "${OUT}/${faafile}_shuf.mmseqs-db" ;
    date
done

if false; then
"${ROOTDIR}/benchmark/src/subsample-half-fasta.py" < "${OUT}/uniref100-2017-01_shuf.faa"        > "${OUT}/down2x-uniref100-2017-01_shuf.faa" 
"${ROOTDIR}/benchmark/src/subsample-half-fasta.py" < "${OUT}/down2x-uniref100-2017-01_shuf.faa" > "${OUT}/down4x-uniref100-2017-01_shuf.faa" 
"${ROOTDIR}/benchmark/src/subsample-half-fasta.py" < "${OUT}/down4x-uniref100-2017-01_shuf.faa" > "${OUT}/down8x-uniref100-2017-01_shuf.faa" 
"${ROOTDIR}/benchmark/src/subsample-half-fasta.py" < "${OUT}/down8x-uniref100-2017-01_shuf.faa" > "${OUT}/down16x-uniref100-2017-01_shuf.faa" 
fi

