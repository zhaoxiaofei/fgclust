#!/usr/bin/env sh
set -evx

INPREF=$1
OUTDIR=$2
SIM=$3

mkdir -p "${OUTDIR}"

if [ $SIM -gt 62 ] ; then cdhitwordsize=5 ; else cdhitwordsize=3 ; fi

function run_mine_with_infastafile_seqid() {
    echo "run_mine_with_infastafile_seqid($1, $2, $3) began at $(date)"
    date; cat "$1.faa"            | bin/fastaseqs-to-distmatrix.out --edsim $3 > "$2-$3.distmatrix"
    date; cat "$2-$3.distmatrix"  | bin/linsetcover.out                        > "$2-$3.ordsetcover"
    date; cat "$2-$3.ordsetcover" | bin/setcover-ords-to-hdrs.out $1.faa       > "$2-$3.hdrsetcover-clu.tsv"
    echo "run_mine_with_infastafile_seqid($1, $2, $3) ended at $(date)"
}

function run_linclust_with_infastafile_seqid() {
    echo "run_linclust_with_infastafile_seqid($1, $2, $3) began at $(date)"
    rm -r "$2-$3.mmseqs-tmpdir" || true
    mkdir "$2-$3.mmseqs-tmpdir"
    rm "$2-$3.mmseqs-clu" || true
    rm "$2-$3.mmseqs-clu.tsv" || true
    date; benchmark/bin/mmseqs linclust         "$1.mmseqs-db" "$2-$3.mmseqs-clu" "$2-$3.mmseqs-tmpdir" --min-seq-id 0.$3
    date; benchmark/bin/mmseqs createtsv        "$1.mmseqs-db" "$1.mmseqs-db" "$2-$3.mmseqs-clu" "$2-$3.mmseqs-clu.tsv"
    echo "run_linclust_with_infastafile_seqid($1, $2, $3) ended at $(date)"
}

function run_cdhit_with_infastafile_seqid() {
    echo "run_cdhit_with_infastafile_seqid($1, $2, $3) began at $(date)"
    date; benchmark/bin/cd-hit -i "$1.faa" -o "$2_cdhit-M0-T0-d0-s80-c$3-n$cdhitwordsize.faa" -M 0 -T 0 -d 0 -s 0.8 -c 0.$3 -n $cdhitwordsize
    date; cat "$2_cdhit-M0-T0-d0-s80-c$3-n$cdhitwordsize.faa.clstr" | benchmark/src/clstr-to-clu-tsv.py > "$2-$3.cdhit-clu.tsv"
    echo "run_cdhit_with_infastafile_seqid($1, $2, $3) ended at $(date)"
}

## run pdb

run_mine_with_infastafile_seqid     "${INPREF}/pdbent-seqres_len-revname-sort" "${OUTDIR}/pdbent-seqres_len-revname-sort" $SIM
run_linclust_with_infastafile_seqid "${INPREF}/pdbent-seqres_len-revname-sort" "${OUTDIR}/pdbent-seqres_len-revname-sort" $SIM
run_cdhit_with_infastafile_seqid    "${INPREF}/pdbent-seqres_len-revname-sort" "${OUTDIR}/pdbent-seqres_len-revname-sort" $SIM

cat "${OUTDIR}/pdbent-seqres_len-revname-sort-${SIM}.hdrsetcover-clu.tsv" | benchmark/src/setcover-hdrs-to-tms.py > "${OUTDIR}/pdbent-seqres_len-revname-sort-${SIM}.hdrsetcover-clu.tms"
cat "${OUTDIR}/pdbent-seqres_len-revname-sort-${SIM}.mmseqs-clu.tsv"      | benchmark/src/setcover-hdrs-to-tms.py > "${OUTDIR}/pdbent-seqres_len-revname-sort-${SIM}.mmseqs-clu.tms"
cat "${OUTDIR}/pdbent-seqres_len-revname-sort-${SIM}.cdhit-clu.tsv"       | benchmark/src/setcover-hdrs-to-tms.py > "${OUTDIR}/pdbent-seqres_len-revname-sort-${SIM}.cdhit-clu.tms"

benchmark/src/gen-pdbent-table.py \
    "${OUTDIR}/pdbent-seqres_len-revname-sort-${SIM}.hdrsetcover-clu.tms" \
    "${OUTDIR}/pdbent-seqres_len-revname-sort-${SIM}.mmseqs-clu.tms" \
    "${OUTDIR}/pdbent-seqres_len-revname-sort-${SIM}.cdhit-clu.tms" > "${OUTDIR}/pdbent-seqres_len-revname-sort-${SIM}.distribution"
benchmark/src/gen-pdbent-cdf.py \
    "${OUTDIR}/pdbent-seqres_len-revname-sort-${SIM}.hdrsetcover-clu.tms" \
    "${OUTDIR}/pdbent-seqres_len-revname-sort-${SIM}.mmseqs-clu.tms" \
    "${OUTDIR}/pdbent-seqres_len-revname-sort-${SIM}.cdhit-clu.tms" > "${OUTDIR}/pdbent-seqres_len-revname-sort-${SIM}.cdf"

# cat "${OUTDIR}/pdbent-seqres_len-revname-sort-$1.hdrsetcover-clu.tms" | grep "^tmscore" | awk '{print substr($2, 0, 4)}' | sort | uniq -c | awk '{print $2"\t"$1}'

## run uniprot

run_mine_with_infastafile_seqid     "${INPREF}/uniref100-2017-03_len-revname-sort" "${OUTDIR}/uniref100-2017-03_len-revname-sort" $SIM
run_linclust_with_infastafile_seqid "${INPREF}/uniref100-2017-03_len-revname-sort" "${OUTDIR}/uniref100-2017-03_len-revname-sort" $SIM
run_cdhit_with_infastafile_seqid    "${INPREF}/uniref100-2017-03_len-revname-sort" "${OUTDIR}/uniref100-2017-03_len-revname-sort" $SIM

