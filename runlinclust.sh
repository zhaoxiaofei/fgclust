#!/usr/bin/env sh
set -evx

if [ $1 -gt 62 ] ; then cdhitwordsize=5 ; else cdhitwordsize=3 ; fi

function run_mine_with_infastafile_seqid() {
    echo "run_mine_with_infastafile_seqid($1, $2) began at $(date)"
    date; cat $1.faa            | bin/fastaseqs-to-distmatrix.out --edsim $2 > $1-$2.distmatrix
    date; cat $1-$2.distmatrix  | bin/linsetcover.out                        > $1-$2.ordsetcover
    date; cat $1-$2.ordsetcover | bin/setcover-ords-to-hdrs.out $1.faa       > $1-$2.hdrsetcover-clu.tsv
    echo "run_mine_with_infastafile_seqid($1, $2) ended at $(date)"
}

function run_linclust_with_infastafile_seqid() {
    echo "run_linclust_with_infastafile_seqid($1, $2) began at $(date)"
    rm -r "$1-$2.mmseqs-tmpdir" || true
    mkdir "$1-$2.mmseqs-tmpdir"
    rm "$1-$2.mmseqs-clu" || true
    rm "$1-$2.mmseqs-clu.tsv" || true
    date; benchmark/bin/mmseqs linclust         "$1.mmseqs-db" "$1-$2.mmseqs-clu" "$1-$2.mmseqs-tmpdir" --min-seq-id 0.$2
    date; benchmark/bin/mmseqs createtsv        "$1.mmseqs-db" "$1.mmseqs-db" "$1-$2.mmseqs-clu" "$1-$2.mmseqs-clu.tsv"
    echo "run_linclust_with_infastafile_seqid($1, $2) ends at $(date)"
}

function run_cdhit_with_infastafile_seqid() {
    echo "run_cdhit_with_infastafile_seqid($1, $2) began at $(date)"
    date; benchmark/bin/cd-hit -i $1.faa -o $1_cdhit-M0-T0-d0-s80-c$2-n$cdhitwordsize.faa -M 0 -T 0 -d 0 -s 0.8 -c 0.$2 -n $cdhitwordsize
    date; cat $1_cdhit-M0-T0-d0-s80-c$2-n$cdhitwordsize.faa.clstr | benchmark/src/clstr-to-clu-tsv.py > $1_cdhit-M0-T0-d0-s80-c$2-n$cdhitwordsize.cdhit-clu.tsv 
    echo "run_cdhit_with_infastafile_seqid($1, $2) ended at $(date)"
}


# run pdb

run_mine_with_infastafile_seqid     benchmark/output/pdbent-seqres_len-revname-sort $1
run_linclust_with_infastafile_seqid benchmark/output/pdbent-seqres_len-revname-sort $1
run_cdhit_with_infastafile_seqid    benchmark/output/pdbent-seqres_len-revname-sort $1

cat benchmark/output/pdbent-seqres_len-revname-sort-$1.hdrsetcover-clu.tsv | benchmark/src/setcover-hdrs-to-tms.py > pdbent-seqres_len-revname-sort-$1.hdrsetcover-clu.tms
cat benchmark/output/pdbent-seqres_len-revname-sort-$1.mmseqs-clu.tsv      | benchmark/src/setcover-hdrs-to-tms.py > pdbent-seqres_len-revname-sort-$1.mmseqs-clu.tms
cat benchmark/output/pdbent-seqres_len-revname-sort-$1.cdhit-clu.tsv       | benchmark/src/setcover-hdrs-to-tms.py > pdbent-seqres_len-revname-sort-$1.cdhit-clu.tms

# run uniprot

run_mine_with_infastafile_seqid     benchmark/output/uniref100-2017-03_len-revname-sort $1
run_linclust_with_infastafile_seqid benchmark/output/uniref100-2017-03_len-revname-sort $1
run_cdhit_with_infastafile_seqid    benchmark/output/uniref100-2017-03_len-revname-sort $1


