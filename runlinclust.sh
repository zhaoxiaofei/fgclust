#!/usr/bin/env sh
set -evx

function run_linclust_with_infastafile_seqid() {
    rm -r "$1-$2.mmseq-tmpdir" || true
    mkdir "$1-$2.mmseq-tmpdir"
    rm "$1-$2.mmseq-clu" || true
    rm "$1-$2.mmseq-clu.tsv" || true
    benchmark/bin/mmseqs linclust         "$1.mmseq-db" "$1-$2.mmseq-clu" "$1-$2.mmseq-tmpdir" --min-seq-id 0.$2
    #benchmark/bin/mmseqs createseqfiledb "$1.mmseq-db" "$1.mmseq-clu" "$1.mmseq-clu_seq"
    #benchmark/bin/mmseqs result2flat     "$1.mmseq-db" "$1.mmseq-db" "$1.mmseq-clu_seq" "$1.mmseq-clu_seq.fasta"
    benchmark/bin/mmseqs createtsv        "$1.mmseq-db" "$1.mmseq-db" "$1-$2.mmseq-clu" "$1-$2.mmseq-clu.tsv"
}

if [ $1 -gt 62 ] ; then cdhitwordsize=5 ; else cdhitwordsize=3 ; fi

#benchmark/bin/mmseqs createdb benchmark/output/pdbent-seqres_len-revname-sort.faa benchmark/output/pdbent-seqres_len-revname-sort.mmseq-db

# run uniprot

run_linclust_with_infastafile_seqid benchmark/output/uniref100-2017-03_len-revname-sort $1


