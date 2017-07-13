#!/usr/bin/env sh
set -evx

ROOTDIR=$(dirname `which $0`)
INPREF=$1
OUTDIR=$2
SIM=$3

mkdir -p "${OUTDIR}"

if [ $SIM -gt 62 ] ; then cdhitwordsize=5 ; else cdhitwordsize=3 ; fi
cdhitestwordsize=10

function resetfile() {
    RETVAL=$?
    if [ "${RETVAL}" == 124 ] ; then 
        echo "TIMEOUT-ERROR-FOR-FILE: $1";
        rm "$1" && touch "$1";
        echo true;
    else
        echo false;
    fi
}

function run_mine_with_infastafile_seqid() {
    time -p {
        echo "run_mine_with_infastafile_seqid($1, $2, $3) began at $(date)"
        date; cat "$1.faa"            | "${ROOTDIR}"/bin/fastaseqs-to-distmatrix.out --edsim $3 > "$2-$3.distmatrix"
        date; cat "$2-$3.distmatrix"  | "${ROOTDIR}"/bin/linsetcover.out                        > "$2-$3.ordsetcover"
        date; cat "$2-$3.ordsetcover" | "${ROOTDIR}"/bin/setcover-ords-to-hdrs.out $1.faa       > "$2-$3.hdrsetcover-clu.tsv"
        echo "run_mine_with_infastafile_seqid($1, $2, $3) ended at $(date)"
    }
}

function run_quaclust_with_infastafile_seqid() {
    time -p {
        echo "run_mmseqsclust_with_infastafile_seqid($1, $2, $3) began at $(date)"
        rm -r "$2-$3.quaclust-tmpdir" || true
        mkdir "$2-$3.quaclust-tmpdir"
        rm "$2-$3.quaclust-clu" || true
        rm "$2-$3.quaclust-clu.tsv" || true
        date; timeout $4 "${ROOTDIR}"/benchmark/bin/mmseqs cluster "$1.mmseqs-db" "$2-$3.quaclust-clu" "$2-$3.quaclust-tmpdir" --min-seq-id 0.$3
        if [ $(resetfile "$2-$3.quaclust-clu.tsv") == true ] ; then quaclust_timeout=1; return 0 ; fi
        date; "${ROOTDIR}"/benchmark/bin/mmseqs createtsv             "$1.mmseqs-db" "$1.mmseqs-db"       "$2-$3.quaclust-clu" "$2-$3.quaclust-clu.tsv"
        echo "run_mmseqsclust_with_infastafile_seqid($1, $2, $3) ended at $(date)"
    }
}

function run_linclust_with_infastafile_seqid() {
    time -p {
        echo "run_linclust_with_infastafile_seqid($1, $2, $3) began at $(date)"
        rm -r "$2-$3.linclust-tmpdir" || true
        mkdir "$2-$3.linclust-tmpdir"
        rm "$2-$3.linclust-clu" || true
        rm "$2-$3.linclust-clu.tsv" || true
        date; "${ROOTDIR}"/benchmark/bin/mmseqs linclust         "$1.mmseqs-db" "$2-$3.linclust-clu" "$2-$3.linclust-tmpdir" --min-seq-id 0.$3
        date; "${ROOTDIR}"/benchmark/bin/mmseqs createtsv        "$1.mmseqs-db" "$1.mmseqs-db" "$2-$3.linclust-clu" "$2-$3.linclust-clu.tsv"
        echo "run_linclust_with_infastafile_seqid($1, $2, $3) ended at $(date)"
    }
}

function run_kclust_with_infastafile_seqid() {
    time -p {
        echo "run_kclust_with_infastafile_seqid($1, $2, $3) began at $(date)"
        date; timeout $4 "${ROOTDIR}"/benchmark/bin/kClust -i "$1.faa" -d "$2-$3.kclust" -s $(python -c "print(round(-0.68506329113924050632 + 0.06025316455696202531*$3, 2))")
        if [ $(resetfile "$2-$3.kclust-clu.tsv") == true ] ; then kclust_timeout=1; return 0 ; fi
        date; cat "$2-$3.kclust/clusters.dmp" | tail -n +2 | "${ROOTDIR}"/bin/setcover-ords-to-hdrs.out "$1.faa" > "$2-$3.kclust-clu.tsv"
        echo "run_kclust_with_infastafile_seqid($1, $2, $3) ended at $(date)"
    }
}

function run_cdhit_with_infastafile_seqid() {
    time -p {
        echo "run_cdhit_with_infastafile_seqid($1, $2, $3) began at $(date)"
        date; timeout $4 "${ROOTDIR}"/benchmark/bin/cd-hit -i "$1.faa" -o "$2_cdhit-M0-T0-d0-s80-c$3-n$cdhitwordsize.faa" -M 0 -T 0 -d 0 -s 0.8 -c 0.$3 -n $cdhitwordsize
        if [ $(resetfile "$2-$3.cdhit-clu.tsv") == true ] ; then cdhit_timeout=1; return 0 ; fi
        date; cat "$2_cdhit-M0-T0-d0-s80-c$3-n$cdhitwordsize.faa.clstr" | "${ROOTDIR}"/benchmark/src/clstr-to-clu-tsv.py > "$2-$3.cdhit-clu.tsv"
        echo "run_cdhit_with_infastafile_seqid($1, $2, $3) ended at $(date)"
    }
}

function run_minenuc_with_infastafile_seqid() {
    time -p {
        echo "run_minenuc_with_infastafile_seqid($1, $2, $3) began at $(date)"
        date; cat "$1.fna"            | "${ROOTDIR}"/bin/fastaseqs-to-distmatrix.out --edsim $3 --isnuc 1 > "$2-$3.distmatrix"
        date; cat "$2-$3.distmatrix"  | "${ROOTDIR}"/bin/linsetcover.out                        > "$2-$3.ordsetcover"
        date; cat "$2-$3.ordsetcover" | "${ROOTDIR}"/bin/setcover-ords-to-hdrs.out $1.fna       > "$2-$3.hdrsetcover-clu.tsv"
        echo "run_minenuc_with_infastafile_seqid($1, $2, $3) ended at $(date)"
    }
}

function run_vsearch_with_infastafile_seqid() {
    time -p {
        echo "run_vsearch_with_infastafile_seqid($1, $2, $3) began at $(date)"
        date; "${ROOTDIR}"/benchmark/bin/vsearch --cluster_smallmem $1.fna --id 0.$3 -uc "$2_vsearch-$3.uc"
        date; cat "$2_vsearch-$3.uc" | "${ROOTDIR}"/benchmark/src/uc-to-clu-tsv.py > "$2_vsearch-$3-clu.tsv"
        echo "run_vsearch_with_infastafile_seqid($1, $2, $3) ended at $(date)"
    }
}

function run_cdhitest_with_infastafile_seqid() {
    if [ $3 -lt 80 ]; then
        rm "$2-$3.cdhitest-clu.tsv" && touch "$2-$3.cdhitest-clu.tsv"
    else
    time -p {
        echo "run_cdhitest_with_infastafile_seqid($1, $2, $3) began at $(date)"
        date; "${ROOTDIR}"/benchmark/bin/cd-hit-est -i "$1.fna" -o "$2_cdhitest-M0-T0-d0-s80-c$3-n$cdhitestwordsize.fna" -M 0 -T 0 -d 0 -s 0.8 -c 0.$3 -n $cdhitestwordsize
        date; cat "$2_cdhitest-M0-T0-d0-s80-c$3-n$cdhitestwordsize.fna.clstr" | "${ROOTDIR}"/benchmark/src/clstr-to-clu-tsv.py > "$2-$3.cdhitest-clu.tsv"
        echo "run_cdhitest_with_infastafile_seqid($1, $2, $3) ended at $(date)"
    }
    fi
}

## run Rfam.seed

run_minenuc_with_infastafile_seqid  "${INPREF}/Rfam.seed_len-revname-sort" "${OUTDIR}/Rfam.seed_len-revname-sort" $SIM 3600
run_vsearch_with_infastafile_seqid  "${INPREF}/Rfam.seed_len-revname-sort" "${OUTDIR}/Rfam.seed_len-revname-sort" $SIM 3600
run_cdhitest_with_infastafile_seqid "${INPREF}/Rfam.seed_len-revname-sort" "${OUTDIR}/Rfam.seed_len-revname-sort" $SIM 3600

cat "${OUTDIR}/Rfam.seed_len-revname-sort-${SIM}.hdrsetcover-clu.tsv" | "${ROOTDIR}"/benchmark/src/pfam-clstr-to-metrics.py
cat "${OUTDIR}/Rfam.seed_len-revname-sort_vsearch-${SIM}-clu.tsv"     | "${ROOTDIR}"/benchmark/src/pfam-clstr-to-metrics.py
cat "${OUTDIR}/Rfam.seed_len-revname-sort-${SIM}.cdhitest-clu.tsv"    | "${ROOTDIR}"/benchmark/src/pfam-clstr-to-metrics.py

## run Pfam-A.seed

run_mine_with_infastafile_seqid     "${INPREF}/Pfam-A.seed_len-revname-sort" "${OUTDIR}/Pfam-A.seed_len-revname-sort" $SIM 3600
run_linclust_with_infastafile_seqid "${INPREF}/Pfam-A.seed_len-revname-sort" "${OUTDIR}/Pfam-A.seed_len-revname-sort" $SIM 3600
run_cdhit_with_infastafile_seqid    "${INPREF}/Pfam-A.seed_len-revname-sort" "${OUTDIR}/Pfam-A.seed_len-revname-sort" $SIM 3600
run_quaclust_with_infastafile_seqid "${INPREF}/Pfam-A.seed_len-revname-sort" "${OUTDIR}/Pfam-A.seed_len-revname-sort" $SIM 3600
run_kclust_with_infastafile_seqid   "${INPREF}/Pfam-A.seed_len-revname-sort" "${OUTDIR}/Pfam-A.seed_len-revname-sort" $SIM 3600

cat "${OUTDIR}/Pfam-A.seed_len-revname-sort-${SIM}.hdrsetcover-clu.tsv" | "${ROOTDIR}"/benchmark/src/pfam-clstr-to-metrics.py
cat "${OUTDIR}/Pfam-A.seed_len-revname-sort-${SIM}.linclust-clu.tsv"    | "${ROOTDIR}"/benchmark/src/pfam-clstr-to-metrics.py
cat "${OUTDIR}/Pfam-A.seed_len-revname-sort-${SIM}.cdhit-clu.tsv"       | "${ROOTDIR}"/benchmark/src/pfam-clstr-to-metrics.py
cat "${OUTDIR}/Pfam-A.seed_len-revname-sort-${SIM}.quaclust-clu.tsv"    | "${ROOTDIR}"/benchmark/src/pfam-clstr-to-metrics.py
cat "${OUTDIR}/Pfam-A.seed_len-revname-sort-${SIM}.kclust-clu.tsv"      | "${ROOTDIR}"/benchmark/src/pfam-clstr-to-metrics.py

## run pdb

run_mine_with_infastafile_seqid     "${INPREF}/pdbent-seqres_len-revname-sort" "${OUTDIR}/pdbent-seqres_len-revname-sort" $SIM 3600
run_linclust_with_infastafile_seqid "${INPREF}/pdbent-seqres_len-revname-sort" "${OUTDIR}/pdbent-seqres_len-revname-sort" $SIM 3600
run_cdhit_with_infastafile_seqid    "${INPREF}/pdbent-seqres_len-revname-sort" "${OUTDIR}/pdbent-seqres_len-revname-sort" $SIM 3600
run_quaclust_with_infastafile_seqid "${INPREF}/pdbent-seqres_len-revname-sort" "${OUTDIR}/pdbent-seqres_len-revname-sort" $SIM 3600
run_kclust_with_infastafile_seqid   "${INPREF}/pdbent-seqres_len-revname-sort" "${OUTDIR}/pdbent-seqres_len-revname-sort" $SIM 3600

cat "${OUTDIR}/pdbent-seqres_len-revname-sort-${SIM}.hdrsetcover-clu.tsv" | "${ROOTDIR}"/benchmark/src/setcover-hdrs-to-tms.py > "${OUTDIR}/pdbent-seqres_len-revname-sort-${SIM}.hdrsetcover-clu.tms"
cat "${OUTDIR}/pdbent-seqres_len-revname-sort-${SIM}.linclust-clu.tsv"    | "${ROOTDIR}"/benchmark/src/setcover-hdrs-to-tms.py > "${OUTDIR}/pdbent-seqres_len-revname-sort-${SIM}.linclust-clu.tms"
cat "${OUTDIR}/pdbent-seqres_len-revname-sort-${SIM}.cdhit-clu.tsv"       | "${ROOTDIR}"/benchmark/src/setcover-hdrs-to-tms.py > "${OUTDIR}/pdbent-seqres_len-revname-sort-${SIM}.cdhit-clu.tms"
cat "${OUTDIR}/pdbent-seqres_len-revname-sort-${SIM}.quaclust-clu.tsv"    | "${ROOTDIR}"/benchmark/src/setcover-hdrs-to-tms.py > "${OUTDIR}/pdbent-seqres_len-revname-sort-${SIM}.quaclust-clu.tms"
cat "${OUTDIR}/pdbent-seqres_len-revname-sort-${SIM}.kclust-clu.tsv"      | "${ROOTDIR}"/benchmark/src/setcover-hdrs-to-tms.py > "${OUTDIR}/pdbent-seqres_len-revname-sort-${SIM}.kclust-clu.tms"

benchmark/src/gen-pdbent-cdf.py \
    "${OUTDIR}/pdbent-seqres_len-revname-sort-${SIM}.hdrsetcover-clu.tms" \
    "${OUTDIR}/pdbent-seqres_len-revname-sort-${SIM}.linclust-clu.tms" \
    "${OUTDIR}/pdbent-seqres_len-revname-sort-${SIM}.cdhit-clu.tms" \
    "${OUTDIR}/pdbent-seqres_len-revname-sort-${SIM}.quaclust-clu.tms" \
    "${OUTDIR}/pdbent-seqres_len-revname-sort-${SIM}.kclust-clu.tms" \
    > "${OUTDIR}/pdbent-seqres_len-revname-sort-${SIM}.cdf"

# cat "${OUTDIR}/pdbent-seqres_len-revname-sort-$1.hdrsetcover-clu.tms" | grep "^tmscore" | awk '{print substr($2, 0, 4)}' | sort | uniq -c | awk '{print $2"\t"$1}'

## run uniprot

timelimit=10000
cdhit_timeout=0
quaclust_timeout=0
kclust_timeout=0
for uniref in uniref100-2011-01 uniref100-2014-01 uniref100-2017-01; do
    run_mine_with_infastafile_seqid     "${INPREF}/${uniref}_len-revname-sort" "${OUTDIR}/${uniref}_len-revname-sort" $SIM $timelimit
    run_linclust_with_infastafile_seqid "${INPREF}/${uniref}_len-revname-sort" "${OUTDIR}/${uniref}_len-revname-sort" $SIM $timelimit
    if [ 0 == $cdhit_timeout    ]; then run_cdhit_with_infastafile_seqid    "${INPREF}/${uniref}_len-revname-sort" "${OUTDIR}/${uniref}_len-revname-sort" $SIM $timelimit; fi
    if [ 0 == $quaclust_timeout ]; then run_quaclust_with_infastafile_seqid "${INPREF}/${uniref}_len-revname-sort" "${OUTDIR}/${uniref}_len-revname-sort" $SIM $timelimit; fi
    if [ 0 == $kclust_timeout   ]; then run_kclust_with_infastafile_seqid   "${INPREF}/${uniref}_len-revname-sort" "${OUTDIR}/${uniref}_len-revname-sort" $SIM $timelimit; fi
    timelimit=$(($timelimit*3))
done


