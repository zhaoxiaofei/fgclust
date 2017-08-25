#!/usr/bin/env sh
set -evx

ROOTDIR=$(dirname `which $0`)
INPREF=$1
OUTDIR=$2
CSVSIM=$3
MINEONLY=$4

mkdir -p "${OUTDIR}"

function resetfile() {
    timeoutret=$?
    set -e
    if [ "${timeoutret}" == 124 ] ; then 
        ISTIMEOUT=true
        echo "TIMEOUT-ERROR-FOR-FILE: $1";
        rm "$1" || true ; touch "$1";
        echo true;
    else
        ISTIMEOUT=false
        echo false;
    fi
}

function run_mine_with_infastafile_seqid() {
    time -p {
        echo "run_mine_with_infastafile_seqid($1, $2, $3, $4) began at $(date)"
        if [ $4 -lt 20 ]; then flag="--alphasize_sim $4"; else flag=""; fi
        date; cat "$1.faa"            | "${ROOTDIR}"/bin/fastaseqs-to-distmatrix.out --edsim $3 $flag > "$2-$3.distmatrix"
        date; cat "$2-$3.distmatrix"  | "${ROOTDIR}"/bin/linsetcover.out                              > "$2-$3.ordsetcover"
        date; cat "$2-$3.ordsetcover" | "${ROOTDIR}"/bin/setcover-ords-to-hdrs.out $1.faa             > "$2-$3.hdrsetcover-clu.tsv"
        echo "run_mine_with_infastafile_seqid($1, $2, $3) ended at $(date)"
    } | tee "$2-$3.mine.time"
}

function run_quaclust_with_infastafile_seqid() {
    time -p {
        echo "run_mmseqsclust_with_infastafile_seqid($1, $2, $3) began at $(date)"
        rm -r "$2-$3.quaclust-tmpdir" || true
        mkdir "$2-$3.quaclust-tmpdir"
        rm "$2-$3.quaclust-clu" || true
        rm "$2-$3.quaclust-clu.tsv" || true
        set +e
        date; timeout $4 "${ROOTDIR}"/benchmark/bin/mmseqs cluster "$1.mmseqs-db" "$2-$3.quaclust-clu" "$2-$3.quaclust-tmpdir" --min-seq-id 0.$3
        resetfile "$2-$3.quaclust-clu.tsv"
        if [ "${ISTIMEOUT}" == true ] ; then quaclust_timeout=true; return 0 ; fi
        date; "${ROOTDIR}"/benchmark/bin/mmseqs createtsv             "$1.mmseqs-db" "$1.mmseqs-db"       "$2-$3.quaclust-clu" "$2-$3.quaclust-clu.tsv"
        echo "run_mmseqsclust_with_infastafile_seqid($1, $2, $3) ended at $(date)"
    } | tee "$2-$3.quaclust.time"
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
    } | tee "$2-$3.linclust.time"
}

function run_kclust_with_infastafile_seqid() {
    time -p {
        echo "run_kclust_with_infastafile_seqid($1, $2, $3) began at $(date)"
        set +e
        date; timeout $4 "${ROOTDIR}"/benchmark/bin/kClust -i "$1.faa" -d "$2-$3.kclust" -s $(python -c "print(round(-0.68506329113924050632 + 0.06025316455696202531*$3, 2))")
        resetfile "$2-$3.kclust-clu.tsv"
        if [ "${ISTIMEOUT}" == true ] ; then kclust_timeout=true; return 0 ; fi
        date; cat "$2-$3.kclust/clusters.dmp" | tail -n +2 | awk '{print $2"\t"$1"\t"1}' | "${ROOTDIR}"/bin/setcover-ords-to-hdrs.out "$1.faa" > "$2-$3.kclust-clu.tsv"
        echo "run_kclust_with_infastafile_seqid($1, $2, $3) ended at $(date)"
    } | tee "$2-$3.kclust.time"
}

function run_cdhit_with_infastafile_seqid() {
    if [ "${cdhit_timeout}" == true ]; then
        echo "run_cdhit_with_infastafile_seqid($1, $2, $3) already timed out in previous iteration."; 
        return 0;
    fi
    cdhitwordsize=5
    if [ $SIM -lt 70 ] ; then cdhitwordsize=4; fi
    if [ $SIM -lt 60 ] ; then cdhitwordsize=3; fi
    time -p {
        echo "run_cdhit_with_infastafile_seqid($1, $2, $3) began at $(date)"
        set +e
        date; timeout $4 "${ROOTDIR}"/benchmark/bin/cd-hit -i "$1.faa" -o "$2_cdhit-M0-T0-d0-s80-c$3-n$cdhitwordsize.faa" -M 0 -T 0 -d 0 -s 0.8 -c 0.$3 -n $cdhitwordsize
        resetfile "$2-$3.cdhit-clu.tsv"
        if [ "${ISTIMEOUT}" == true ] ; then cdhit_timeout=true; return 0 ; fi
        date; cat "$2_cdhit-M0-T0-d0-s80-c$3-n$cdhitwordsize.faa.clstr" | "${ROOTDIR}"/benchmark/src/clstr-to-clu-tsv.py > "$2-$3.cdhit-clu.tsv"
        echo "run_cdhit_with_infastafile_seqid($1, $2, $3) ended at $(date)"
    } | tee "$2-$3.cdhit.time" 
}

function run_minenuc_with_infastafile_seqid() {
    time -p {
        echo "run_minenuc_with_infastafile_seqid($1, $2, $3) began at $(date)"
        date; cat "$1.fna"            | "${ROOTDIR}"/bin/fastaseqs-to-distmatrix.out --edsim $3 --isnuc 1 > "$2-$3.distmatrix"
        date; cat "$2-$3.distmatrix"  | "${ROOTDIR}"/bin/linsetcover.out                        > "$2-$3.ordsetcover"
        date; cat "$2-$3.ordsetcover" | "${ROOTDIR}"/bin/setcover-ords-to-hdrs.out $1.fna       > "$2-$3.hdrsetcover-clu.tsv"
        echo "run_minenuc_with_infastafile_seqid($1, $2, $3) ended at $(date)"
    } | tee "$2-$3.minenuc.time"
}

function run_vsearch_with_infastafile_seqid() {
    time -p {
        echo "run_vsearch_with_infastafile_seqid($1, $2, $3) began at $(date)"
        date; "${ROOTDIR}"/benchmark/bin/vsearch --cluster_smallmem $1.fna --id 0.$3 -uc "$2_vsearch-$3.uc"
        date; cat "$2_vsearch-$3.uc" | "${ROOTDIR}"/benchmark/src/uc-to-clu-tsv.py > "$2_vsearch-$3-clu.tsv"
        echo "run_vsearch_with_infastafile_seqid($1, $2, $3) ended at $(date)"
    } | tee "$2-$3.vsearch.time"
}

function run_cdhitest_with_infastafile_seqid() {
    if [ $3 -lt 80 ]; then
        rm "$2-$3.cdhitest-clu.tsv" || true ; touch "$2-$3.cdhitest-clu.tsv"
        return 0
    fi
    cdhitestwordsize=10
    time -p {
        echo "run_cdhitest_with_infastafile_seqid($1, $2, $3) began at $(date)"
        date; "${ROOTDIR}"/benchmark/bin/cd-hit-est -i "$1.fna" -o "$2_cdhitest-M0-T0-d0-s80-c$3-n$cdhitestwordsize.fna" -M 0 -T 0 -d 0 -s 0.8 -c 0.$3 -n $cdhitestwordsize
        date; cat "$2_cdhitest-M0-T0-d0-s80-c$3-n$cdhitestwordsize.fna.clstr" | "${ROOTDIR}"/benchmark/src/clstr-to-clu-tsv.py > "$2-$3.cdhitest-clu.tsv"
        echo "run_cdhitest_with_infastafile_seqid($1, $2, $3) ended at $(date)"
    } | tee "$2-$3.cdhitest.time"
}

for SIM in $(echo $CSVSIM | sed "s/,/ /g"); do
    
    function gen_fam_metrics() {
        cat "$1.tsv" | "${ROOTDIR}"/benchmark/src/pfam-clstr-to-metrics.py | tee "$1.fam-metrics"
    }

    ## run Rfam.seed

    if [[ $4 == *"mine"* ]];    then run_minenuc_with_infastafile_seqid  "${INPREF}/Rfam.seed_shuf" "${OUTDIR}/Rfam.seed_shuf" $SIM 3334; fi
    if [[ $4 == *"vsearch"* ]]; then run_vsearch_with_infastafile_seqid  "${INPREF}/Rfam.seed_shuf" "${OUTDIR}/Rfam.seed_shuf" $SIM 3334; fi
    if [[ $4 == *"cdhit"* ]];   then run_cdhitest_with_infastafile_seqid "${INPREF}/Rfam.seed_shuf" "${OUTDIR}/Rfam.seed_shuf" $SIM 3334; fi

    if [[ $4 == *"mine"* ]];    then gen_fam_metrics "${OUTDIR}/Rfam.seed_shuf-${SIM}.hdrsetcover-clu"; fi
    if [[ $4 == *"vsearch"* ]]; then gen_fam_metrics "${OUTDIR}/Rfam.seed_shuf_vsearch-${SIM}-clu"    ; fi
    if [[ $4 == *"cdhit"* ]];   then gen_fam_metrics "${OUTDIR}/Rfam.seed_shuf-${SIM}.cdhitest-clu"   ; fi

    ## run Pfam-A.seed

    if [[ $4 == *"mine"* ]];     then run_mine_with_infastafile_seqid     "${INPREF}/Pfam-A.seed_shuf" "${OUTDIR}/Pfam-A.seed_shuf" $SIM 3334; fi
    if [[ $4 == *"linclust"* ]]; then run_linclust_with_infastafile_seqid "${INPREF}/Pfam-A.seed_shuf" "${OUTDIR}/Pfam-A.seed_shuf" $SIM 3334; fi
    if [[ $4 == *"quaclust"* ]]; then run_quaclust_with_infastafile_seqid "${INPREF}/Pfam-A.seed_shuf" "${OUTDIR}/Pfam-A.seed_shuf" $SIM 3334; fi
    if [[ $4 == *"cdhit"* ]];    then run_cdhit_with_infastafile_seqid    "${INPREF}/Pfam-A.seed_shuf" "${OUTDIR}/Pfam-A.seed_shuf" $SIM 3334; fi
    if [[ $4 == *"kclust"* ]];   then run_kclust_with_infastafile_seqid   "${INPREF}/Pfam-A.seed_shuf" "${OUTDIR}/Pfam-A.seed_shuf" $SIM 3334; fi
    
    if [[ $4 == *"mine"* ]];     then gen_fam_metrics "${OUTDIR}/Pfam-A.seed_shuf-${SIM}.hdrsetcover-clu"; fi
    if [[ $4 == *"linclust"* ]]; then gen_fam_metrics "${OUTDIR}/Pfam-A.seed_shuf-${SIM}.linclust-clu"   ; fi
    if [[ $4 == *"quaclust"* ]]; then gen_fam_metrics "${OUTDIR}/Pfam-A.seed_shuf-${SIM}.quaclust-clu"   ; fi
    if [[ $4 == *"cdhit"* ]];    then gen_fam_metrics "${OUTDIR}/Pfam-A.seed_shuf-${SIM}.cdhit-clu"      ; fi
    if [[ $4 == *"kclust"* ]];   then gen_fam_metrics "${OUTDIR}/Pfam-A.seed_shuf-${SIM}.kclust-clu"     ; fi
    
    # skip the rest if sim is either 60 or 80 
    if [ "60,80" == *"$SIM"* ]; then continue; fi

    ## run pdb
    
    if [[ $4 == *"mine"* ]];     then run_mine_with_infastafile_seqid     "${INPREF}/pdbent-seqres_shuf" "${OUTDIR}/pdbent-seqres_shuf" $SIM 3334; fi
    if [[ $4 == *"linclust"* ]]; then run_linclust_with_infastafile_seqid "${INPREF}/pdbent-seqres_shuf" "${OUTDIR}/pdbent-seqres_shuf" $SIM 3334; fi
    if [[ $4 == *"quaclust"* ]]; then run_quaclust_with_infastafile_seqid "${INPREF}/pdbent-seqres_shuf" "${OUTDIR}/pdbent-seqres_shuf" $SIM 3334; fi
    if [[ $4 == *"cdhit"* ]];    then run_cdhit_with_infastafile_seqid    "${INPREF}/pdbent-seqres_shuf" "${OUTDIR}/pdbent-seqres_shuf" $SIM 3334; fi
    if [[ $4 == *"kclust"* ]];   then run_kclust_with_infastafile_seqid   "${INPREF}/pdbent-seqres_shuf" "${OUTDIR}/pdbent-seqres_shuf" $SIM 3334; fi
    
    (
        flock -e 200
        cat "${OUTDIR}/pdbent-seqres_shuf-${SIM}."*clu.tsv | "${ROOTDIR}"/benchmark/src/update_memtable.py "${OUTDIR}/pdbent-seqres_shuf_memtable.tsv"
        function clu_tsv_to_tms() {
            cat "$1" | "${ROOTDIR}"/benchmark/src/setcover-hdrs-to-tms.py "${OUTDIR}/pdbent-seqres_shuf_memtable.tsv" > "$1".tms;
            benchmark/src/gen-pdbent-cdf.py "$1".tms | tee "$1".cdf;
        }
        if [[ $4 == *"mine"* ]];     then clu_tsv_to_tms "${OUTDIR}/pdbent-seqres_shuf-${SIM}.hdrsetcover-clu.tsv" ; fi
        if [[ $4 == *"linclust"* ]]; then clu_tsv_to_tms "${OUTDIR}/pdbent-seqres_shuf-${SIM}.linclust-clu.tsv"    ; fi
        if [[ $4 == *"quaclust"* ]]; then clu_tsv_to_tms "${OUTDIR}/pdbent-seqres_shuf-${SIM}.quaclust-clu.tsv"    ; fi
        if [[ $4 == *"cdhit"* ]];    then clu_tsv_to_tms "${OUTDIR}/pdbent-seqres_shuf-${SIM}.cdhit-clu.tsv"       ; fi
        if [[ $4 == *"kclust"* ]];   then clu_tsv_to_tms "${OUTDIR}/pdbent-seqres_shuf-${SIM}.kclust-clu.tsv"      ; fi 
    ) 200> "${OUTDIR}/pdbent-seqres_shuf_memtable.flockfile"

    # cat "${OUTDIR}/pdbent-seqres_shuf-$1.hdrsetcover-clu.tms" | grep "^tmscore" | awk '{print substr($2, 0, 4)}' | sort | uniq -c | awk '{print $2"\t"$1}'
    
    ## run uniprot

    timelimit=20000
    cdhit_timeout=false
    quaclust_timeout=false
    kclust_timeout=false
    #for uniref in down32x-uniref100-2017-01 down16x-uniref100-2017-01 down8x-uniref100-2017-01 down4x-uniref100-2017-01 down2x-uniref100-2017-01 uniref100-2017-01; do
    for uniref in uniref100-02 uniref100-12 uniref100-2011-01 uniref100-2014-01 uniref100-2017-01; do
        
        if [[ $4 == *"mine"* ]];     then run_mine_with_infastafile_seqid     "${INPREF}/${uniref}_shuf" "${OUTDIR}/${uniref}_shuf" $SIM $timelimit; fi
        if [[ $4 == *"linclust"* ]]; then run_linclust_with_infastafile_seqid "${INPREF}/${uniref}_shuf" "${OUTDIR}/${uniref}_shuf" $SIM $timelimit; fi
        if [[ $4 == *"quaclust"* ]]; then run_quaclust_with_infastafile_seqid "${INPREF}/${uniref}_shuf" "${OUTDIR}/${uniref}_shuf" $SIM $timelimit; fi
        if [[ $4 == *"cdhit"* ]];    then run_cdhit_with_infastafile_seqid    "${INPREF}/${uniref}_shuf" "${OUTDIR}/${uniref}_shuf" $SIM $timelimit; fi
        if [[ $4 == *"kclust"* ]];   then run_kclust_with_infastafile_seqid   "${INPREF}/${uniref}_shuf" "${OUTDIR}/${uniref}_shuf" $SIM $timelimit; fi
        timelimit=$(($timelimit*2))
    done
done

