#!/usr/bin/env sh
set -evx

if [ -n "$SGE_TASK_ID" ]; then
    INPREF="${qsubINPREF}"
    OUTDIR="${qsubOUTDIR}"
    CSVSIM="${qsubCSVSIM}"
    PROG="${qsubPROG}" #$4
else
    ROOTDIR=$(dirname `which $0`)
    INPREF="$1"
    OUTDIR="$2"
    CSVSIM="$3"
    PROG="$4"
fi

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
    { time -p {
        echo "run_mine_with_infastafile_seqid($1, $2, $3, $4) eval-began-at $(date)"
        date; cat "$1.faa" "$1.fna"   | "${FGCLUST}"/fastaseqs-to-distmatrix.out --sim-perc $3 > "$2-$3.distmatrix" || true
        date; cat "$2-$3.distmatrix"  | "${FGCLUST}"/linsetcover.out                           > "$2-$3.ordsetcover"
        date; cat "$2-$3.ordsetcover" | "${FGCLUST}"/setcover-ords-to-hdrs.out $1.faa          > "$2-$3.hdrsetcover-clu.tsv"
        echo "run_mine_with_infastafile_seqid($1, $2, $3) eval-ended-at $(date)"
    } } 2>&1 | tee "$2-$3.mine.time"
}

function run_quaclust_with_infastafile_seqid() {
    { time -p {
        echo "run_mmseqsclust_with_infastafile_seqid($1, $2, $3) eval-began-at $(date)"
        rm -r "$2-$3.quaclust-tmpdir" || true
        mkdir "$2-$3.quaclust-tmpdir"
        rm "$2-$3.quaclust-clu" || true
        rm "$2-$3.quaclust-clu.tsv" || true
        set +e
        date; timeout $4 "${ROOTDIR}"/benchmark/bin/mmseqs cluster "$1.mmseqs-db" "$2-$3.quaclust-clu" "$2-$3.quaclust-tmpdir" --min-seq-id 0.$3
        resetfile "$2-$3.quaclust-clu.tsv"
        if "${ISTIMEOUT}"; then return 0 ; fi
        date; "${ROOTDIR}"/benchmark/bin/mmseqs createtsv             "$1.mmseqs-db" "$1.mmseqs-db"       "$2-$3.quaclust-clu" "$2-$3.quaclust-clu.tsv"
        echo "run_mmseqsclust_with_infastafile_seqid($1, $2, $3) eval-ended-at $(date)"
    } } 2>&1 | tee "$2-$3.quaclust.time"
}

function run_linclust_with_infastafile_seqid() {
    { time -p {
        echo "run_linclust_with_infastafile_seqid($1, $2, $3) eval-began-at $(date)"
        rm -r "$2-$3.linclust-tmpdir" || true
        mkdir "$2-$3.linclust-tmpdir"
        rm "$2-$3.linclust-clu" || true
        rm "$2-$3.linclust-clu.tsv" || true
        date; "${ROOTDIR}"/benchmark/bin/mmseqs linclust         "$1.mmseqs-db" "$2-$3.linclust-clu" "$2-$3.linclust-tmpdir" --min-seq-id 0.$3
        date; "${ROOTDIR}"/benchmark/bin/mmseqs createtsv        "$1.mmseqs-db" "$1.mmseqs-db" "$2-$3.linclust-clu" "$2-$3.linclust-clu.tsv"
        echo "run_linclust_with_infastafile_seqid($1, $2, $3) eval-ended-at $(date)"
    } } 2>&1 | tee "$2-$3.linclust.time"
}

function run_kclust_with_infastafile_seqid() {
    { time -p {
        echo "run_kclust_with_infastafile_seqid($1, $2, $3) eval-began-at $(date)"
        set +e
        date; timeout $4 "${ROOTDIR}"/benchmark/bin/kClust -M 111000 -i "$1.faa" -d "$2-$3.kclust" -s $(python -c "print(round(-0.68506329113924050632 + 0.06025316455696202531*$3, 2))")
        resetfile "$2-$3.kclust-clu.tsv"
        if "${ISTIMEOUT}"; then return 0 ; fi
        date; cat "$2-$3.kclust/clusters.dmp" | tail -n +2 | awk '{print $2"\t"$1"\t"1}' | "${ROOTDIR}"/bin/setcover-ords-to-hdrs.out "$2-$3.kclust/db_sorted.fas" > "$2-$3.kclust-clu.tsv"
        echo "run_kclust_with_infastafile_seqid($1, $2, $3) eval-ended-at $(date)"
    } } 2>&1 | tee "$2-$3.kclust.time"
    if grep -Fq " eval-ended-at " "$2-$3.kclust.time"; then kclustnext=true; else kclustnext=false; fi
}

function run_cdhit_with_infastafile_seqid() {
    cdhitwordsize=5
    if [ $SIM -lt 70 ] ; then cdhitwordsize=4; fi
    if [ $SIM -lt 60 ] ; then cdhitwordsize=3; fi
    { time -p {
        echo "run_cdhit_with_infastafile_seqid($1, $2, $3) eval-began-at $(date)"
        set +e
        date; timeout $4 "${ROOTDIR}"/benchmark/bin/cd-hit -i "$1.faa" -o "$2_cdhit-M0-T0-d0-s80-c$3-n$cdhitwordsize.faa" -M 0 -T 0 -d 0 -s 0.8 -c 0.$3 -n $cdhitwordsize
        resetfile "$2-$3.cdhit-clu.tsv"
        if "${ISTIMEOUT}"; then return 0 ; fi
        date; cat "$2_cdhit-M0-T0-d0-s80-c$3-n$cdhitwordsize.faa.clstr" | "${ROOTDIR}"/benchmark/src/clstr-to-clu-tsv.py > "$2-$3.cdhit-clu.tsv"
        echo "run_cdhit_with_infastafile_seqid($1, $2, $3) eval-ended-at $(date)"
    } } 2>&1 | tee "$2-$3.cdhit.time" 
    if grep -Fq " eval-ended-at " "$2-$3.cdhit.time"; then cdhitnext=true; else cdhitnext=false; fi
}

function run_vsearch_with_infastafile_seqid() {
    { time -p {
        echo "run_vsearch_with_infastafile_seqid($1, $2, $3) eval-began-at $(date)"
        date; "${ROOTDIR}"/benchmark/bin/vsearch --cluster_smallmem $1.fna --id 0.$3 -uc "$2_vsearch-$3.uc"
        date; cat "$2_vsearch-$3.uc" | "${ROOTDIR}"/benchmark/src/uc-to-clu-tsv.py > "$2_vsearch-$3-clu.tsv"
        echo "run_vsearch_with_infastafile_seqid($1, $2, $3) eval-ended-at $(date)"
    } } 2>&1 | tee "$2-$3.vsearch.time"
}

function run_cdhitest_with_infastafile_seqid() {
    if [ $3 -lt 80 ]; then
        rm "$2-$3.cdhitest-clu.tsv" || true ; touch "$2-$3.cdhitest-clu.tsv"
        return 0
    fi
    cdhitestwordsize=10
    { time -p {
        echo "run_cdhitest_with_infastafile_seqid($1, $2, $3) eval-began-at $(date)"
        date; "${ROOTDIR}"/benchmark/bin/cd-hit-est -i "$1.fna" -o "$2_cdhitest-M0-T0-d0-s80-c$3-n$cdhitestwordsize.fna" -M 0 -T 0 -d 0 -s 0.8 -c 0.$3 -n $cdhitestwordsize
        date; cat "$2_cdhitest-M0-T0-d0-s80-c$3-n$cdhitestwordsize.fna.clstr" | "${ROOTDIR}"/benchmark/src/clstr-to-clu-tsv.py > "$2-$3.cdhitest-clu.tsv"
        echo "run_cdhitest_with_infastafile_seqid($1, $2, $3) eval-ended-at $(date)"
    } } 2>&1 | tee "$2-$3.cdhitest.time"
}

for SIM in $(echo $CSVSIM | sed "s/,/ /g"); do
    FGCLUST="${OUTDIR}/${SIM}-fgclust-bin/"
    if [[ "${PROG}" == *"mine"* ]]; then mkdir -p "${FGCLUST}"; cp "${ROOTDIR}/bin/"*.out "${FGCLUST}"; fi
done

for SIM in $(echo $CSVSIM | sed "s/,/ /g"); do
    FGCLUST="${OUTDIR}/${SIM}-fgclust-bin/"
 
    function gen_fam_metrics() {
        cat "$1.tsv" | "${ROOTDIR}"/benchmark/src/pfam-clstr-to-metrics.py | tee "$1.fam-metrics"
    }

    ## run Rfam.seed

    if [[ "${PROG}" == *"mine"* ]];    then run_mine_with_infastafile_seqid     "${INPREF}/Rfam.seed_shuf" "${OUTDIR}/Rfam.seed_shuf" $SIM 3334; fi
    if [[ "${PROG}" == *"vsearch"* ]]; then run_vsearch_with_infastafile_seqid  "${INPREF}/Rfam.seed_shuf" "${OUTDIR}/Rfam.seed_shuf" $SIM 3334; fi
    if [[ "${PROG}" == *"cdhit"* ]];   then run_cdhitest_with_infastafile_seqid "${INPREF}/Rfam.seed_shuf" "${OUTDIR}/Rfam.seed_shuf" $SIM 3334; fi

    if [[ "${PROG}" == *"mine"* ]];    then gen_fam_metrics "${OUTDIR}/Rfam.seed_shuf-${SIM}.hdrsetcover-clu"; fi
    if [[ "${PROG}" == *"vsearch"* ]]; then gen_fam_metrics "${OUTDIR}/Rfam.seed_shuf_vsearch-${SIM}-clu"    ; fi
    if [[ "${PROG}" == *"cdhit"* ]];   then gen_fam_metrics "${OUTDIR}/Rfam.seed_shuf-${SIM}.cdhitest-clu"   ; fi

    ## run Pfam-A.seed

    if [[ "${PROG}" == *"mine"* ]];     then run_mine_with_infastafile_seqid     "${INPREF}/Pfam-A.seed_shuf" "${OUTDIR}/Pfam-A.seed_shuf" $SIM 3334; fi
    if [[ "${PROG}" == *"linclust"* ]]; then run_linclust_with_infastafile_seqid "${INPREF}/Pfam-A.seed_shuf" "${OUTDIR}/Pfam-A.seed_shuf" $SIM 3334; fi
    if [[ "${PROG}" == *"quaclust"* ]]; then run_quaclust_with_infastafile_seqid "${INPREF}/Pfam-A.seed_shuf" "${OUTDIR}/Pfam-A.seed_shuf" $SIM 3334; fi
    if [[ "${PROG}" == *"cdhit"* ]];    then run_cdhit_with_infastafile_seqid    "${INPREF}/Pfam-A.seed_shuf" "${OUTDIR}/Pfam-A.seed_shuf" $SIM 3334; fi
    if [[ "${PROG}" == *"kclust"* ]];   then run_kclust_with_infastafile_seqid   "${INPREF}/Pfam-A.seed_shuf" "${OUTDIR}/Pfam-A.seed_shuf" $SIM 3334; fi
    
    if [[ "${PROG}" == *"mine"* ]];     then gen_fam_metrics "${OUTDIR}/Pfam-A.seed_shuf-${SIM}.hdrsetcover-clu"; fi
    if [[ "${PROG}" == *"linclust"* ]]; then gen_fam_metrics "${OUTDIR}/Pfam-A.seed_shuf-${SIM}.linclust-clu"   ; fi
    if [[ "${PROG}" == *"quaclust"* ]]; then gen_fam_metrics "${OUTDIR}/Pfam-A.seed_shuf-${SIM}.quaclust-clu"   ; fi
    if [[ "${PROG}" == *"cdhit"* ]];    then gen_fam_metrics "${OUTDIR}/Pfam-A.seed_shuf-${SIM}.cdhit-clu"      ; fi
    if [[ "${PROG}" == *"kclust"* ]];   then gen_fam_metrics "${OUTDIR}/Pfam-A.seed_shuf-${SIM}.kclust-clu"     ; fi

done

for SIM in $(echo $CSVSIM | sed "s/,/ /g"); do
    FGCLUST="${OUTDIR}/${SIM}-fgclust-bin/"

    # skip the rest if sim is either 60 or 80 
    if [[ "60,80" == *"$SIM"* ]]; then continue; fi

    ## run pdb 
    (
        if [[ "${PROG}" == *"mine"* ]];     then run_mine_with_infastafile_seqid     "${INPREF}/pdbent-seqres_shuf" "${OUTDIR}/pdbent-seqres_shuf" $SIM 3334; fi
        if [[ "${PROG}" == *"linclust"* ]]; then run_linclust_with_infastafile_seqid "${INPREF}/pdbent-seqres_shuf" "${OUTDIR}/pdbent-seqres_shuf" $SIM 3334; fi
        if [[ "${PROG}" == *"quaclust"* ]]; then run_quaclust_with_infastafile_seqid "${INPREF}/pdbent-seqres_shuf" "${OUTDIR}/pdbent-seqres_shuf" $SIM 3334; fi
        if [[ "${PROG}" == *"cdhit"* ]];    then run_cdhit_with_infastafile_seqid    "${INPREF}/pdbent-seqres_shuf" "${OUTDIR}/pdbent-seqres_shuf" $SIM 3334; fi
        if [[ "${PROG}" == *"kclust"* ]];   then run_kclust_with_infastafile_seqid   "${INPREF}/pdbent-seqres_shuf" "${OUTDIR}/pdbent-seqres_shuf" $SIM 3334; fi
        flock -e 200
        cat "${OUTDIR}/pdbent-seqres_shuf-${SIM}."*clu.tsv | "${ROOTDIR}"/benchmark/src/update_memtable.py "${OUTDIR}/pdbent-seqres_shuf_memtable.tsv"
        function clu_tsv_to_tms() {
            cat "$1" | "${ROOTDIR}"/benchmark/src/setcover-hdrs-to-tms.py "${OUTDIR}/pdbent-seqres_shuf_memtable.tsv" > "$1".tms;
            benchmark/src/gen-pdbent-cdf.py "$1".tms | tee "$1".cdf;
        }
        if [[ "${PROG}" == *"mine"* ]];     then clu_tsv_to_tms "${OUTDIR}/pdbent-seqres_shuf-${SIM}.hdrsetcover-clu.tsv" ; fi
        if [[ "${PROG}" == *"linclust"* ]]; then clu_tsv_to_tms "${OUTDIR}/pdbent-seqres_shuf-${SIM}.linclust-clu.tsv"    ; fi
        if [[ "${PROG}" == *"quaclust"* ]]; then clu_tsv_to_tms "${OUTDIR}/pdbent-seqres_shuf-${SIM}.quaclust-clu.tsv"    ; fi
        if [[ "${PROG}" == *"cdhit"* ]];    then clu_tsv_to_tms "${OUTDIR}/pdbent-seqres_shuf-${SIM}.cdhit-clu.tsv"       ; fi
        if [[ "${PROG}" == *"kclust"* ]];   then clu_tsv_to_tms "${OUTDIR}/pdbent-seqres_shuf-${SIM}.kclust-clu.tsv"      ; fi 
    ) 200> "${OUTDIR}/pdbent-seqres_shuf_memtable.flockfile"

    ## run uniprot

    timelim=$((3600*50))
    if [[ "${PROG}" == *"mine"* ]]    ; then minenext=true    ; else minenext=false    ; fi 
    if [[ "${PROG}" == *"linclust"* ]]; then linclustnext=true; else linclustnext=false; fi 
    if [[ "${PROG}" == *"quaclust"* ]]; then quaclustnext=true; else quaclustnext=false; fi 
    if [[ "${PROG}" == *"cdhit"* ]]   ; then cdhitnext=true   ; else cdhitnext=false   ; fi 
    if [[ "${PROG}" == *"kclust"* ]]  ; then kclustnext=true  ; else kclustnext=false  ; fi 
    #for uniref in down32x-uniref100-2017-01 down16x-uniref100-2017-01 down8x-uniref100-2017-01 down4x-uniref100-2017-01 down2x-uniref100-2017-01 uniref100-2017-01; do
    for uniref in uniref100-02 uniref100-12 uniref100-2011-01 uniref100-2014-01 uniref100-2017-01; do    
        if $minenext    ; then run_mine_with_infastafile_seqid     "${INPREF}/${uniref}_shuf" "${OUTDIR}/${uniref}_shuf" $SIM $timelim; fi
        if $linclustnext; then run_linclust_with_infastafile_seqid "${INPREF}/${uniref}_shuf" "${OUTDIR}/${uniref}_shuf" $SIM $timelim; fi
        if $quaclustnext; then run_quaclust_with_infastafile_seqid "${INPREF}/${uniref}_shuf" "${OUTDIR}/${uniref}_shuf" $SIM $timelim; fi
        if $cdhitnext;    then run_cdhit_with_infastafile_seqid    "${INPREF}/${uniref}_shuf" "${OUTDIR}/${uniref}_shuf" $SIM $timelim; fi
        if $kclustnext;   then run_kclust_with_infastafile_seqid   "${INPREF}/${uniref}_shuf" "${OUTDIR}/${uniref}_shuf" $SIM $timelim; fi
    done
done

