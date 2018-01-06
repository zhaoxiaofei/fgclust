#!/usr/bin/env sh

BINROOT=$(dirname `which $0`)
if [ $# -lt 2 ]; then
    echo "Usage: $0 <input-fasta-file> <output-fasta-file> <optional-parameters>"
    echo "  <input-fasta-file>  : the input fasta, fastq, or gzipped-fastq file that contains either only nucleotide sequences or only protein sequences"
    echo "  <output-fasta-file> : the output fasta file that contains centroid sequences"
    echo "Note: "
    echo "  the sequence similarity cutoff (--sim-perc) is by default 50 for protein, 70 for RNA, and 90 for DNA."
    echo "  the type of input sequences (which can be of protein, RNA, or DNA) in <input-fasta-file> is by default automatically detected."
    echo "  <output-fasta-file>.clu is the output tsv file that contains the three following fields per line: "
    echo "                          the representative centroid, the member covered by the centroid, and the sequence similarity between the centroid and the member"
    echo "  if <output-fasta-file> ends with the extension .keep, then all intermediate files are kept instead of being removed by default"
    if [ "$1" == "-help" ]; then 
        echo "The <optional-parameters> for processing sequences."
        echo "" | "${BINROOT}/procseqs.out" --help
        echo "The <optional-parameters> for computing distance matrix:"
        echo "" | "${BINROOT}/fastaseqs-to-distmatrix.out" --help
        echo "The <optional-parameters> for clustering:"
        echo "" | "${BINROOT}/linsetcover.out" --help
    else
        echo "For more detail about <optional-parameters>, please enter"
        echo "  $0 -help"
    fi
    exit 0
fi

set -evx

INFILE="$1"
OUTFASTA="$2"
OUTCLUST="${2}.clu"

OUTDIR="${OUTFASTA}_tmpdir"

mkdir -p "${OUTDIR}/bin/" || true
cp "${BINROOT}/"*.out "${OUTDIR}/bin/"
BINROOT="${OUTDIR}/bin/"

if [[ "$INFILE" == *".gz" ]]; then catcmd=zcat; else catcmd=cat; fi

$catcmd "${INFILE}"                 | "${BINROOT}/procseqs.out"                "${@:3}"              > "${OUTDIR}/shuf.fasta"
cat     "${OUTDIR}/shuf.fasta"      | "${BINROOT}/fastaseqs-to-distmatrix.out" "${@:3}"              > "${OUTDIR}/distmatrix.tsv"
cat     "${OUTDIR}/distmatrix.tsv"  | "${BINROOT}/linsetcover.out"             "${@:3}"              > "${OUTDIR}/ordsetcover.tsv"
cat     "${OUTDIR}/ordsetcover.tsv" | "${BINROOT}/setcover-ords-to-hdrs.out"  "${OUTDIR}/shuf.fasta" > "${OUTCLUST}"
cat     "${OUTDIR}/ordsetcover.tsv" | "${BINROOT}/setcover-ords-to-fasta.out" "${OUTDIR}/shuf.fasta" > "${OUTFASTA}"           

if [[ "${OUTFASTA}" != *".keep" ]]; then
    rm -r "${OUTDIR}"
fi

