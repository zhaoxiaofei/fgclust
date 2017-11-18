#!/usr/bin/env sh
#set -evx

BINROOT=$(dirname `which $0`)
if [ $# -lt 4 ]; then
    echo "Usage: $0 <input-fasta-file> <output-fasta-file> <output-cluster-file> <sequence-similarity> <optional-params>"
    echo "  input-fasta-file    : the input fasta, fastq, or gzipped-fastq file that contains either only nucleotide sequences or only protein sequences"
    echo "  output-fasta-file   : the output fasta file that contains centroid sequences"
    echo "  output-cluster-file : the output tsv file that contains the three following fields per line: "
    echo "                        the representative centroid, the member covered by the centroid, and the sequence similarity between the centroid and the member"
    echo "  sequence-similarity : the percent sequence similarity between 1 and 99"
    echo "Optional parameters for sorting"
    echo "  --israndom          : sort by decreasing sequence length instead of shuffling randomly"
    echo "Note: "
    echo "  the alphabet (which can be of either amino-acid or nucleotide) of the sequences in input-fasta-file is automatically detected"
    echo "Optional parameters for clustering:"
    "${BINROOT}/fastaseqs-to-distmatrix.out" --help
    exit 0
fi

INFILE="$1"
OUTFASTA="$2"
OUTCLUST="$3"
SIM="$4"

OUTDIR="${OUTFASTA}_tmpdir"

mkdir -p "${OUTDIR}" || true

if [[ "$INFILE" == *.gz ]]; then catcmd=zcat; else catcmd=cat; fi

$catcmd "${INFILE}"                 | "${BINROOT}/len-revname-sort.out"        "${@:5}"                     > "${OUTDIR}/shuf.fasta"
cat     "${OUTDIR}/shuf.fasta"      | "${BINROOT}/fastaseqs-to-distmatrix.out" --sim-perc "${SIM}" "${@:5}" > "${OUTDIR}/distmatrix.tsv"
cat     "${OUTDIR}/distmatrix.tsv"  | "${BINROOT}/linsetcover.out"                                          > "${OUTDIR}/ordsetcover.tsv"
cat     "${OUTDIR}/ordsetcover.tsv" | "${BINROOT}/setcover-ords-to-hdrs.out"  "${OUTDIR}/shuf.fasta"        > "${OUTCLUST}"
cat     "${OUTDIR}/ordsetcover.tsv" | "${BINROOT}/setcover-ords-to-fasta.out" "${OUTDIR}/shuf.fasta"        > "${OUTFASTA}"           

rm -r "${OUTDIR}"

