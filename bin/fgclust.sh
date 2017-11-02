#!/usr/bin/env sh
set -evx
if [ $# -lt 4 ]; then
    echo "Usage: $0 <input-fasta-file> <output-fasta-file> <output-cluster-file> <sequence-similarity>"
    echo "  input-fasta-file    : the input fasta file that contains either only nucleotide sequences or only protein sequences"
    echo "  output-fasta-file   : the output fasta file that contains centroid sequences"
    echo "  output-cluster-file : the output tsv file that contains the three following fields per line: "
    echo "                        the representative centroid, the member covered by the centroid, and the sequence similarity between the centroid and the member"
    echo "  sequence-similarity : the percent sequence similarity between 1 and 99"
    echo "Note: "
    echo "  the alphabet (which can be of either amino-acid or nucleotide) of the sequences in input-fasta-file is automatically detected"
    exit 0
fi

INFILE="$1"
OUTFASTA="$2"
OUTCLUST="$3"
SIM="$4"

OUTDIR="${OUTFASTA}_tmpdir"

BINROOT=$(dirname `which $0`)

mkdir -p "${OUTDIR}" || true

cat "${INFILE}"                 | "${BINROOT}/len-revname-sort.out" --israndom 1                 > "${OUTDIR}/shuf.fasta"
cat "${OUTDIR}/shuf.fasta"      | "${BINROOT}/fastaseqs-to-distmatrix.out" --sim-perc "${SIM}"   > "${OUTDIR}/distmatrix.tsv"
cat "${OUTDIR}/distmatrix.tsv"  | "${BINROOT}/linsetcover.out"                                   > "${OUTDIR}/ordsetcover.tsv"
cat "${OUTDIR}/ordsetcover.tsv" | "${BINROOT}/setcover-ords-to-hdrs.out"  "${OUTDIR}/shuf.fasta" > "${OUTCLUST}"
cat "${OUTDIR}/ordsetcover.tsv" | "${BINROOT}/setcover-ords-to-fasta.out" "${OUTDIR}/shuf.fasta" > "${OUTFASTA}"           

rm -r "${OUTDIR}"

