#!/usr/bin/env sh

BINROOT=$(dirname `which $0`)
if [ $# -lt 3 ]; then
    echo "Usage: $0 <input-fasta-file> <output-fasta-file> <sequence-similarity> <optional-params>"
    echo "  input-fasta-file    : the input fasta, fastq, or gzipped-fastq file that contains either only nucleotide sequences or only protein sequences"
    echo "  output-fasta-file   : the output fasta file that contains centroid sequences"
    echo "  sequence-similarity : the percent sequence similarity between 1 and 99"
    echo "Optional parameters for sorting"
    echo "  --israndom          : sort by decreasing sequence length instead of shuffling randomly"
    echo "Note: "
    echo "  the alphabet (which can be of either amino-acid or nucleotide) of the sequences in input-fasta-file is automatically detected"
    echo "  <output-fasta-file>.clu is the output tsv file that contains the three following fields per line: "
    echo "                          the representative centroid, the member covered by the centroid, and the sequence similarity between the centroid and the member"
    echo "Optional parameters for clustering:"
    "${BINROOT}/fastaseqs-to-distmatrix.out" --help
    exit 0
fi

set -evx

INFILE="$1"
OUTFASTA="$2"
SIM="$3"

OUTDIR="${OUTFASTA}_tmpdir"

mkdir -p "${OUTDIR}" || true
cp "${BINROOT}/"*.out "${OUTDIR}"

if [[ "$INFILE" == *.gz ]]; then catcmd=zcat; else catcmd=cat; fi

$catcmd "${INFILE}"                 | "${OUTDIR}/len-revname-sort.out"        "${@:4}"                     > "${OUTDIR}/shuf.fasta"
cat     "${OUTDIR}/shuf.fasta"      | "${OUTDIR}/fastaseqs-to-distmatrix.out" --sim-perc "${SIM}" "${@:4}" > "${OUTDIR}/distmatrix.tsv"
cat     "${OUTDIR}/distmatrix.tsv"  | "${OUTDIR}/linsetcover.out"                                          > "${OUTDIR}/ordsetcover.tsv"
cat     "${OUTDIR}/ordsetcover.tsv" | "${OUTDIR}/setcover-ords-to-hdrs.out"  "${OUTDIR}/shuf.fasta"        > "${OUTFASTA}.clu"
cat     "${OUTDIR}/ordsetcover.tsv" | "${OUTDIR}/setcover-ords-to-fasta.out" "${OUTDIR}/shuf.fasta"        > "${OUTFASTA}"           

rm -r "${OUTDIR}"

