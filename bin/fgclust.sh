#!/usr/bin/env sh

if [ $# -lt 2 ]; then
    echo "Usage: $0 <input-fasta-file> <output-fasta-file> <optional-parameters>"
    echo "  <input-fasta-file>  : the input fasta, fastq, or gzipped-fastq file that contains either only nucleotide sequences or only protein sequences"
    echo "  <output-fasta-file> : the output fasta file that contains centroid sequences"
    echo "The <optional-parameters> for sorting"
    echo "  --israndom          : sort by decreasing sequence length instead of shuffling randomly"
    echo "Note: "
    echo "  the alphabet (which can be of either amino-acid or nucleotide) of the sequences in input-fasta-file is automatically detected"
    echo "  <output-fasta-file>.clu is the output tsv file that contains the three following fields per line: "
    echo "                          the representative centroid, the member covered by the centroid, and the sequence similarity between the centroid and the member"
    echo "  if <output-fasta-file> ends with the extension .keep, then all intermediate files are kept instead of bein removed by default"
    echo "The <optional-parameters> for clustering:"
    echo "" | "${BINROOT}/fastaseqs-to-distmatrix.out" --help
    exit 0
fi

set -evx

BINROOT=$(dirname `which $0`)
INFILE="$1"
OUTFASTA="$2"
OUTCLUST="${2}.clu"

OUTDIR="${OUTFASTA}_tmpdir"

mkdir -p "${OUTDIR}/bin/" || true
cp "${BINROOT}/"*.out "${OUTDIR}/bin/"
BINROOT="${OUTDIR}/bin/"

if [[ "$INFILE" == *".gz" ]]; then catcmd=zcat; else catcmd=cat; fi

$catcmd "${INFILE}"                 | "${BINROOT}/len-revname-sort.out"        "${@:3}"              > "${OUTDIR}/shuf.fasta"
cat     "${OUTDIR}/shuf.fasta"      | "${BINROOT}/fastaseqs-to-distmatrix.out" "${@:3}"              > "${OUTDIR}/distmatrix.tsv"
cat     "${OUTDIR}/distmatrix.tsv"  | "${BINROOT}/linsetcover.out"                                   > "${OUTDIR}/ordsetcover.tsv"
cat     "${OUTDIR}/ordsetcover.tsv" | "${BINROOT}/setcover-ords-to-hdrs.out"  "${OUTDIR}/shuf.fasta" > "${OUTCLUST}"
cat     "${OUTDIR}/ordsetcover.tsv" | "${BINROOT}/setcover-ords-to-fasta.out" "${OUTDIR}/shuf.fasta" > "${OUTFASTA}"           

if [[ "${OUTFASTA}" != *".keep" ]]; then
    rm -r "${OUTDIR}"
fi

