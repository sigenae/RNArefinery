#!/bin/sh
# Sample shell script to perform filtering
echo "Reading refinery config..." >&2
source refinery.properties
#
# parse command line args
#
while [[ $# > 1 ]]
do
key="$1"
case $key in
    -s|--searchpath)
    SEARCHPATH="$2"
    shift # past argument
    ;;
    --default)
    DEFAULT=YES
    ;;
    *)
            # unknown option
    ;;
esac
shift # past argument or value
done
echo SEARCH PATH USED "$SEARCHPATH" >&2
#
# 0.1. making sure FPKM3 exists
#
if [ ! -f $SEARCHPATH/transcripts_fpkm_3.fa ]; then
    echo "File not found, breaking!"
    exit 10
fi
echo "  $SEARCHPATH/transcripts_fpkm_3.fa found..." >&2
#
# 0.2. make the filter folder
#
mkdir -p $SEARCHPATH/filters
#
# 0.3. make data and llists
#
echo "  making raw data and llist..." >&2
cat ${SEARCHPATH}/transcripts_fpkm_3.fa | sed "s/>/>${SEARCHPATH}_/g" >${SEARCHPATH}/filters/raw_contigs.fa
makellist ${SEARCHPATH}/filters/raw_contigs.fa >${SEARCHPATH}/filters/raw_contigs.llist
#
# 0.4. run blat using the MRNA5 data
#
#echo "  running BLAT on the NEW CDNA/MRNA? SET..." >&2
#blat ${reference5_mrna} ${SEARCHPATH}/filters/raw_contigs.fa ${SEARCHPATH}/filters/blat_cdna5.psl
#cat ${SEARCHPATH}/filters/blat_cdna5.psl | ${refinery_blatparser} > ${SEARCHPATH}/filters/blat_cdna5.best.tsv
#
# 0.5 filter those with less than 20% overlap on MRNA template
#
#echo "  filtering raw contigs using NEW MRNA SET alignment..." >&2
#Rscript ${refinery_bin}/filter_square75.R ${SEARCHPATH}/filters/raw_contigs.llist ${SEARCHPATH}/filters/blat_cdna5.best.tsv ${SEARCHPATH}/filters/blat_cdna5.keepers.list
#include_mf ${SEARCHPATH}/filters/raw_contigs.fa ${SEARCHPATH}/filters/contigs_after_cdna5.fa ${SEARCHPATH}/filters/blat_cdna5.keepers.list
#makellist ${SEARCHPATH}/filters/contigs_after_cdna5.fa >${SEARCHPATH}/filters/contigs_after_cdna5.llist
#
# 0.6. run blat using the PEP5 data
#
#echo "  running BLAT on NEW PEPTIDE SET..." >&2
#blat -t=dnax ${SEARCHPATH}/filters/contigs_after_cdna5.fa -q=prot ${reference5_protein} ${SEARCHPATH}/filters/blat_cdna5_pep5.psl
#cat ${SEARCHPATH}/filters/blat_cdna5_pep5.psl | ${refinery_blatparser} -f t > ${SEARCHPATH}/filters/blat_cdna5_pep5.best.tsv
#
# 0.7. re-filter by the REFSEQ PEP data
#
echo "  filtering contigs using BLAT REFSEQ PEP alignment..." >&2
Rscript ${refinery_bin}/filter_square75_iverse.R ${SEARCHPATH}/filters/contigs_after_cdna5.llist ${SEARCHPATH}/filters/blat_cdna5_pep5.best.tsv ${SEARCHPATH}/filters/blat_cdna5_pep5_keepers.list
include_mf ${SEARCHPATH}/filters/raw_contigs.fa ${SEARCHPATH}/filters/contigs_after_cdna5_pep5.fa ${SEARCHPATH}/filters/blat_cdna5_pep5_keepers.list
makellist ${SEARCHPATH}/filters/contigs_after_cdna5_pep5.fa >${SEARCHPATH}/filters/contigs_after_cdna5_pep5.llist
