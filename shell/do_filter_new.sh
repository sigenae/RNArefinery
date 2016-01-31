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
# 0.3. run blat using the data
#
echo "  running BLAT on NEW MRNA SET..." >&2
blat ${reference5_mrna} ${SEARCHPATH}/filters/raw_contigs.fa ${SEARCHPATH}/filters/blat_raw_genome5_mrna.psl
cat ${SEARCHPATH}/filters/blat_raw_genome5_mrna.psl | ${refinery_blatparser} > ${SEARCHPATH}/filters/blat_raw_genome5_mrna.best.tsv
#
# 0.4 filter those with less than 20% overlap on CDS
#
echo "  filtering raw contigs using NEW MRNA SET alignment..." >&2
Rscript ${refinery_bin}/filter_square75.R ${SEARCHPATH}/filters/raw_contigs.llist ${SEARCHPATH}/filters/blat_raw_genome5_mrna.best.tsv ${SEARCHPATH}/filters/genome5_mrna_keepers.list
include_mf ${SEARCHPATH}/filters/raw_contigs.fa ${SEARCHPATH}/filters/contigs_after_genome5_mrna.fa ${SEARCHPATH}/filters/genome5_mrna_keepers.list
makellist ${SEARCHPATH}/filters/contigs_after_genome5_mrna.fa >${SEARCHPATH}/filters/contigs_after_genome5_mrna.llist
