#!/bin/sh
# A sample shell script to perform filtering
#
echo "Reading refinery config..." >&2
source refinery.properties
#
# parse command line args, the search path, relative to the current one is accepted
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
# uncoment below to clean residues (if re-running the filter for example)
# rm $SEARCHPATH/filters/*
#
# 0.3. make data and llists (adds the sample name to FASTA sequence names)
#
echo "  making raw data and llist..." >&2
cat ${SEARCHPATH}/transcripts_fpkm_3.fa | sed "s/>/>${SEARCHPATH}_/g" >${SEARCHPATH}/filters/raw_contigs.fa
makellist ${SEARCHPATH}/filters/raw_contigs.fa >${SEARCHPATH}/filters/raw_contigs.llist
#
# 0.4. run blat using the CDS data
#
echo "  running BLAT on ENSEMBL CDS..." >&2
blat ${reference_cds} ${SEARCHPATH}/filters/raw_contigs.fa ${SEARCHPATH}/filters/blat_raw_cds.psl
cat ${SEARCHPATH}/filters/blat_raw_cds.psl | ${refinery_blatparser} > ${SEARCHPATH}/filters/blat_raw_cds.best.tsv
#
# 0.5. re-filter by the CDS template
#
echo "  filtering raw contigs using BLAT CDS alignment..." >&2
Rscript ${refinery_bin}/filter_square75.R ${SEARCHPATH}/filters/raw_contigs.llist ${SEARCHPATH}/filters/blat_raw_cds.best.tsv ${SEARCHPATH}/filters/filter0_cds_keepers.list
include_mf ${SEARCHPATH}/filters/raw_contigs.fa ${SEARCHPATH}/filters/contigs_after_cds.fa ${SEARCHPATH}/filters/filter0_cds_keepers.list
makellist ${SEARCHPATH}/filters/contigs_after_cds.fa >${SEARCHPATH}/filters/contigs_after_cds.llist
#
# 0.4. run blat using the CDNA data
#
echo "  running BLAT on ENSEMBL CDNA..." >&2
blat ${reference_cdna} ${SEARCHPATH}/filters/contigs_after_cds.fa ${SEARCHPATH}/filters/blat_cds_cdna.psl
cat ${SEARCHPATH}/filters/blat_cds_cdna.psl | ${refinery_blatparser} > ${SEARCHPATH}/filters/blat_cds_cdna.best.tsv
#
# 0.5. re-filter by the CDNA template
#
echo "  filtering contigs using BLAT CDNA alignment..." >&2
Rscript ${refinery_bin}/filter_square75.R ${SEARCHPATH}/filters/contigs_after_cds.llist ${SEARCHPATH}/filters/blat_cds_cdna.best.tsv ${SEARCHPATH}/filters/filter1_cdna_keepers.list
include_mf ${SEARCHPATH}/filters/raw_contigs.fa ${SEARCHPATH}/filters/contigs_after_cdna.fa ${SEARCHPATH}/filters/filter1_cdna_keepers.list
makellist ${SEARCHPATH}/filters/contigs_after_cdna.fa >${SEARCHPATH}/filters/contigs_after_cdna.llist
#
# 0.6. run blat using the REFSEQ PEP data
#
echo "  running BLAT on REFSEQ PEP..." >&2
blat -t=dnax ${SEARCHPATH}/filters/contigs_after_cdna.fa -q=prot ${reference_pep_sequence} ${SEARCHPATH}/filters/blat_cdna_refseq_pep.psl
cat ${SEARCHPATH}/filters/blat_cdna_refseq_pep.psl | ${refinery_blatparser} -f t > ${SEARCHPATH}/filters/blat_cdna_refseq_pep.best.tsv
#
#blastall -a 3 -e .01 -d ${reference_pep} -p blastx -i ${SEARCHPATH}/filters/contigs_after_cdna.fa | gzip -9 - > ${SEARCHPATH}/filters/contigs_after_refseq_pep.blastx.out.gz
#zcat ${SEARCHPATH}/filters/contigs_after_refseq_pep.blastx.out.gz | ${refinery_bin}/blast_to_table.pl | ${refinery_bin}/hit_table_sorter.pl > ${SEARCHPATH}/filters/contigs_after_refseq_pep.best.tsv
#
# 0.7. re-filter by the REFSEQ PEP data
#
echo "  filtering contigs using BLAT REFSEQ PEP alignment..." >&2
Rscript ${refinery_bin}/ffilter_square75_iverse.R ${SEARCHPATH}/filters/contigs_after_cdna.llist ${SEARCHPATH}/filters/blat_cdna_refseq_pep.best.tsv ${SEARCHPATH}/filters/filter2_refseq_pep_keepers.list
include_mf ${SEARCHPATH}/filters/raw_contigs.fa ${SEARCHPATH}/filters/contigs_after_refseq_pep.fa ${SEARCHPATH}/filters/filter2_refseq_pep_keepers.list
makellist ${SEARCHPATH}/filters/contigs_after_refseq_pep.fa >${SEARCHPATH}/filters/contigs_after_refseq_pep.llist