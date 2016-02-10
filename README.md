### RNA Refinery

# A wrapper pipeline designed to find un-annotated transcripts. The main [R](https://cran.r-project.org) script runs:
1. [SRADb](https://www.bioconductor.org/packages/release/bioc/html/SRAdb.html)-backed search identifying publicly available datasets for an organism of interest.
2. [Aspera/FTP](http://www.ncbi.nlm.nih.gov/books/NBK242625/)-based retrieval of the data.
3. deNovo assembly of transcriptome from SRA runs using [DRAP](http://www.sigenae.org/drap/index.html) (an Oases/Trinity-based assembler).
4. Multi-step filtering process based on publicly available resources (ENSEMBL & NCBI -backed cDNA/CDS/mRNA/protein/ncRNA datasets).
5. Blast-based annotation of sequences made it through those filters.
