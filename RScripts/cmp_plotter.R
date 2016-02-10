#
# Comparing RNA Refinery results when using two genome versions (V4 ENSEMBL and V5 NCBI).
#
require(stringr)
require(reshape)
require(plyr)
require(dplyr)
require(data.table)
#
require(RSQLite)
require(SRAdb)
require(RMySQL)
#
require(ggplot2)
require(gridExtra)
require(scales)
require(Cairo)
###################################################################################################
# 0.0 define the formatter for plots
#
percent <- function(x, digits = 2, format = "f", ...) {
  paste0(formatC(100 * x, format = format, digits = digits, ...), "%")
}
###################################################################################################
# 0.1 Set the working folder
#
archive_dir <- "/Users/psenin/RNARefinery"
setwd(archive_dir)
###################################################################################################
# FUNCTION: list all the folders within the directory
#
list_dirs <-
  function(path = ".", pattern = NULL, all.dirs = FALSE, full.names = FALSE, ignore.case = FALSE) {
    # use full.names=TRUE to pass to file.info
    all <-
      list.files(path, pattern, all.dirs, full.names = TRUE, recursive = FALSE, ignore.case)
    dirs <- all[file.info(all)$isdir]
    if (isTRUE(full.names))
      return(dirs)
    else
      return(basename(dirs))
  }
###################################################################################################
# 0.2 GENOME V4 folders
#
dirs4 <- data.frame(folder = paste(archive_dir, "/genome4/",
                                   list_dirs(paste(archive_dir, "/genome4", sep = "")), sep = ""))
dirs4$accession <- str_match(dirs4$folder, ".*/(.*)")[,2]
###################################################################################################
# 0.3 GENOME V5 folders
#
dirs5 <- data.frame(folder = paste(archive_dir, "/genome5/",
                                   list_dirs(paste(archive_dir, "/genome5", sep = "")), sep = ""))
dirs5$accession <- str_match(dirs5$folder, ".*/(.*)")[,2]
###################################################################################################
# 0.4 Get connected to local DB copy and get some SRA info
#
sqlfile <-
  "/Users/psenin/RNARefinery/SRAmetadb.sqlite"
sra_con <- dbConnect(SQLite(),sqlfile)
#
# define the SRA query function
#
sra_info <- function(run_accession) {
  rs <-
    dbGetQuery(sra_con,paste(
      "SELECT sra.* FROM sra WHERE run_accession=",shQuote(run_accession),sep = ""
    ))
  rs
}
dirs4$reads <-
  daply(dirs4, .(accession), function(x) {
    sra_info(x$accession)$spots
  })
dirs5$reads <-
  daply(dirs4, .(accession), function(x) {
    sra_info(x$accession)$spots
  })
###################################################################################################
# 0.5 Get the tissue data from the local file
#
tissues <- read.csv("/Users/psenin/GitHub/RNArefinery/RScripts/data/tissue.csv")
dirs4 <- merge(dirs4, tissues, by = c("accession"))
dirs5 <- merge(dirs5, tissues, by = c("accession"))
###################################################################################################
# populate other data, for this define a FIND FILE function
#
find_file <-
  function(path = ".", filename = "", recursive = FALSE) {
    print(paste("searching ", filename, " in ", path))
    all <- list.files(path, filename, full.names = TRUE,
                      recursive = recursive, ignore.case = TRUE
    )
    all <- all[!file.info(all)$isdir]
    if (length(all) > 0) {
      return(all)
    }else {
      return(NA)
    }
  }
###################################################################################################
# 0.6 Get the raw (assembled, FPKM 3) contigs LLIST
#
dirs4$contigs_llist <-
  daply(dirs4, .(accession), function(x) {
    find_file(paste(x$folder, "/", sep = ""), "raw_contigs.llist", recursive = T)
  })
dirs5$contigs_llist <-
  daply(dirs5, .(accession), function(x) {
    find_file(paste(x$folder, "/", sep = ""), "raw_contigs.llist", recursive = T)
  })
###################################################################################################
# 0.7 Get the contig stats -- it is the same procedure for both
#
assemblys_contigs_llist <- dlply(dirs4, .(accession), function(x){
  print(paste(x$contigs_llist))
  list = fread(input = as.character(x$contigs_llist), header = F)
  setnames(list, c("contig_name", "length"))
  list
})
df = ldply(assemblys_contigs_llist,function(x){length(x$contig_name)})
setnames(df,c("accession","contigs"))
dirs4 <- merge(dirs4, df)
dirs5 <- merge(dirs5, df)
###################################################################################################
# 1.0 PLOT THE ASSEMBLYS UNDER THE TEST
#
plot_assembly_by_tissue <- ggplot(data = dirs4, aes(x = reads, y = contigs,
    color = paste(accession, tissue, sep = ","))) + ggtitle("Assembly stats: reads vs contigs") +
    theme_bw() + geom_jitter(lwd = 2, alpha = 0.7) +
  guides(color = guide_legend(nrow = 4, byrow = TRUE)) +
  theme(legend.position = "bottom", legend.text = element_text(size = 8),
        legend.title = element_blank(), legend.key.size = unit(0.3, "lines"))
plot_assembly_by_tissue
#
Cairo(width = 300, height = 300, file = "rnarefinery_plots_assemblys_stats.pdf", type = "pdf",
      pointsize = 20, bg = "transparent", canvas = "white", units = "px", dpi = 60 )
print(plot_assembly_by_tissue)
dev.off()
###################################################################################################
# 2.0 Get the TSV filenames for CDS/CDNA alignment
#     since with v.4.0 genome we did use CDS/CDNA
#     and with v 5.0 we used only CDNA, we try to bring these at the same page
#
dirs4$cds4_hits <-
  daply(dirs4, .(accession), function(x) {
    find_file(paste(x$folder, "/", sep = ""), "blat_raw_cds.best.tsv", recursive = T)
  })
#
dirs4$cdna4_hits <-
  daply(dirs4, .(accession), function(x) {
    find_file(paste(x$folder, "/", sep = ""), "blat_cds_cdna.best.tsv", recursive = T)
  })
#
dirs5$cdna5_hits <-
  daply(dirs5, .(accession), function(x) {
    find_file(paste(x$folder, "/", sep = ""), "blat_cdna5.best.tsv", recursive = T)
  })

###################################################################################################
# 2.1 Define functions to read the LLIST (sequence stats file) and BEST.TSV (BLAT output file)
#
read_llist <- function(fname) {
  list = fread(input = as.character(fname), header = F)
  setnames(list, c("contig_name", "length"))
  list
}
#
read_tsv <- function(fname){
  tsv = fread(input = as.character(fname), header = T)
  setnames(tsv, gsub("%","",names(tsv)))
  tsv
}
###################################################################################################
# 2.2 Plot the CDS & CDNA v.4 alignment
#
plots_tqcover_cdna4 = list()
for (i in 1:length(dirs4$accession)) {
  #i=2
  acc = dirs4$accession[i]
  #
  incontigs_llist = read_llist(dirs4[dirs4$accession == acc,]$contigs_llist)
  tsv_cds = read_tsv(dirs4[dirs4$accession == acc,]$cds4_hits)
  tsv_cdna = read_tsv(dirs4[dirs4$accession == acc,]$cdna4_hits)
  # at this point we have CDS and CDNA hits as data frames, as well as the raw contigs data
  # remember, we excluded all that matches more than 75%:
  # exclude_list <- tsv[(tsv$"%identity" >= 75.0) &
  #        ((tsv$"%tCoverage" >= 75.0) | (tsv$"%qCoverage" >= 75.0)),]$tName
  # so.... NO HITS set
  no_hits_cds = setdiff(incontigs_llist$contig_name, tsv_cds$qName)
  no_hits_cdna = setdiff(no_hits_cds, tsv_cdna$qName)
  #
  # for those that hit, we want to see the best hit on the plot, i.e. we need to chose the best
  # hit from CDS and CDNA
  #
  tmp <- base::merge(tsv_cds, tsv_cdna, by = c("qName"), all = TRUE)
  tmp[, "max_tCov"] <- apply(select(tmp, contains("tCoverage")), 1, function(x){ max(x, na.rm = T)})
  tmp[, "max_qCov"] <- apply(select(tmp, contains("qCoverage")), 1, function(x){ max(x, na.rm = T)})
  tmp[, "max_identity"] <-
    apply(select(tmp, contains("identity.")), 1, function(x){ max(x, na.rm = T)})
  # voila -- let's plot

  title = paste("CDS & CDNA v.4: ", dirs4$accession[i], ", ",  dirs4$tissue[i],"\n",
            length(no_hits_cdna), " out of ",length(incontigs_llist$contig_name),
            " contigs were not aligned (",
            percent(length(no_hits_cdna)/length(incontigs_llist$contig_name)),")",sep = "")
  #
  plots_tqcover_cdna4[[i]] = ggplot(data = tmp, aes(x = max_tCov, y = max_qCov,
        colour = max_identity)) + geom_jitter(alpha = 0.5) + geom_density2d() +
    theme(legend.position = "bottom") +
    ggtitle(title) + scale_x_continuous(limits = c(0,100)) + scale_y_continuous(limits = c(0,100)) +
    geom_hline(aes(yintercept = 75), color = "red") +
    geom_vline(aes(xintercept = 75), color = "red") +
    scale_colour_gradientn(name = "Identity:  ",limits = c(50,100),
      colours = c("red","yellow","green","lightblue","darkblue"),
      breaks = c(50,75,100),labels = c("low(50)","medium(75)","high(100)"),
      guide = guide_colorbar(
        title.theme = element_text(size = 14, angle = 0),title.vjust = 1,
        barheight = 0.6, barwidth = 6, label.theme = element_text(size = 10, angle = 0)
      )
    )
}
Cairo(
  width = 800, height = 1200,
  file = "cdna4_alignment.pdf", type = "pdf", pointsize = 20,
  bg = "transparent", canvas = "white", units = "px", dpi = 60
)
do.call(grid.arrange,  plots_tqcover_cdna4)
dev.off()

###################################################################################################
# 2.3 plot the CDNA5 alignment
#
plots_tqcover_cdna5 = list()
for (i in 1:length(dirs5$accession)) {
  #i=2
  acc = dirs5$accession[i]
  #
  incontigs_llist = read_llist(dirs5[dirs5$accession == acc,]$contigs_llist)
  tsv = read_tsv(dirs5[dirs5$accession == acc,]$cdna5_hits)
  #
  no_hits = setdiff(incontigs_llist$contig_name, tsv$qName)
  title = paste("Gnomon mRNA v.5: ", dirs5$accession[i], ", ",  dirs5$tissue[i],"\n",length(no_hits),
                " out of ",length(incontigs_llist$contig_name), " contigs were not aligned (",
                percent(length(no_hits)/length(incontigs_llist$contig_name)),")",sep = "")
  #
  plots_tqcover_cdna5[[i]] = ggplot(data = tsv, aes(x = tCoverage, y = qCoverage,
      colour = identity)) + geom_jitter(alpha = 0.5) + geom_density2d() +
    theme(legend.position = "bottom") + ggtitle(title) + scale_x_continuous(limits = c(0,100)) +
    scale_y_continuous(limits = c(0, 100)) +
    geom_hline(aes(yintercept = 75), color = "red") +
    geom_vline(aes(xintercept = 75), color = "red") +
    scale_colour_gradientn(
      name = "Identity:  ",limits = c(50,100),
      colours = c("red","yellow","green","lightblue","darkblue"),
      breaks = c(50,75,100),labels = c("low(50)","medium(75)","high(100)"),
      guide = guide_colorbar(
        title.theme = element_text(size = 14, angle = 0),title.vjust = 1,
        barheight = 0.6, barwidth = 6, label.theme = element_text(size = 10, angle = 0)
      )
    )
}
Cairo(
  width = 800, height = 1200,
  file = "cdna5_alignment.pdf", type = "pdf", pointsize = 20,
  bg = "transparent", canvas = "white", units = "px", dpi = 60
)
do.call(grid.arrange,  plots_tqcover_cdna5)
dev.off()

###################################################################################################
# 3.0 Get the TSV filenames for both genome versions PEPTIDE alignments
#
dirs4$pep4_cdna_llist <-
  daply(dirs4, .(accession), function(x) {
    find_file(paste(x$folder, "/", sep = ""), "contigs_after_cdna.llist", recursive = T)
  })
dirs4$pep4_hits <-
  daply(dirs4, .(accession), function(x) {
    find_file(paste(x$folder, "/", sep = ""), "blat_cdna_refseq_pep.best.tsv", recursive = T)
  })
dirs5$pep5_cdna_llist <-
  daply(dirs5, .(accession), function(x) {
    find_file(paste(x$folder, "/", sep = ""), "contigs_after_cdna5.llist", recursive = T)
  })
dirs5$pep5_hits <-
  daply(dirs5, .(accession), function(x) {
    find_file(paste(x$folder, "/", sep = ""), "blat_cdna5_pep5.best.tsv", recursive = T)
  })
###################################################################################################
# 3.1. Plot the PEP4 alignment
#
plots_tqcover_pep4 = list()
llist_count=rep(0,length(dirs4$accession))
for (i in 1:length(dirs4$accession)) {
  #i=2
  acc = dirs4$accession[i]
  #
  incontigs_llist = read_llist(dirs4[dirs4$accession == acc,]$pep4_cdna_llist)
  llist_count[i] = length(incontigs_llist$contig_name)
  tsv = read_tsv(dirs4[dirs4$accession == acc,]$pep4_hits)
  #
  no_hits = setdiff(incontigs_llist$contig_name, tsv$tName)
  title = paste("Refseq PEP v.4: ", dirs4$accession[i], ", ",  dirs4$tissue[i],"\n",length(no_hits),
                " out of ",length(incontigs_llist$contig_name), " contigs were not aligned (",
                percent(length(no_hits)/length(incontigs_llist$contig_name)),")",sep = "")
  #
  plots_tqcover_pep4[[i]] = ggplot(data = tsv, aes(x = qCoverage, y = tCoverage,
       colour = identity)) + geom_jitter(alpha = 0.5) + geom_density2d() +
    theme(legend.position = "bottom") + ggtitle(title) + scale_x_continuous(limits = c(0,100)) +
    scale_y_continuous(limits = c(0,100))  +
    geom_hline(aes(yintercept = 75), color = "red") +
    geom_vline(aes(xintercept = 75), color = "red") +
    scale_colour_gradientn(
      name = "Identity:  ",limits = c(50,100),
      colours = c("red","yellow","green","lightblue","darkblue"),
      breaks = c(50,75,100),labels = c("low(50)","medium(75)","high(100)"),
      guide = guide_colorbar(
        title.theme = element_text(size = 14, angle = 0),title.vjust = 1,
        barheight = 0.6, barwidth = 6, label.theme = element_text(size = 10, angle = 0)
      )
    )
}
Cairo(
  width = 800, height = 1200,
  file = "PEP4_alignment.pdf", type = "pdf", pointsize = 20,
  bg = "transparent", canvas = "white", units = "px", dpi = 60
)
do.call(grid.arrange,  plots_tqcover_pep4)
dev.off()
dirs4$contigs_after_cdna <- llist_count
###################################################################################################
# 3.2 plot the PEP5 alignment
#
plots_tqcover_pep5 = list()
llist_count=rep(0,length(dirs5$accession))
for (i in 1:length(dirs5$accession)) {
  #i=2
  acc = dirs5$accession[i]
  #
  incontigs_llist = read_llist(dirs5[dirs5$accession == acc,]$pep5_cdna_llist)
  llist_count[i] = length(incontigs_llist$contig_name)
  tsv = read_tsv(dirs5[dirs5$accession == acc,]$pep5_hits)
  #
  no_hits = setdiff(incontigs_llist$contig_name, tsv$tName)
  title = paste("Gnomon Prot v.5: ", dirs5$accession[i], ", ",  dirs5$tissue[i],"\n",length(no_hits),
                " out of ",length(incontigs_llist$contig_name), " contigs were not aligned (",
                percent(length(no_hits)/length(incontigs_llist$contig_name)),")",sep = "")
  #
  plots_tqcover_pep5[[i]] = ggplot(data = tsv, aes(x = qCoverage, y = tCoverage,
       colour = identity)) + geom_jitter(alpha = 0.5) + geom_density2d() +
    theme(legend.position = "bottom") + ggtitle(title) + scale_x_continuous(limits = c(0, 100)) +
    scale_y_continuous(limits = c(0, 100))  +
    geom_hline(aes(yintercept = 75), color = "red") +
    geom_vline(aes(xintercept = 75), color = "red") +
    scale_colour_gradientn(
      name = "Identity:  ",limits = c(50,100),
      colours = c("red","yellow","green","lightblue","darkblue"),
      breaks = c(50,75,100),labels = c("low(50)","medium(75)","high(100)"),
      guide = guide_colorbar(
        title.theme = element_text(size = 14, angle = 0),title.vjust = 1,
        barheight = 0.6, barwidth = 6, label.theme = element_text(size = 10, angle = 0)
      )
    )
}
Cairo(
  width = 800, height = 1200,
  file = "PEP5_alignment.pdf", type = "pdf", pointsize = 20,
  bg = "transparent", canvas = "white", units = "px", dpi = 60
)
do.call(grid.arrange,  plots_tqcover_pep5)
dev.off()
dirs5$contigs_after_cdna <- llist_count
###################################################################################################
# 4.0 Find Genome and Refseq DNA files
#
dirs4$pep4_llist <-
  daply(dirs4, .(accession), function(x) {
    find_file(paste(x$folder, "/", sep = ""), "contigs_after_refseq_pep.llist", recursive = T)
  })
dirs4$refseq_dna4_hits <-
  daply(dirs4, .(accession), function(x) {
    find_file(paste(x$folder, "/", sep = ""), "blat_pep_refseq_dna.best.tsv", recursive = T)
  })
dirs4$genome4_hits <-
  daply(dirs4, .(accession), function(x) {
    find_file(paste(x$folder, "/", sep = ""), "blat_reseq_dna_genome.best.tsv", recursive = T)
  })
dirs5$pep5_llist <-
  daply(dirs5, .(accession), function(x) {
    find_file(paste(x$folder, "/", sep = ""), "contigs_after_cdna5_pep5.llist", recursive = T)
  })
dirs5$genome5_hits <-
  daply(dirs5, .(accession), function(x) {
    find_file(paste(x$folder, "/", sep = ""), "blat_cdna5_pep5_genome5.best.tsv", recursive = T)
  })

###################################################################################################
# 4.1 Plot the Refseq DNA and genome alignment
#
llist_count=rep(0,length(dirs4$accession))
plots_tqcover_genome4 = list()
for (i in 1:length(dirs4$accession)) {
  #i=2
  acc = dirs4$accession[i]
  #
  incontigs_llist = read_llist(dirs4[dirs4$accession == acc,]$pep4_llist)
  llist_count[i] = length(incontigs_llist$contig_name)
  tsv_refseq_dna4 = read_tsv(dirs4[dirs4$accession == acc,]$refseq_dna4_hits)
  tsv_genome4 = read_tsv(dirs4[dirs4$accession == acc,]$genome4_hits)
  # at this point we have CDS and CDNA hits as data frames, as well as the raw contigs data
  # remember, we excluded all that matches more than 75%:
  # exclude_list <- tsv[(tsv$"%identity" >= 75.0) &
  #        ((tsv$"%tCoverage" >= 75.0) | (tsv$"%qCoverage" >= 75.0)),]$tName
  # so.... NO HITS set
  no_hits_dna = setdiff(incontigs_llist$contig_name, tsv_refseq_dna4$qName)
  no_hits_genome = setdiff(no_hits_dna, tsv_genome4$qName)
  #
  # for those that hit, we want to see the best hit on the plot, i.e. we need to chose the best
  # hit from CDS and CDNA
  #
  tmp <- base::merge(tsv_refseq_dna4, tsv_genome4, by = c("qName"), all = TRUE)
  tmp[, "max_tCov"] <- apply(select(tmp, contains("tCoverage")), 1, function(x){ max(x, na.rm = T)})
  tmp[, "max_qCov"] <- apply(select(tmp, contains("qCoverage")), 1, function(x){ max(x, na.rm = T)})
  tmp[, "max_identity"] <-
    apply(select(tmp, contains("identity.")), 1, function(x){ max(x, na.rm = T)})
  # voila -- let's plot

  title = paste("Refseq DNA & Genome v.4: ", dirs4$accession[i], ", ",  dirs4$tissue[i],"\n",
                length(no_hits_genome), " out of ",length(incontigs_llist$contig_name),
                " contigs were not aligned (",
                percent(length(no_hits_genome)/length(incontigs_llist$contig_name)),")",sep = "")
  #
  plots_tqcover_genome4[[i]] = ggplot(data = tmp, aes(x = max_tCov, y = max_qCov,
    colour = max_identity)) + geom_jitter(alpha = 0.5) + geom_density2d() +
    theme(legend.position = "bottom") +
    ggtitle(title) + scale_x_continuous(limits = c(0,100)) + scale_y_continuous(limits = c(0,100)) +
    geom_hline(aes(yintercept = 75), color = "red") +
    geom_vline(aes(xintercept = 75), color = "red") +
    scale_colour_gradientn(name = "Identity:  ",limits = c(50,100),
      colours = c("red","yellow","green","lightblue","darkblue"),
      breaks = c(50,75,100),labels = c("low(50)","medium(75)","high(100)"),
      guide = guide_colorbar(
      title.theme = element_text(size = 14, angle = 0),title.vjust = 1,
      barheight = 0.6, barwidth = 6, label.theme = element_text(size = 10, angle = 0)
     )
    )
}
Cairo(
  width = 800, height = 1200,
  file = "genome4_alignment.pdf", type = "pdf", pointsize = 20,
  bg = "transparent", canvas = "white", units = "px", dpi = 60
)
do.call(grid.arrange,  plots_tqcover_genome4)
dev.off()
dirs4$contigs_after_pep <- llist_count

###################################################################################################
# 4.2 plot the Genome 5 alignment
#
llist_count=rep(0,length(dirs5$accession))
plots_tqcover_genome5 = list()
for (i in 1:length(dirs5$accession)) {
  #i=2
  acc = dirs5$accession[i]
  #
  incontigs_llist = read_llist(dirs5[dirs5$accession == acc,]$pep5_llist)
  llist_count[i] = length(incontigs_llist$contig_name)
  tsv = read_tsv(dirs5[dirs5$accession == acc,]$genome5_hits)
  #
  no_hits = setdiff(incontigs_llist$contig_name, tsv$qName)
  title = paste("NCBI Genome v.5: ", dirs5$accession[i], ", ",  dirs5$tissue[i],"\n",length(no_hits),
                " out of ",length(incontigs_llist$contig_name), " contigs were not aligned (",
                percent(length(no_hits)/length(incontigs_llist$contig_name)),")",sep = "")
  #
  plots_tqcover_genome5[[i]] = ggplot(data = tsv, aes(x = tCoverage, y = qCoverage,
                                                    colour = identity)) + geom_jitter(alpha = 0.5) + geom_density2d() +
    theme(legend.position = "bottom") + ggtitle(title) + scale_x_continuous(limits = c(0,100)) +
    scale_y_continuous(limits = c(0, 100)) +
    geom_hline(aes(yintercept = 75), color = "red") +
    geom_vline(aes(xintercept = 75), color = "red") +
    scale_colour_gradientn(
      name = "Identity:  ",limits = c(50,100),
      colours = c("red","yellow","green","lightblue","darkblue"),
      breaks = c(50,75,100),labels = c("low(50)","medium(75)","high(100)"),
      guide = guide_colorbar(
        title.theme = element_text(size = 14, angle = 0),title.vjust = 1,
        barheight = 0.6, barwidth = 6, label.theme = element_text(size = 10, angle = 0)
      )
    )
}
Cairo(
  width = 800, height = 1200,
  file = "genome5_alignment.pdf", type = "pdf", pointsize = 20,
  bg = "transparent", canvas = "white", units = "px", dpi = 60
)
do.call(grid.arrange,  plots_tqcover_genome5)
dev.off()
dirs5$contigs_after_pep <- llist_count

###################################################################################################
# 5.0 Let's look onto differences....
#
#
# 5.1 Get the table
dirs4$genome4_llist <-
  daply(dirs4, .(accession), function(x) {
    find_file(paste(x$folder, "/", sep = ""), "contigs_after_genome2.llist", recursive = T)
  })
llist_count=rep(0,length(dirs4$accession))
for (i in 1:length(dirs4$accession)) {
  acc = dirs4$accession[i]
  incontigs_llist = read_llist(dirs4[dirs4$accession == acc,]$genome4_llist)
  llist_count[i] = length(incontigs_llist$contig_name)
}
dirs4$contigs_after_genome <- llist_count
#
dirs5$genome5_llist <-
  daply(dirs5, .(accession), function(x) {
    find_file(paste(x$folder, "/", sep = ""), "contigs_after_cdna5_pep5_genome5.llist", recursive = T)
  })
llist_count=rep(0,length(dirs5$accession))
for (i in 1:length(dirs5$accession)) {
  acc = dirs5$accession[i]
  incontigs_llist = read_llist(dirs5[dirs5$accession == acc,]$pep5_llist)
  llist_count[i] = length(incontigs_llist$contig_name)
}
dirs5$contigs_after_genome <- llist_count

names(dirs4)
select(dirs4, accession, reads, tissue, contigs, contigs_after_cdna, contigs_after_pep,
       contigs_after_genome)

names(dirs5)
select(dirs5, accession, reads, tissue, contigs, contigs_after_cdna, contigs_after_pep,
       contigs_after_genome)


cmp_tsv = read_tsv("v4_vs_v5.best.tsv")
hist(cmp_tsv$tCoverage)
