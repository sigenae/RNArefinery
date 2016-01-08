#
# read the three tables
#
data_rna <- read.table("../sandbox/centroids.blastn.rna.best.taxonized.tsv",header = T)
data_dna <- read.table("../sandbox/centroids.blastn.nr.best.taxonized.tsv",header = T)
data_prot <- read.table("../sandbox/centroids.blastx.protein.table.best.taxonized.tsv",header = T)
#
# check for NAs
#
which(is.na(data_rna$organism))
which(is.na(data_prot$organism))
which(is.na(data_dna$organism))
#
# get the common to three tables list
#
common_queries <- intersect(data_dna$query_name, data_prot$query_name)
common_queries <- intersect(common_queries, data_rna$query_name)
#
# build a common table
#
common_table <- merge(
  data_rna[data_rna$query_name %in% common_queries,],
  data_dna[data_dna$query_name %in% common_queries,],
  by = c("query_name")
)
common_table <- merge(
  common_table,
  data_prot[data_prot$query_name %in% common_queries,],
  by = c("query_name")
)
#
# see its head
#
head(common_table)
head(common_table[,c(1,2,11,20)])
head(common_table[,c(1,10,19,28)])
#
# work on the taxonomical assignments of common names
#
common_org_list = unique(c(as.character(unlist(common_table[,10])),
            as.character(unlist(common_table[,19])),
            as.character(unlist(common_table[,28]))))
#
taxonomy_for_id <- function(id) {
  r_search <- entrez_search(db = "taxonomy", term = id)
  org <- entrez_summary(db = "taxonomy", id = r_search$ids)
  # NCBI bans IPs which issue more than 3 requests per second
  Sys.sleep(0.33)
  as.data.frame(org)
}

common_org_taxonomy <- alply(common_org_list, 1, function(x){
  res <- try(taxonomy_for_id(x), silent = TRUE)
  if (inherits(res,'try-error')) {
    return(c(x, NA))
  } else {
    return(c(x,res))
  }
  c(x, NA)
}, .progress = progress_text(char = ":"))

str(common_org_taxonomy[[50]])

common_org_table <- ldply(common_org_taxonomy, function(x){
  tmp = data.frame(x, stringsAsFactors = F)
  data.frame(taxid = tmp$taxid, division = tmp$division, 
             genbankdivision = tmp$genbankdivision, genus = tmp$genus, 
             species = tmp$species, subspecies = tmp$subsp, 
             scientific_name = tmp$scientificname,
             common_name = tmp$commonname)
  })

head(common_org_table)

write.table(common_org_table[,-1], "../sandbox/taxonomy_table.tsv", sep = "\t", col.names = T, 
            row.names = F, quote = T)
