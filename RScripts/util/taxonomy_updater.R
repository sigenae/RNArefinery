require(rentrez)
require(XML)
require(plyr)
require(dplyr)
require(stringr)
require(data.table)
#
# load the data
#
options(echo = TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
print(paste("arguments supplied:",args))
#
# [0.1] take care about the input file names
hits <- args[1]
database_name <- args[2]
out_hits <- args[3]
#
# [0.2] read the input
#
hits <- fread(hits, sep = "\t")
setnames(hits, c("query_name","hit_name","description","identity","bits","query_length",
                 "alignment_length","gaps","evalue"))
head(hits, 3)
#
# [0.3] define the query function
#
organism_for_id <- function(id) {
  r_search <- entrez_search(db = database_name, term = id)
  org <- entrez_summary(db = database_name, id = r_search$ids)$organism
  # NCBI bans IPs which issue more than 3 requests per second
  Sys.sleep(0.33)
  org
}
#
# [0.4] run the update
#
# Caveat:
# because of this
# "Erreur dans list_to_dataframe(res, attr(.data, "split_labels"), .id, id_as_factor) : 
# Results do not have equal lengths"
# 
# we shall use lists
#
organisms <- alply(unique(hits$hit_name), 1, function(x){
  t <- str_match(x, "(?<=\\|)(.*?)(?=\\|)")[1,2]
  if (!any(is.na(t))) {
    res <- try(organism_for_id(t), silent = TRUE)
    if (inherits(res,'try-error')) {
      return(c(x, NA))
    } else {
      return(c(x,res))
    }
  }
  c(x, NA)
}, .progress = progress_text(char = ":"))
org_table <- ldply(organisms, function(x){data.frame(hit_id = x[1], organism = x[2], 
                                                     stringsAsFactors = F)})[,2:3]
row.names(org_table) <- org_table$hit_id
#
# [0.5] merge the retrieved data with the table
#
hits$organism <- daply(hits, .(query_name), function(x){org_table[paste(x$hit_name),]$organism})
#
# [0.6] save the results
#
write.table(hits, out_hits, sep = "\t", col.names = T, row.names = F, quote = T)
