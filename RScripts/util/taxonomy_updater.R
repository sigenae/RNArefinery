require(rentrez)
require(XML)
require(plyr)
require(dplyr)
require(stringr)
require(data.table)
#
# load the data
#
hits <- fread("../PerlCode/test/table.txt", sep = "\t")
setnames(hits, c("query_name","hit_name","description","identity","bits","query_length",
                 "alignment_length","gaps","evalue"))
head(hits, 3)
#
# define the query function
#
organism_for_id <- function(id) {
  r_search <- entrez_search(db = "nuccore", term = id)
  org <- entrez_summary(db = "nuccore", id = r_search$ids)$organism
  # NCBI bans IPs which issue more than 3 requests per second
  Sys.sleep(0.5)
  org
}
#
# run the update
#
organisms <- adply(hits$hit_name, 1, function(x){
  t <- str_match(x, "(?<=\\|)(.*?)(?=\\|)")[1,2]
  if (!any(is.na(t))) {
    res <- try(organism_for_id(t),silent = TRUE)
    if (inherits(res,'try-error')) {
      return(c(x, NA))
    } else {
      return(c(x,res))
    }
  }
  c(x, NA)
}, .progress = progress_text(char = ":"))
#
# merg the retrieved data with the table
#
organisms = select(organisms, c(2, 3))
setnames(organisms, c("hit_name","organism"))
hits = merge(hits, organisms, by = c("hit_name"))
#
# save the results
#
write.table(hits, "hits_with_affiliation.txt", sep = "\t", col.names = T,row.names = F, quote=T)
