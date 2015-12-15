require(rentrez)
require(XML)
require(plyr)
require(dplyr)
require(stringr)
require(data.table)
#
hits <- fread("../PerlCode/test/test.txt", sep = "\t")
#
organism_for_id <- function(id) {
  r_search <- entrez_search(db = "nuccore", term = id)
  org <- entrez_summary(db = "nuccore", id = r_search$ids)$organism
  # NCBI bans IPs which issue more than 3 requests per second
  Sys.sleep(0.5)
  org
}
#
#
#
organisms <- adply(hits$V2, 1, function(x){
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
