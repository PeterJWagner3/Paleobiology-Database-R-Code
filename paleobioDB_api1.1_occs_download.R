# paleobiodb occurrences download
# duplicate ealry_age and late_age column names need to be culled

pbdb_occs <- function(taxa, show = c("ident", "phylo", "time", "paleoloc", "geo"), interval = NA) {
  require(stringr)
  require(readr)

  taxa <- str_c(taxa, collapse = ",")
  show <- str_c(show, collapse = ",")
  interval <- str_c(interval, collapse = ",")
  if(is.na(interval)==T) {
   pbdb_req <- str_c("http://paleobiodb.org/data1.1/occs/list.txt?base_name=",
        taxa, "&show=", show, "&limit=all")
  }
  else{
    pbdb_req <- str_c("http://paleobiodb.org/data1.1/occs/list.txt?base_name=",
                      taxa, "&show=", show,"&interval=", interval, "&limit=all")
  }
  out <- read_csv(pbdb_req, col_types = list(subgenus_name = col_character(),
                                             subgenus_reso = col_character(),
                                             genus_reso = col_character(),
                                             species_reso = col_character(),
                                             reid_no = col_integer(),
                                             tectonic_setting = col_character()))

  # drop duplicated early_age and late_age columns from base and time
  drop <- which(duplicated(names(out)) == T)
  out <- out[,-drop]

  out
}