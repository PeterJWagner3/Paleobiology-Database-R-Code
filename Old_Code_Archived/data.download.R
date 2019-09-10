# short script to download the data with the paleobioDB api
# can't decide if I want to use Bapst's script or the ones in paleobioDB package
# I don't want to use either! just use api directly
library(magrittr)
library(dplyr)
library(stringr)
library(readr)
library(fields)
library(tidyr)

# taxon list --------------------------

# uncomment to load all marine phyla, but chordata and plants
# all_taxa <- read_csv("https://paleobiodb.org/data1.1/taxa/list.txt?id=69296&rel=all_taxa")
#
## filter to only animal phyla minus chordata
# phyla <- all_taxa %>% filter(rank == "phylum") %>% select(taxon_name) %>%
         # filter(!grepl("^Chordata|phyta$", taxon_name))


# taxon_list <- phyla$taxon_name
taxon_list = c("Brachiopoda")

# Load timescale ----------------------
time_scale <- read_csv("https://paleobiodb.org/data1.1/intervals/list.txt?scale=1")
stage <- time_scale %>% filter(level == 5) %>% arrange(early_age)


## Marine environments ---------------
reef_env <- c("reef, buildup or bioherm","perireef or subreef","intrashelf/intraplatform reef","platform/shelf-margin reef","slope/ramp reef","basin reef")
carbonate_env <- c("carbonate indet.","peritidal","shallow subtidal indet.","open shallow subtidal","lagoonal/restricted shallow subtidal","sand shoal","reef, buildup or bioherm","perireef or subreef","intrashelf/intraplatform reef","platform/shelf-margin reef","slope/ramp reef","basin reef","deep subtidal ramp","deep subtidal shelf","deep subtidal indet.","offshore ramp","offshore shelf","offshore indet.","slope","basinal (carbonate)","basinal (siliceous)")
siliciclastic_env <- c("marine indet.","marginal marine indet.","coastal indet.","estuary/bay","lagoonal","paralic indet.","delta plain","interdistributary bay","delta front","prodelta","deltaic indet.","foreshore","shoreface","transition zone/lower shoreface","offshore","submarine fan","basinal (siliciclastic)","deep-water indet.")
marine_env <- c(carbonate_env,siliciclastic_env)

# Load data and filter to only stage resolved occs ---------------------------
occs.raw <- pbdb_occs(taxa = taxon_list)

# filter to only marine environments -----------------------------------------
occs.marine <- occs.raw %>% filter(cx_int_no %in% stage$interval_no) %>%
  filter(environment %in% marine_env)

# filter out taxonomically uncertain occurrences
# no "sp." "cf." "?", etc
# no unresolved genera
# and select useful variables ------------------
occs.keep <- occs.marine %>% filter(species_reso == "" |
                                   is.na(species_reso)) %>%
                          filter(species_name != "sp.") %>%
                          filter(genus_reso == "" | genus_reso == "n. gen.") %>%
                          filter(genus != "") %>%
                          filter(class != "") %>%
                          select(occurrence_no,
                                 taxon_name,
                                 matched_name,
                                 genus_name,
                                 genus,
                                 species_name,
                                 family,
                                 order,
                                 class,
                                 phylum,
                                 cx_int_no,
                                 early_age,
                                 late_age,
                                 paleolng,
                                 paleolat,
                                 environment) %>%
                          mutate(clgen = str_c(class, genus, sep = '_'))


# for genera and interval, measure gcd of genus, and age of interval -------
data.int <- occs.keep %>% group_by(clgen, cx_int_no) %>%
  summarise(early = max(early_age), late = min(late_age),
    gcd = max(rdist.earth(matrix(c(paleolat, paleolng), ncol = 2), miles = F)),
    oc.num = length(late_age)) %>%
  separate(clgen, c("class", "genus"), sep = '_', remove = F)

# Origination and extinction times for each genus ------------
gen.fl <- data.int %>% group_by(clgen) %>% summarise(FA = max(early), LA = min(late))

# make holocene end = 0 ---------------------
gen.fl[which(is.na(gen.fl$LA)),]$LA <- 0

# join genus ranges with interval samples --
data.par <- full_join(data.int, gen.fl, by = "clgen")

# set holocene end = 0 ---------------------
data.par[which(is.na(data.par$late)),]$late <- 0
