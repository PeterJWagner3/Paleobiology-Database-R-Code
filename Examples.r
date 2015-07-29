install.packages("paleobioDB", dependencies=TRUE)
library(paleobioDB)

# download all occurrences entered from oldest date to youngest ate.  Returns: 1) List of species without taxonomic data; 2) list of species with taxonomic data plus latest opinion information.
taxon<-"Trilobita"
oldest<-"Cambrian"
youngest<-"Ordovician"
start_date<-"2015-01-01"
end_date<-"2015-12-31"
file_format<-".xls"
# this will return all Paleozoic gastropod entered in January-April of 2015, with separate .xls files for species with and without taxonomic data. 
paleodb_occurrence_curation(taxon,oldest,youngest,start_date,end_date,file_format)

# Find genera with gaps in their range ≥ specified amount (signif_gap).  Returns list giving latest occurrence on old end of range and earliest occurrence on the young end fo the range
taxon<-"Cypraeidae"
oldest<-"Mesozoic"
youngest<-"Cenozoic"
file_format<-".xls"
signif_gap<-20
# this will give you all gaps of 20+ million years among cypraeid genera from the Mesozoic through the Cenozoic in the PaleoDB, with the output given as .xls.
paleodb_find_gaps(taxon,oldest,youngest,signif_gap,file_format)

# Find all of the different ways in which a single rock unit (Formation + Member) have been entered.  Divides combinations by different stage names and zones.  
taxon<-c("Porifera","Anthozoa","Brachiopoda","Bryozoa","Mollusca","Arthropoda","Hemichordata","Chordata","Echinodermata")
oldest<-"Paleozoic"
youngest<-"Cenozoic"
file_format<-".xls"
# this will give you all formations+members containing sponges, corals, brachiopods, etc., with all different zone and local stage assignments.  
paleodb_vett_formations(taxon,oldest,youngest,file_format)


# read a list of species and retrieve occurrence & collection data for those species.
taxon_list_name<-"Hone_Species.txt"
taxon_set<-"Hone_taxa"
file_format<-".xls"
# After reading the file "Hone_species.txt," this will retrieve occurrence and collection information for those taxa, plus a list of those taxa lacking occurrences.  
get_paleodb_occurrence_for_taxon_list(taxon_set, taxon_list_name, file_format)

# find genus names combined with a species name within some suprageneric taxon.  Good for "Alloymid n. gen. stirtoni"
clade<-"Allomyidae"
species<-"stirtoni"
sub=TRUE
# if you have an occurrence in which someone puts a known species into an indeterminate genus within a family, subfamily, etc., then this will look up species within that family to see what matches might exist.  
find_species_name_within_higher_taxon(clade,species,sub)

# download occurrence, collection and/or taxonomic data for some taxonomic group from some interval of time.
taxon<-"Caviidae"
oldest<-"Cenozoic"
youngest<-"Cenozoic"
file_format<-".xls"
get_records<-TRUE
get_localities<-TRUE
get_taxonomy<-TRUE
# if "get_records," "get_localities" and "get_taxonomy" all are set to 'TRUE' then this basically gives you everything about Paleozoic gastropods. 
paleodb_data_download(taxon,oldest,youngest,get_records,get_localities, get_taxonomy,file_format)

# download all occurrence, collection and taxonomic data for some taxonomic group
taxon<-"Cathaysiorthidae"
file_format<-".xls"
# This will just give you everything we've got on Trilobites
get_a_group(taxon,file_format)

oldest<-"Ordovician"
taxon<-"Leptaena"
youngest<-"Silurian"
start_date<-"1998-01-01"
end_date<-"2015-06-01"
file_format<-".xls"
get_records<-TRUE
get_localities<-TRUE
get_abundances<-FALSE
get_taxonomy<-TRUE
taxon_level<-"species"
# this gets occurrences, localities, taxonomy and abundances: with options for all.  (It does not mean much wihtout occurrences, however!)
paleodb_data_download(taxon,oldest,youngest,get_records,get_localities,get_abundances,get_taxonomy,file_format)
# this gets only occurrences with abundances given as specimens or individuals
paleodb_abundance_data_download(taxon,oldest,youngest,get_localities,file_format)
paleodb_taxonomic_data_download(taxon,oldest,youngest,taxon_level,file_format)

taxon<-"Leptaena"
oldest<-"Ordovician"
youngest<-"Silurian"
start_date<-"1998-01-01"
end_date<-"2015-06-01"
country<-"China"
file_format<-".xls"
get_records<-TRUE
get_localities<-TRUE
get_abundances<-FALSE
get_taxonomy<-TRUE
# this filters occurrences for particular countries or sets of countries
paleodb_data_download_by_country(taxon,oldest,youngest,get_records,get_localities,get_abundances,get_taxonomy,country,file_format)

taxon<-"Cetacea"
file_format<-".xls"
# this finds all taxa for which the last opinion places it in an invalid (junior synonym, nomen whateverum, etc.) taxon.
find_taxa_assigned_to_invalid_higher_taxon(taxon,file_format)
