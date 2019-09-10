#install.packages("paleobioDB", dependencies=TRUE)
#install.packages("magrittr", dependencies=TRUE)
#install.packages("stringr", dependencies=TRUE)
#install.packages("BioPhysConnectoR", dependencies=TRUE)
library(paleobioDB)
library(magrittr)
library(stringr)
library(BioPhysConnectoR)
# data curation scripts

# Routine to find species with occurrences but no taxonomic data
paleodb_occurrence_curation <- function(taxon,oldest,youngest,start_date,end_date,file_format) {
# taxonomy curation scripts
# Arguments:
# 	taxon: proper taxonomic name
# 	oldest: oldest geological interval from which you want new records
# 	youngest: youngest geological interval from which you want new records
# 	start_date: the oldest entered records that you want to include
# 	end_date: the most recent records that you want to include
# 	file_format: the end tag on the output files: '.xls' for Excel, '.txt', '.csv" for comma-delimited
# get temporal information for search
strat_names <- pbdb_intervals(limit="all") #read all intervals at once
#sn <- read.table("http://paleobiodb.org/data1./intervals/list.txt?scale=1", sep=',', header=T)

# MAKE SURE WE RECOGNIZE THE OLDEST AND YOUNGEST AGES!!!
if(oldest %in% strat_names$nam) {
	start_ma <- strat_names$eag[match(oldest, strat_names$nam)]
} else 
{ print("Error! ",oldest," is not a term in our chronostratigraphic lexicon!")}

if(youngest %in% strat_names$nam) {
	end_ma <- strat_names$lag[match(youngest, strat_names$nam)]
} else 
{ print("Error! ",youngest," is not a term in our chronostratigraphic lexicon!")}

#http <- paste("http://paleobiodb.org/data1.1/occs/list.txt?base_name=",taxon,"&max_ma=",start_ma,"&min_ma=",end_ma,"&show=time,loc,crmod&limit=all&modified_since=",start_date,"&modified_before=",end_date,sep="")
#http <- paste("http://paleobiodb.org/data1.1/occs/list.txt?base_name=",taxon,"&max_ma=",start_ma,"&min_ma=",end_ma,"&show=time,loc,crmod&limit=all&modified_since=",start_date,"&modified_before=",end_date,sep="")
#            http://paleobiodb.org/data1.2/occs/list.tsv?base_name=Cypraeoidea&interval=Cretaceous,Quaternary&occs_created_after=1998-01-01&occs_created_before=2017-01-01&show=ref
http <- paste("http://paleobiodb.org/data1.2/occs/list.txt?base_name=",taxon,"&interval=",oldest,",",youngest,"&occs_created_after=",start_date,"&occs_modified_before=",end_date,"&show=ref",sep="")
raw_occurrences <- read.table(http, sep=',', header=T)
recordsspc <- subset(raw_occurrences,raw_occurrences$identified_rank=="species")
recordssubspc <- subset(raw_occurrences,raw_occurrences$identified_rank=="subspecies")
occurrences <- rbind(recordsspc,recordssubspc)
finds <- dim(occurrences)[1]
# get occurrence data: if this can be modified to make it species-only, then do that.
occurrences$identified_name <- paleodb_species_occurrence_name_cleaning(occurrences$identified_name)
occurrences <- clear_na_from_matrix(occurrences,"")

for (i in 1:finds)	{
	wrds <- length(strsplit(as.character(occurrences$difference[i])," ")[[1]])
	if (wrds>3)	{
		state <- strsplit(as.character(occurrences$difference[i])," ")[[1]]
		if (state[wrds-2]=="species" && (state[wrds-1]=="not" && state[wrds]=="entered"))
			occurrences$difference[i] <- as.character("species not entered")
		}
	}

#Separate species without known authors.  
#fooA <- subset(occurrences,occurrences$difference=="species not entered")
#fooB <- subset(occurrences,occurrences$difference=="subjective synonym of, species not entered")
#fooC <- subset(occurrences,occurrences$difference=="replaced by, species not entered")
#foo <- rbind(fooA,fooB,fooC)
foo <- subset(occurrences,occurrences$difference=="species not entered")
dummy <- order(foo$identified_name)
author_needed <- foo[dummy,]	# occurrence information ONLY for species without author information
unentered_taxa <- unique(author_needed$identified_name)
unentered_taxa_ttl <- length(unique(author_needed$identified_name))
unentered_taxa_finds <- vector(length=unentered_taxa_ttl)
unentered_taxa_main_ref <- vector(length=unentered_taxa_ttl)
for (t in 1:unentered_taxa_ttl)	{
	unentered_taxa_finds[t] <- length(subset(author_needed$identified_name,author_needed$identified_name==unentered_taxa[t]))
	if (unentered_taxa_finds[t]==1)	{
		unentered_taxa_main_ref[t] <- as.character(subset(author_needed$primary_reference,author_needed$identified_name==unentered_taxa[t]))
		}	else {
		all_refs <- subset(author_needed$primary_reference,author_needed$identified_name==unentered_taxa[t])
		refs <- unique(all_refs)	# list of references for this species
		rr <- length(refs)		# number of references for this species
		if (rr==1)	{			# if only one, then it's easy
			unentered_taxa_main_ref[t] <- as.character(refs[1])
			} else {
			mr <- mxr <- 0	# counters to find the ref with most occurrences for this species
			for (r in 1:rr)	{
				if (length(subset(all_refs,all_refs==refs[r]))>mxr)	{
					mxr <- length(subset(all_refs,all_refs==refs[r]))
					mr <- r
					}
				}
			unentered_taxa_main_ref[t] <- as.character(refs[mr])	# tally the dominant reference
			}
		}
	}
# put together the "needy species" for output.
species_needy <- cbind(unentered_taxa,unentered_taxa_finds,unentered_taxa_main_ref)
colnames(species_needy) <- c("Species","Records","Primary Reference")

# separate species with known authors  
foo1 <- subset(occurrences,occurrences$identified_rank=="species")
# read through species.  Separate senior species from replaced ones.
foo1 <- subset(foo1,foo1$difference!="species not entered")
foo2 <- subset(foo1,foo1$difference=="subjective synonym of")
foo1 <- subset(foo1,foo1$difference!="subjective synonym of")
foo2 <- rbind(foo2,subset(foo1,foo1$difference=="objective synonym of"))
foo1 <- subset(foo1,foo1$difference!="objective synonym of")
foo2 <- rbind(foo2,subset(foo1,foo1$difference=="replaced by"))
foo1 <- subset(foo1,foo1$difference!="replaced by")
dummy1 <- order(foo1$accepted_name)	# "good" species
dummy2 <- order(foo2$identified_name)	# replaced speces (synonyms, homonyms, etc.)
# put together list of relevant names: accepted ones and replaced ones
if (length(dummy2)>0)	{
	author_known <- rbind(foo1[dummy1,],foo2[dummy2,])
	entered_taxa <- rbind(cbind(unique(as.character(foo1[dummy1,]$accepted_name)),unique(foo1[dummy1,]$accepted_no)),cbind(unique(as.character(foo2[dummy2,]$identified_name)),unique(foo2[dummy2,]$identified_no)))
	entered_junior <- length(unique(foo2[dummy2,]$identified_no))	# number of senior names
	} else {
	author_known <- foo1[dummy1,]
	entered_taxa <- cbind(unique(as.character(foo1[dummy1,]$accepted_name)),unique(foo1[dummy1,]$accepted_no))
	entered_junior <- 0
	}

### USE httpT to get all of the relevant taxon names.
### Use httpO to get the latest opinion.
### THEN... read through occurrences to get numbers of finds & main occurrence.
### Look up author and latest opinion from this table using MATCH
httpT <- paste("http://paleobiodb.org/data1.2/taxa/list.txt?base_name=",taxon,"&rank=species,subspecies&variant=all&rel=all_children&show=ref,refattr",sep="")
#             http://paleobiodb.org/data1.2/taxa/opinions.txt?base_name=Cypraeidae&rank=species,subspecies&op_type=all
httpO <- paste("http://paleobiodb.org/data1.2/taxa/opinions.txt?base_name=",taxon,"&rank=species,subspecies&op_type=all&show=ref,refattr",sep="")
taxon_information <- read.table(httpT, sep=',', header=T)
taxon_opinions <- read.table(httpO, sep=',', header=T)
taxon_opinions_final <- subset(taxon_opinions,taxon_opinions$opinion_type=="class")
entered_senior <- length(unique(foo1[dummy1,]$accepted_no))	# number of senior names
entered_taxa_ttl <- entered_senior+entered_junior
entered_taxa_finds <- vector(length=entered_taxa_ttl)
entered_taxa_main_ref <- vector(length=entered_taxa_ttl)
entered_taxa_last_opinion <- vector(length=entered_taxa_ttl)
entered_taxa_opinions_ttl <- vector(length=entered_taxa_ttl)

for (t in 1:entered_taxa_ttl)	{
	# find most common reference for occurrences
	if (t<=entered_senior)	{
		entered_taxa_finds[t] <- length(subset(author_known$accepted_name,author_known$accepted_name==entered_taxa[t,1]))
		} else {
		entered_taxa_finds[t] <- length(subset(author_known$identified_name,author_known$identified_name==entered_taxa[t,1]))
		}
	if (entered_taxa_finds[t]==1)	{
		if (t<=entered_senior)	{
			entered_taxa_main_ref[t] <- as.character(subset(author_known$primary_reference,author_known$accepted_name==entered_taxa[t,1]))
			} else	{
			entered_taxa_main_ref[t] <- as.character(subset(author_known$primary_reference,author_known$identified_name==entered_taxa[t,1]))
			}	# end case of replaced name
		} else {
		if (t<=entered_senior)	{
			all_refs <- subset(author_known$primary_reference,author_known$accepted_name==entered_taxa[t,1])
			} else {
			all_refs <- subset(author_known$primary_reference,author_known$identified_name==entered_taxa[t,1])
			}	# end case of replaced name
		refs <- unique(all_refs)	# list of references for this species
		rr <- length(refs)		# number of references for this species
		if (rr==1)	{			# if only one, then it's easy
			entered_taxa_main_ref[t] <- as.character(refs[1])
			} else {
			mr <- mxr <- 0	# counters to find the ref with most occurrences for this species
			for (r in 1:rr)	{
				if (length(subset(all_refs,all_refs==refs[r]))>mxr)	{
					mxr <- length(subset(all_refs,all_refs==refs[r]))
					mr <- r
					}
				}
			entered_taxa_main_ref[t] <- as.character(refs[mr])	# tally the dominant reference
			}	# end case of multiple references
		}	# end case of multiple finds
	# get taxonomic opinions
	paleodb_taxon_no <- as.numeric(entered_taxa[t,2])
	# get original taxon name & information
	tx <- match(paleodb_taxon_no,taxon_information$taxon_no)
	# if name has changed, then original number & current number will be different
	if (taxon_information$orig_no[tx]!=taxon_information$taxon_no[tx])	tx <- match(taxon_information$orig_no[tx],taxon_information$taxon_no)
	if(t==1)	{
		entered_taxon_info <- taxon_information[tx,]
		colnames(entered_taxon_info) <- colnames(taxon_information)
		entered_taxon_opinions <- subset(taxon_opinions,taxon_opinions$orig_no==taxon_information$orig_no[tx])
		} else {
		entered_taxon_info <- rbind(entered_taxon_info,taxon_information[tx,])
		entered_taxon_opinions <- rbind(entered_taxon_opinions,subset(taxon_opinions,taxon_opinions$orig_no==taxon_information$orig_no[tx]))
		}	# concatenate new information onto growing entered_taxon_info matrix
#	paleodb_taxon_no <- as.numeric(entered_taxa[t,2])
	# get original taxon name & information
#	tx <- match(paleodb_taxon_no,taxon_information$taxon_no)
	entered_taxa_opinions_ttl[t] <- dim(subset(taxon_opinions,taxon_opinions$orig_no==taxon_information$orig_no[tx]))[1]
	txo <- match(taxon_information$orig_no[tx],taxon_opinions_final$orig_no)
	entered_taxa_last_opinion[t] <- as.character(taxon_opinions_final$primary_reference[txo])
	}
species_entered <- cbind(entered_taxon_info$orig_no,entered_taxon_info$taxon_name,entered_taxon_info$accepted_no,entered_taxon_info$accepted_name,entered_taxa_finds,entered_taxon_info$ref_author,entered_taxon_info$ref_pubyr,entered_taxa_opinions_ttl,entered_taxa_last_opinion)
colnames(species_entered) <- c("orig_no","taxon_name","accepted_no","accepted_name","records","original_author(s))","ref_pubyr","total_opinions","current_opinion")

# set up file name.  If just one time unit used, then use that name instead of oldest and youngest time unit names
if (oldest!=youngest)	{
	time_range <- paste(oldest,"-",youngest,sep="")
	} else {
		time_range <- oldest
	}
filename1 <- paste("Entered",taxon,time_range,start_date,"to",end_date,sep="_")
filename2 <- paste("Unentered",taxon,time_range,start_date,"to",end_date,sep="_")
filename3 <- paste("Entered_Opinions",taxon,time_range,start_date,"to",end_date,sep="_")
filename1 <- paste(filename1,file_format,sep="")
filename2 <- paste(filename2,file_format,sep="")
filename3 <- paste(filename3,file_format,sep="")
# tab delimited or comma delimited
if (file_format!=".csv")	{
	write.table(species_entered,filename1,sep="\t",eol="\n",col.names=TRUE,row.names=FALSE)
	write.table(species_needy,filename2,sep="\t",eol="\n",col.names=TRUE,row.names=FALSE)
	write.table(entered_taxon_opinions,filename3,sep="\t",eol="\n",col.names=TRUE,row.names=FALSE)
	}
if (file_format==".csv")	{
	write.table(species_entered,filename1,sep=",",eol="\n",col.names=TRUE,row.names=FALSE)
	write.table(species_needy,filename2,sep=",",eol="\n",col.names=TRUE,row.names=FALSE)
	write.table(entered_taxon_opinions,filename3,sep=",",eol="\n",col.names=TRUE,row.names=FALSE)
	}
}

# Routine to find gaps in occurrences for taxa
paleodb_find_gaps <- function(taxon,oldest,youngest,signif_gap,file_format)	{
	# get temporal information for search
	strat_names <- pbdb_intervals(limit="all") #read all intervals at once
	
	if(oldest %in% strat_names$nam) {
		start_ma <- strat_names$eag[match(oldest, strat_names$nam)]
	} else 
	{ print("Error! ",oldest," is not a term in our chronostratigraphic lexicon!")}
	
	if(youngest %in% strat_names$nam) {
		end_ma <- strat_names$lag[match(youngest, strat_names$nam)]
	} else 
	{ print("Error! ",youngest," is not a term in our chronostratigraphic lexicon!")}
	
	# get genus names
	httpG <- paste("http://paleobiodb.org/data1.1/taxa/list.txt?name=",taxon,"&rel=all_children&show=attr&rank=genus,subgenus&limit=9999",sep="")
	genera <- read.table(httpG, sep=',', header=T)
	httpF <- paste("http://paleobiodb.org/data1.1/occs/list.txt?base_name=",taxon,"&max_ma=",start_ma,"&min_ma=",end_ma,"&show=time,loc,crmod&limit=all",sep="")
	occurrences  <- read.table(httpF, sep=',', header=T)
	
	# remove n. sp., aff., etc.
  occurrences$taxon_name <- paleodb_species_occurrence_name_cleaning(occurrences$taxon_name)
  occurrences$matched_name <- paleodb_species_occurrence_name_cleaning(occurrences$matched_name)
  occurrences$late_age <- clear_na_from_vector(occurrences$late_age,0)
	
	# sometimes taxon_name comes up NA; replace it with matched_name
  occurrences <- clear_matrix_na_with_another_cell_value(occurrences,match("taxon_name",colnames(occurrences)),match("matched_name",colnames(occurrences)))
  
  mid_age <- (occurrences$early_age+occurrences$late_age)/2
  occurrences <- cbind(occurrences,mid_age)
  finds <- length(occurrences$matched_name)
  found_genus <- vector(length=finds)
	# separate genus name
  for (i in 1:finds)  {
		found_genus[i] <- strsplit(as.character(occurrences$matched_name[i])," ")[[1]][1]
  		}
	occurrences <- cbind(occurrences,found_genus)
	gap_ma <- vector(length=finds)
	occurrences <- cbind(occurrences,gap_ma)
	collection_no_older <- occurrences$collection_no
	matched_name_older <- occurrences$matched_name
	early_age_older <- occurrences$early_age
	late_age_older <- occurrences$late_age
	early_interval_older <- occurrences$early_interval
	occurrences <- cbind(occurrences,collection_no_older,matched_name_older,early_age_older,late_age_older,early_interval_older)
	sortdata <- c("found_genus","mid_age")
	# flip-flop age so that we can get oldest members first
	occurrences$mid_age <- -1*occurrences$mid_age
	occurrences <- occurrences[do.call("order",occurrences[sortdata]),]
	occurrences$mid_age <- -1*occurrences$mid_age
	
	gaps <- 1
	k <- 0
	for (i in 1:(finds-1))  {
		j <- i+1
		k <- k+1
		if ((occurrences$taxon_rank[i]=="genus" || occurrences$taxon_rank[i]=="subgenus") || occurrences$taxon_rank[i]=="species")	{
			if (as.character(occurrences$found_genus[i])==as.character(occurrences$found_genus[j]) && (occurrences$late_age[i]-occurrences$early_age[j])>signif_gap)  {
				gap_length <- occurrences$late_age[j]-occurrences$early_age[i]
				if (gaps==1)  {
					occurrences$gap_ma[i] <- as.numeric(occurrences$late_age[j])-as.numeric(occurrences$early_age[i])
					occurrences$collection_no_older[i] <- occurrences$collection_no[j]
					# add "sp." to genus or subgenus only ids.
					if (strsplit(as.character(occurrences$matched_name_older[i])," ")[[1]][1]==1)	{
						occurrences$matched_name_older[i] <- paste(occurrences$matched_name[j],"sp.",sep=" ")
					} else if (strsplit(as.character(occurrences$matched_name_older[i])," ")[[1]][1]==2 && occurrences$matched_rank[j]=="subgenus")	{
						occurrences$matched_name_older[i] <- paste(occurrences$matched_name[j],"sp.",sep=" ")	
					}	else {
						occurrences$matched_name_older[i] <- occurrences$matched_name[j]
					}
					occurrences$early_age_older[i] <- occurrences$early_age[j]
					occurrences$late_age_older[i] <- occurrences$late_age[j]
					occurrences$early_interval_older[i] <- occurrences$early_interval[j]   	
					genus_gaps <- occurrences[i,]
					gaps <- gaps+1
				} else {
					occurrences$gap_ma[i] <- as.numeric(occurrences$late_age[j])-as.numeric(occurrences$early_age[i])
					occurrences$gap_ma[i] <- as.numeric(occurrences$late_age[j])-as.numeric(occurrences$early_age[i])
					occurrences$collection_no_older[i] <- occurrences$collection_no[j]
					occurrences$matched_name_older[i] <- occurrences$matched_name[j]
					occurrences$early_age_older[i] <- occurrences$early_age[j]
					occurrences$late_age_older[i] <- occurrences$late_age[j]
					occurrences$early_interval_older[i] <- occurrences$early_interval[j]   	
					genus_gaps <- rbind(genus_gaps,occurrences[i,])
					gaps <- gaps+1
				}
				last_gap <- j
				#		Debugging routine
				#    flurb <- paste(k,i,gaps,genus_gaps$gap_ma[gaps-1],genus_gaps$found_genus[gaps-1],occurrences$found_genus[i],sep=" ")
				#    print(flurb)
			}
		}  
	}
	output <- as.character(genus_gaps$found_genus)
	chch <- "genus"
	output <- cbind(output,as.character(genus_gaps$matched_name_older))
	chch <- cbind(chch,"taxon_name_younger")
	output <- cbind(output,genus_gaps$collection_no_older)
	chch <- cbind(chch,"collection_no_younger")
	output <- cbind(output,as.character(genus_gaps$early_interval_older))
	chch <- cbind(chch,"interval_younger")
	output <- cbind(output,as.character(genus_gaps$early_age_older))
	chch <- cbind(chch,"early_age_younger")
	output <- cbind(output,as.character(genus_gaps$late_age_older))
	chch <- cbind(chch,"late_age_younger")
	output <- cbind(output,as.character(genus_gaps$matched_name))
	chch <- cbind(chch,"taxon_name_older")
	output <- cbind(output,genus_gaps$collection_no)
	chch <- cbind(chch,"collection_no_older")
	output <- cbind(output,as.character(genus_gaps$early_interval))
	chch <- cbind(chch,"interval_older")
	output <- cbind(output,as.character(genus_gaps$early_age))
	chch <- cbind(chch,"early_age_older")
	output <- cbind(output,as.character(genus_gaps$late_age))
	colnames(output) <- cbind(chch,"late_age_older")
	#output <- cbind(output,genus_gaps$collection_no_older,as.character(genus_gaps$early_interval_older),as.character(genus_gaps$early_age_older),as.character(genus_gaps$late_age_older),as.character(genus_gaps$matched_name),genus_gaps$collection_no,as.character(genus_gaps$early_interval),as.character(genus_gaps$early_age),as.character(genus_gaps$late_age))
	
	#colnames(output) <- cbind("genus","taxon_name_older","collection_no_older","interval_older","early_age_older","late_age_older","taxon_name_younger","collection_no_younger","interval_younger","early_age_younger","late_age_younger")
	
	if (oldest!=youngest)	{
		time_range <- paste(oldest,"-",youngest,sep="")
	} else {
		time_range <- oldest	
	}
	filename <- paste("Suspicious_Gaps",taxon,time_range,sep="_")
	filename <- paste(filename,file_format,sep="")
	if (file_format!=".csv")	
		write.table(output,filename,sep="\t",eol="\n",row.names=FALSE)
	if (file_format==".csv")	
		write.table(output1,filename1,sep=",",eol="\n",row.names=FALSE)
}

# Collect information about rock units to make stratigraphy consistent
paleodb_vett_formations <- function(taxon,oldest,youngest,today,file_format)	{
	# get temporal information for search
strat_names <- pbdb_intervals(limit="all") #read all intervals at once
	
if(oldest %in% strat_names$nam) {
	start_ma <- strat_names$eag[match(oldest, strat_names$nam)]
	} else {
	print("Error! ",oldest," is not a term in our chronostratigraphic lexicon!")
	}
if(youngest %in% strat_names$nam) {
	end_ma <- strat_names$lag[match(youngest, strat_names$nam)]
	} else {
	print("Error! ",youngest," is not a term in our chronostratigraphic lexicon!")
	}
	
request <- paste("base_name=",as.character(taxon[1]),sep="")
gr <- length(taxon)
if (gr>1)	{
	for (t in 2:gr)	{
		request <- paste(request,"&base_name=",sep="")
		request <- paste(request,as.character(taxon[t]),sep="")
		}
	}
	
#	http <- paste("http://paleobiodb.org/data1.1/colls/list.txt?base_name=",taxon,"&max_ma=",start_ma,"&min_ma=",end_ma,"&show=time,stratext,crmod&limit=all",sep="")
http <- paste("http://paleobiodb.org/data1.1/colls/list.txt?",request,"&max_ma=",start_ma,"&min_ma=",end_ma,"&show=time,stratext,crmod&limit=all",sep="")
collections  <- read.table(http, sep=',', header=T)
collections <- clear_na_from_matrix(collections, "")
lists <- dim(collections)[1]
mid_age <- (collections$early_age+collections$late_age)/2
collections <- cbind(collections,mid_age)
collections$formation <- clean_rock_unit_names(collections$formation,TRUE)
collections$member <- clean_rock_unit_names(collections$member,FALSE)
	
collections$zone <- gsub("\"", "",collections$zone)
collections$zone <- gsub("\\?","",collections$zone)
collections$zone <- gsub("  "," ",collections$zone)
collections$zone <- gsub(" Zone","",collections$zone)
collections$zone <- gsub(" zone","",collections$zone)
	
# remove quotation debris
for (i in 1:dim(collections)[2])	{
	collections[,i] <- gsub("â€œ","",collections[,i])
	collections[,i] <- gsub("â€\u009d","",collections[,i])
	}
	
rock_unit <- paste(collections$formation,collections$member,sep=" - ")
rock_unit <- gsub(", NA", "",rock_unit)
rock_unit <- gsub(" \\?" ,"",rock_unit)
collections <- cbind(collections,rock_unit)

sortdata <- c("rock_unit","mid_age","early_interval","zone")
collections <- collections[do.call("order",collections[sortdata]),]

reduced_collections <- cbind(as.character(collections$formation),as.character(collections$member),as.character(collections$early_interval),as.character(collections$late_interval),as.character(collections$zone),as.numeric(collections$early_age),as.numeric(collections$late_age),as.numeric(collections$mid_age))
colnames(reduced_collections) <- cbind("formation","member","early_interval","late_interval","zone","early_age","late_age","mid_age")
b <- 1
while (collections$formation[b]=="" && collections$member[b]=="")
	b <- b+1
	
orphans <- collections[1:b,]

reduced_collections[,7] <- clear_na_from_vector(reduced_collections[,7],0)
reduced_collections[,8] <- clear_na_from_vector(reduced_collections[,8],0)
formation_info <- reduced_collections[b,]
k <- 1
for (c in (b+1):lists)	{
	d <- c-1
	if (sum(as.numeric(reduced_collections[c,]!=reduced_collections[d,]))>0 && reduced_collections[c,1]!=" - ")	{
		formation_info <- rbind(formation_info,reduced_collections[c,])
		k <- k+1		
		}	# end case of rock units with different information
	}	# end search for unique lists
	
if (oldest!=youngest)	{
	time_range <- paste(oldest,"-",youngest,sep="")
	} else {
	time_range <- oldest
	}
	
filename1 <- paste("Rock_Unit_Information",taxon,time_range,sep="_")
filename1 <- paste(filename1,file_format,sep="")
if (file_format!=".csv")	
	write.table(formation_info,filename1,sep="\t",eol="\n",row.names=FALSE)
if (file_format==".csv")	
	write.table(formation_info,filename1,sep=",",eol="\n",row.names=FALSE)

filename2 <- paste("Formationless_Collections",taxon,time_range,sep="_")
filename2 <- paste(filename2,file_format,sep="")
if (file_format!=".csv")	
	write.table(orphans,filename2,sep="\t",eol="\n",row.names=FALSE)
if (file_format==".csv")	
	write.table(orphans,filename2,sep=",",eol="\n",row.names=FALSE)
}

# Routine to read a list of taxon names and get their occurrence data
get_paleodb_occurrence_for_taxon_list <- function(taxon_set, taxon_list_name, file_format)	{
#taxon_list_name <- "Hone_Species.txt"
#taxon_set <- "Hone"
#file_format <- ".xls"
taxon_list <- read.table(file=taxon_list_name, sep="\t", header=FALSE, stringsAsFactors=TRUE)
species <- dim(taxon_list)[1]
samp <- 0
taxon_records <- vector(length=species)
for (sp in 3:species) {
	spc <- as.character(taxon_list[sp,1])
	spc <- gsub(" ","%20",spc)
  http <- paste("http://paleobiodb.org/data1.1/occs/list.txt?base_name=",spc,"&show=time,loc&limit=all",sep="")
  taxon_finds <- read.table(http, sep=',', header=T)
	taxon_records[sp] <- dim(taxon_finds)[1]
  if (taxon_records[sp]>0)  {    
    if (samp==0)  {
      occurrences <- read.table(http, sep=',', header=T)
      } else {
      occurrences <- rbind(occurrences,taxon_finds)  
      }
    samp <- samp+taxon_records[sp]
		}
  }

localities <- unique(occurrences$collection_no)
localities <- sort(localities,FALSE)
locales <- length(localities)
for (l in 1:locales)  {
  httpC <- paste("http://www.paleobiodb.org/data1.1/colls/list.txt?id=",localities[l],"&show=loc,time,strat,stratext",sep="")
  if (l==1) {
    collection_info <- read.table(httpC, sep=',', header=T)
    } else {
      nc <- read.table(httpC, sep=',', header=T)
      collection_info <- rbind(collection_info,nc)
    }
  }

recordless_taxa <- subset(taxon_list, taxon_records==0)

filename1 <- paste("Recordless",taxon_set,sep="_")
filename1 <- paste(filename1,file_format,sep="")

filename2 <- paste(taxon_set,"Localities",sep="_")
filename2 <- paste(filename2,file_format,sep="")

filename3 <- paste(taxon_set,"Records",sep="_")
filename3 <- paste(filename3,file_format,sep="")

if (file_format!=".csv")	{
    write.table(recordless_taxa,filename1,sep="\t",eol="\n",row.names=FALSE, col.names=FALSE)
    write.table(collection_info,filename2,sep="\t",eol="\n",row.names=FALSE)
    write.table(occurrences,filename3,sep="\t",eol="\n",row.names=FALSE)
  }
if (file_format==".csv")	{
	write.table(recordless_taxa,filename1,sep=",",eol="\n",row.names=FALSE)
	write.table(collection_info,filename2,sep=",",eol="\n",row.names=FALSE)
	write.table(occurrences,filename3,sep=",",eol="\n",row.names=FALSE)
 	}
}

# Routine to find species with given name (e.g., "smithi") within a suprageneric taxon
find_species_name_within_higher_taxon <- function(clade,species,sub)	{
httpM <- paste("http://www.paleobiodb.org/data1.1/taxa/list.txt?name=",clade,"&rel=all_children&show=attr,app",sep="")
members <- read.table(httpM, sep=',', header=T)
member_species <- subset(members, members$rank=="species")
if (sub==TRUE)
	member_species <- rbind(member_species,subset(members, members$rank=="subspecies"))

sp_ct <- dim(member_species)[1]

matches <- 0
if (sp_ct>0){
  for (s in 1:sp_ct)	{
  	words <- length(strsplit(as.character(member_species$taxon_name[s])," ")[[1]])
  	species_name <- strsplit(as.character(member_species$taxon_name[s])," ")[[1]][words]
  	if (species_name==species)
  		matches <- matches+1
  	}

  if (matches>0)	{	
  	poss_species <- vector(length=matches)
  	m <- 0
  	for (s in 1:sp_ct)	{
  		words <- length(strsplit(as.character(member_species$taxon_name[s])," ")[[1]])
  		species_name <- strsplit(as.character(member_species$taxon_name[s])," ")[[1]][words]
  		if (species_name==species)	{
  			m <- m+1
  			poss_species[m] <- as.character(member_species$taxon_name[s])
  			}
  		}
  	} else {
  	  poss_species <- vector(length=1)
  	  poss_species[1] <- paste("no matches")
  	  
  	}
  } else {
    poss_species <- vector(length=1)
    poss_species[1] <- paste(clade," not in database")
  }

return(poss_species)
}

# Routine to find taxa last assigned to a higher taxon that now is defunct
find_taxa_assigned_to_invalid_higher_taxon <- function(taxon,file_format)	{
httpO <- paste("http://paleobiodb.org/data1.2/taxa/opinions.txt?base_name=",taxon,"&order=hierarchy&limit=all",sep="")
opinions <- read.table(httpO, sep=',', header=T)
#opinions[,13] <- clear_na_from_vector(opinions[,13],"")
#opinions[,13] <- as.numeric(opinions[,13])
taxa <- dim(opinions)[1]
statuses <- vector(length=max(opinions$orig_no))
parent_status <- vector(length=taxa)
parent <- vector(length=max(opinions$orig_no))
parent_no <- vector(length=max(opinions$orig_no))
authority <- vector(length=max(opinions$orig_no))
probs <- 0
for (t in 1:taxa) {
  tn <- as.numeric(opinions$child_spelling_no[t])
  statuses[tn] <- as.character(opinions$status[t])
  tn <- as.numeric(opinions$orig_no[t])
  statuses[tn] <- as.character(opinions$status[t])
  parent[tn] <- as.character(opinions$parent_name[t])
  parent_no[tn] <- as.numeric(opinions$parent_no[t])
  authority[tn] <- paste(as.character(opinions$author[t]),as.character(opinions$pubyr[t]),sep=" ")
  if (t>1 && opinions$status[t]=="belongs to") {
    ht <- as.numeric(opinions$parent_no[t])
    if (statuses[ht]=="belongs to") {
      parent_status[t] <- "No worries"
      } else if (statuses[ht]=="replaced by")  {
      pt <- as.numeric(opinions$parent_no[t])
      parent_status[t] <- paste(as.character(opinions$parent_name[t]),"has been replaced by",parent[pt],"by",authority[ht],sep=" ")
      probs <- probs+1
      } else if (statuses[ht]=="subjective synonym of")  {
      pt <- as.numeric(opinions$parent_no[t])
      parent_status[t] <- paste(as.character(opinions$parent_name[t]),"was subjectively synonymized with",parent[pt],"by",authority[ht],sep=" ")
      probs <- probs+1
      } else if (statuses[ht]=="objective synonym of")  {
      pt <- as.numeric(opinions$parent_no[t])
      parent_status[t] <- paste(as.character(opinions$parent_name[t]),"was objectively synonymized with",parent[pt],"by",authority[ht],sep=" ")
      probs <- probs+1
      } else if (statuses[ht]=="nomen dubium" || (statuses[ht]=="nomen nudum" || statuses[ht]=="nomen oblitum")) {
        parent_status[t] <- paste(as.character(opinions$parent_name[t]),"is now a",statuses[ht],"according to",authority[ht],sep=" ")
        probs <- probs+1
      } else if (statuses[ht]=="invalid subgroup of") {
      parent_status[t] <- paste(as.character(opinions$parent_name[t]),"is now an invalid subgroup of",parent[pt],"according to",authority[ht],sep=" ")
      probs <- probs+1
      } else {
      parent_status[t] <- paste(as.character(opinions$parent_name[t])," is ",statuses[ht],sep="")
      probs <- probs+1
      }
    }
  }
output <- matrix(0,probs,3)
prc <- 0
for (t in 2:taxa) {
  if (opinions$status[t]=="belongs to" && parent_status[t]!="No worries"){
    prc <- prc+1
    output[prc,1] <- as.numeric(opinions$orig_no[t])
    output[prc,2] <- as.character(opinions$taxon_name[t])
    output[prc,3] <- as.character(parent_status[t])
#    print(c(as.numeric(opinions$orig_no[t]),as.character(opinions$taxon_name[t]),as.character(parent_status[t])))
    }
  }
colnames(output) <- c("orig_no","taxon_name","issue")
problem_children <- paste(taxon,"Problem_Children",sep="_")
problem_children <- paste(problem_children,file_format,sep="")

if (file_format!=".csv")  {
  write.table(output,problem_children,sep="\t",eol="\n",row.names=FALSE, col.names=TRUE)
  }
if (file_format==".csv")	{
  write.table(output,problem_children,sep=",",eol="\n",row.names=FALSE,col.names=TRUE)
  }

}

find_taxa_assigned_to_invalid_higher_taxa <- function(taxon,file_format=".xls")	{

httpO <- paste("http://paleobiodb.org/data1.2/taxa/opinions.txt?base_name=",taxon,"&order=hierarchy&limit=all",sep="")
opinions <- read.csv(httpO)
x <- opinions

val <- as.numeric(sapply(opinions$parent_no, function(x) which(x==opinions$orig_no)))

name_valid <- data.frame(id=opinions$orig_no, taxon_name=opinions$taxon_name, issue=paste(opinions$taxon_name[val],opinions$status[val],opinions$parent_name[val],"according to",opinions$author[val],opinions$pubyr[val]))

output <- name_valid[which(opinions$status[val] != "belongs to"),]

opinion_status <- opinions$status[match(output$taxon_name, opinions$taxon_name)]

output <- output[which(opinion_status=="belongs to"),]
  
problem_children <- paste(taxon,"Problem_Children",sep="_")
problem_children <- paste(problem_children,file_format,sep="")
  
if (file_format!=".csv")  {
	write.table(output,problem_children,sep="\t",eol="\n",row.names=FALSE, col.names=TRUE)
  }
if (file_format==".csv")        {
	write.csv(output,problem_children, row.names=F)
  }
  
}

find_oldest_species <- function(taxon,earliest,latest)	{
http <- paste("http://paleobiodb.org/data1.2/occs/list.tsv?base_name=",taxon,"&taxon_reso=species&interval=",earliest,",",latest,sep="")
occurrences <- read.table(http, sep='\t', header=T)
lb_oldest <- max(occurrences$max_ma)
ub_oldest <- max(occurrences$min_ma)
records <- dim(occurrences)[1]
elders <- 0
for (f in 1:records)	{
	if (occurrences$max_ma[f]>ub_oldest)	{
		elders <- elders+1
		if (elders==1)	{
			finwes <- occurrences[f,]
			}	else {
			finwes <- rbind(finwes,occurrences[f,])
			}
		}
	}
methuselahs <- paste("Oldest_",taxon,"_species.xls",sep="")
write.table(finwes,methuselahs,sep="\t",eol="\n",row.names=FALSE, col.names=TRUE)
}

# routine to find any subgenera that have subgenera....
find_inconsistent_genus_and_subgenus_ranks <- function(taxon)	{
#taxon: a higher taxon that includes genera and subgenera
tranks <- "genus,subgenus"
http <- paste("http://paleobiodb.org/data1.2/taxa/list.txt?base_name=",taxon,sep="")
http <- paste(http,"&rank=",tranks,sep="")
http <- paste(http,"&show=attr,parent",sep="")
taxonomy <- read.table(http, sep=',', header=T)
genera <- subset(taxonomy,taxonomy$taxon_rank=="genus")
subgenera <- subset(taxonomy,taxonomy$taxon_rank=="subgenus")
subgenus_names <- str_replace_all(word(subgenera$taxon_name,2), "[[:punct:]]", "")
subgenus_genera <- word(subgenera$taxon_name,1)
combined <- cbind(subgenus_genera,subgenus_names)
genus_names <- as.character(genera$taxon_name)
uniq_subgenus_genera <- unique(subgenus_genera)
gsg <- length(uniq_subgenus_genera)
np <- prob <- 0
prob_gen <- ""
for (g in 1:gsg)	{
	n <- match(uniq_subgenus_genera[g],genus_names)
	if (is.na(n)) {
		prob <- prob+1
		mtch <- 0
		if (prob==1)	{
			prob_gen <- uniq_subgenus_genera[g]
			} else prob_gen <- c(prob_gen,uniq_subgenus_genera[g])
		}	else {
		np <- np+1
		mtch <- genera$taxon_no[n]
		}
	if (g==1)	{
		matchgen_no <- mtch
		} else matchgen_no <- c(matchgen_no,mtch)
	}
if (prob>0)	{
	output <- vector(length=1)
	outlines <- 1
	for (p in 1:prob)	{
		orig_mtch <- match(prob_gen[p],combined[,2])
		if (is.na(orig_mtch))	{
			xxx <- clear_na_from_vector(match(combined[,1],prob_gen[p]),0)
			for (i in 1:length(xxx))	xxx[i] <- i*(xxx[i])
			yyy <- subset(xxx,xxx>0)
			sg <- paste(combined[yyy,1]," (",combined[yyy,2],")",sep="")
			nn <- match(sg,taxonomy$taxon_name)
			new_parent <- taxonomy$parent_name[nn]
			line1 <- paste(combined[yyy,2]," is assigned to: ",combined[yyy,1]," (",combined[yyy,2],")",sep="")
			line2 <- paste("however, ",combined[yyy,1]," is either synonymized into or corrected as ",new_parent,".",sep="")
			print(line1)
			if (outlines==1)	{
				output[outlines] <- line1
				}	else {
				outlines <- outlines+1
				output <- c(output,line1)
				}
			outlines <- outlines+1
			dum <- paste("     ",line2,sep="")
			print(dum)
			dum <- paste("",line2,sep="\t")
			output <- c(output,dum)
			}	else {
			cur_comb <- paste(combined[orig_mtch,1]," (",combined[orig_mtch,2],")",sep="")
			sgn <- match(cur_comb,subgenera$taxon_name)
			author <- str_replace_all(as.character(subgenera$taxon_attr[sgn]), "[[:punct:]]", "")
			if (is.na(author))	{
				line1 <- paste(combined[orig_mtch,2],"(author unknown)")
				} else if (author=="")	{
				line1 <- paste(combined[orig_mtch,2],"(author unknown)")
				} else {
				line1 <- paste(combined[orig_mtch,2],author)
				}
			line1 <- paste(line1," currently is ",combined[orig_mtch,1]," (",combined[orig_mtch,2],"), yet we also have:",sep="")
			print(line1)
			if (outlines==1)	{
#				output[outlines,1] <- "st"
				output[outlines] <- line1
				}	else {
				outlines <- outlines+1
#				dum[1] <- "st"
#				dum[2] <- line1
				output <- c(output,line1)
				}
			xxx <- clear_na_from_vector(match(combined[,1],prob_gen[p]),0)
			ttl <- sum(xxx)
			for (i in 1:length(xxx))	xxx[i] <- i*(xxx[i])
			yyy <- subset(xxx,xxx>0)
			if (length(yyy)==1)	{
				dy <- yyy
				yyy <- vector(length=1)
				yyy[1] <- dy
				}
			for (j in 1:length(yyy))	{
				line2 <- paste(combined[yyy[j],1]," (",combined[yyy[j],2],")",sep="")
				dum <- paste("     ",line2,sep="")
				print(dum)
				dum <- paste("",line2,sep="\t")
				outlines <- outlines+1
				output <- c(output,dum)
				}
			}
		} # end problem cases
	outfile <- paste(taxon,"_genus_and_subgenus_conflicts.txt",file_format,sep="")
	}	else print("All genera with subgenera are last ranked as genera.")
}

find_genera_with_no_occurrences("Littorinimorpha","xls")
find_genera_with_no_occurrences <- function(higher_taxon,file_format)	{
#higher_taxon: a higher taxon that includes genera and subgenera
#tranks: lower taxon ranks to include; write 'genus,subgenus' to get genera & subgenera
#file_format:
#	csv: comma-delimited
#	txt: tab-delimited text
#	tab: tab-delimited text
#	xls: tab-delimited text that usually will be opened directly by Excel
tranks <- "genus,subgenus"
http <- paste("http://paleobiodb.org/data1.2/taxa/list.txt?base_name=",higher_taxon,sep="")
http <- paste(http,"&rank=",tranks,sep="")
http <- paste(http,"&show=attr,parent",sep="")
taxonomy <- read.table(http, sep=',', header=T)
absent <- subset(taxonomy,taxonomy$n_occs==0)
invalids <- subset(absent,absent$difference!="")
valids <- subset(absent,absent$difference=="")
subgenera <- subset(taxonomy,taxonomy$taxon_rank=="subgenus")
subgenus_name <- str_replace_all(word(subgenera$taxon_name,2), "[[:punct:]]", "")
subgenera <- cbind(subgenera,subgenus_name)
cases <- dim(valids)[1]
ex <- 0
i <- 1
# look for subgenera such as Cardium (Cardium) where Cardium has occurrences
while (i<=cases)	{
	if (as.character(valids$taxon_rank[i])=="subgenus")	{
		sgn <- str_replace_all(word(valids$taxon_name[i],2), "[[:punct:]]", "")
		gen <- word(valids$taxon_name[i],1)
		gen_no <- match(sgn,taxonomy$taxon_name)
		if (is.na(gen_no)==FALSE)	{
			doh <- paste(i-1,cases-1,ex+1,as.character(valids$taxon_name[i]),sep=" ")
			if (i==1)	{
				valids <- valids[2:cases,]
				}	else if (i<cases) {
				valids <- rbind(valids[1:(i-1),],valids[(i+1):cases,])
				} else valids <- valids[1:(i-1),]
			i <- i-1
			cases <- cases-1
			ex <- ex+1
			}	
		# end case where we have a subgenus
		} else {
		sgn_no <- match(valids$taxon_name[i],subgenera$subgenus_name)
		if (is.na(sgn_no)==FALSE)	{
			if (subgenera$n_occs[sgn_no]>0)	{
				if (i==1)	{
					valids <- valids[2:cases,]
					}	else if (i<cases) {
					valids <- rbind(valids[1:(i-1),],valids[(i+1):cases,])
					} else valids <- valids[1:(i-1),]
				i <- i-1
				cases <- cases-1
				ex <- ex+1
				}
			} 	
		}	# end case where we have a genus	
	i <- i+1
	}

if (cases>0)	{
	filename <- paste(higher_taxon,"_Missing_from_PaleoDB.",file_format,sep="")
	write.table(valids,filename,sep="\t",eol="\n",col.names=TRUE,row.names=FALSE)
	}	else print("?!?!  All of the genera have occurrences!!! High Five!!!!")
}
http://paleobiodb.org/data1.2/taxa/opinions.txt?base_name=Zinolia&show=full
httpO <- "http://paleobiodb.org/data1.2/taxa/opinions.txt?base_name=Trilemna&show=full"
i <- 2
higher_taxon <- "Littorinimorpha"
file_format <- "xls"
find_genera_with_no_occurrences(higher_taxon,file_format)
genera <- subset(taxonomy,taxonomy$taxon_rank=="genus")
subgenera <- subset(taxonomy,taxonomy$taxon_rank=="subgenus")
subgenus_names <- str_replace_all(word(subgenera$taxon_name,2), "[[:punct:]]", "")
subgenus_genera <- word(subgenera$taxon_name,1)
combined <- cbind(subgenus_genera,subgenus_names)
genus_names <- as.character(genera$taxon_name)
uniq_subgenus_genera <- unique(subgenus_genera)
gsg <- length(uniq_subgenus_genera)
np <- prob <- 0
prob_gen <- ""
for (g in 1:gsg)	{
	n <- match(uniq_subgenus_genera[g],genus_names)
	if (is.na(n)) {
		prob <- prob+1
		mtch <- 0
		if (prob==1)	{
			prob_gen <- uniq_subgenus_genera[g]
			} else prob_gen <- c(prob_gen,uniq_subgenus_genera[g])
		}	else {
		np <- np+1
		mtch <- genera$taxon_no[n]
		}
	if (g==1)	{
		matchgen_no <- mtch
		} else matchgen_no <- c(matchgen_no,mtch)
	}
if (prob>0)	{
	output <- vector(length=1)
	outlines <- 1
	for (p in 1:prob)	{
		orig_mtch <- match(prob_gen[p],combined[,2])
		if (is.na(orig_mtch))	{
			xxx <- clear_na_from_vector(match(combined[,1],prob_gen[p]),0)
			for (i in 1:length(xxx))	xxx[i] <- i*(xxx[i])
			yyy <- subset(xxx,xxx>0)
			sg <- paste(combined[yyy,1]," (",combined[yyy,2],")",sep="")
			nn <- match(sg,taxonomy$taxon_name)
			new_parent <- taxonomy$parent_name[nn]
			line1 <- paste(combined[yyy,2]," is assigned to: ",combined[yyy,1]," (",combined[yyy,2],")",sep="")
			line2 <- paste("however, ",combined[yyy,1]," is either synonymized into or corrected as ",new_parent,".",sep="")
			print(line1)
			if (outlines==1)	{
				output[outlines] <- line1
				}	else {
				outlines <- outlines+1
				output <- c(output,line1)
				}
			outlines <- outlines+1
			dum <- paste("     ",line2,sep="")
			print(dum)
			dum <- paste("",line2,sep="\t")
			output <- c(output,dum)
			}	else {
			cur_comb <- paste(combined[orig_mtch,1]," (",combined[orig_mtch,2],")",sep="")
			sgn <- match(cur_comb,subgenera$taxon_name)
			author <- str_replace_all(as.character(subgenera$taxon_attr[sgn]), "[[:punct:]]", "")
			if (is.na(author))	{
				line1 <- paste(combined[orig_mtch,2],"(author unknown)")
				} else if (author=="")	{
				line1 <- paste(combined[orig_mtch,2],"(author unknown)")
				} else {
				line1 <- paste(combined[orig_mtch,2],author)
				}
			line1 <- paste(line1," currently is ",combined[orig_mtch,1]," (",combined[orig_mtch,2],"), yet we also have:",sep="")
			print(line1)
			if (outlines==1)	{
#				output[outlines,1] <- "st"
				output[outlines] <- line1
				}	else {
				outlines <- outlines+1
#				dum[1] <- "st"
#				dum[2] <- line1
				output <- c(output,line1)
				}
			xxx <- clear_na_from_vector(match(combined[,1],prob_gen[p]),0)
			ttl <- sum(xxx)
			for (i in 1:length(xxx))	xxx[i] <- i*(xxx[i])
			yyy <- subset(xxx,xxx>0)
			if (length(yyy)==1)	{
				dy <- yyy
				yyy <- vector(length=1)
				yyy[1] <- dy
				}
			for (j in 1:length(yyy))	{
				line2 <- paste(combined[yyy[j],1]," (",combined[yyy[j],2],")",sep="")
				dum <- paste("     ",line2,sep="")
				print(dum)
				dum <- paste("",line2,sep="\t")
				outlines <- outlines+1
				output <- c(output,dum)
				}
			}
		} # end problem cases
	outfile <- paste(taxon,"_genus_and_subgenus_conflicts.",file_format,sep="")
	if (file_format=="csv")	{
		flsp <- ","
		} else flsp <- "\t"
	write(output,outfile,)
	}	else print("All genera with subgenera are last ranked as genera.")
}


else {
			httpO <- paste("http://paleobiodb.org/data1.2/taxa/opinions.txt?base_name=",as.character(valids$taxon_name[i]),"&rank=genus,subgenus&show=full",sep="")
			taxon_opinions <- read.table(httpO, sep=',', header=T)
			ops <- length(taxon_opinions$spelling_reason)
			if (taxon_opinions$spelling_reason=="misspelling")	{
				if (i==1)	{
					valids <- valids[2:cases,]
					}	else if (i<cases) {
					valids <- rbind(valids[1:(i-1),],valids[(i+1):cases,])
					} else valids <- valids[1:(i-1),]
				i <- i-1
				cases <- cases-1
				ex <- ex+1
				}
			}


else {
			httpO <- paste("http://paleobiodb.org/data1.2/taxa/opinions.txt?base_name=",sgn,"&rank=genus,subgenus&show=full",sep="")
			taxon_opinions <- read.table(httpO, sep=',', header=T)
			if (taxon_opinions$spelling_reason=="misspelling")	{
				if (i==1)	{
					valids <- valids[2:cases,]
					}	else if (i<cases) {
					valids <- rbind(valids[1:(i-1),],valids[(i+1):cases,])
					} else valids <- valids[1:(i-1),]
				i <- i-1
				cases <- cases-1
				ex <- ex+1
				}
			}