#install.packages("paleobioDB", dependencies=TRUE)
library(paleobioDB)
library(stringr)

clear_na_from_matrix <- function(data, replacement)  {
size <- dim(data)
for (i in 1:size[1])	{
	for (j in 1:size[2]) if (is.na(data[i,j]))	data[i,j]<-replacement
	}
return(data)
}

clear_na_from_vector <- function(data, replacement)	{
size <- length(data)
for (i in 1:size[1])	if (is.na(data[i]))	data[i]<-replacement
return(data)
}

download_occurrence_data <- function(taxa,onset="Cambrian",end="Holocene",save_files=TRUE,output_type=".txt") {
taxa <- paste(taxa, collapse = ",")
http <- paste("http://paleobiodb.org/data1.2/occs/list.csv?base_name=",taxa,"&interval=",onset,",",end,"&show=refattr&limit=all",sep = "")
fetch <- RCurl::getURL(http)
all_finds <- utils::read.csv(text = fetch, header = TRUE, stringsAsFactors=TRUE)
species_finds <- rbind(subset(all_finds,all_finds$identified_rank=="species"),subset(all_finds,all_finds$identified_rank=="subspecies"))
cleaned_names <- sapply(as.character(species_finds$identified_name),scourgify_occurrence_identifications)
species_finds$identified_name <- cleaned_names

if (save_files)	{
	if (output_type==".csv")	{
		sepr <- ","
		}	else	{
		sepr <- "\t"
		}
	if (onset!=end)	{
		timespan <- paste(onset,"-",end,sep="")
		}	else	timespan <- onset
	output1 <- paste(timespan,"_",taxa,"_Occurrences",output_type,sep="")
	write.table(species_finds,file=output1,sep = sepr,row.names = FALSE)
	}
return(species_finds)
}

download_collection_data <- function(taxa,onset="Cambrian",end="Holocene",standardize_members=TRUE,save_files=TRUE,output_type=".txt") {
taxa <- paste(taxa, collapse = ",")
http <- paste("http://paleobiodb.org/data1.2/colls/list.csv?base_name=",taxa,"&interval=",onset,",",end,"&show=loc,paleoloc,strat,stratext,refattr",sep="")
fetch <- RCurl::getURL(http)
collections <- utils::read.csv(text = fetch, header = TRUE, stringsAsFactors=TRUE)
ttl_coll <- dim(collections)[1]

# clean up rock unit names
#clean_groups <- scourgify_rock_unit_names(named_rock_unit=collections$stratgroup,delete_rock_type=TRUE)
named_rock_unit <- collections$stratgroup
clean_groups <- sapply(as.character(named_rock_unit),scourgify_rock_unit_names,delete_rock_type=TRUE)
collections$stratgroup <- clean_groups
named_rock_unit <- collections$formation
clean_formations <- sapply(as.character(named_rock_unit),scourgify_rock_unit_names,delete_rock_type=TRUE)
collections$formation <- clean_formations
named_rock_unit <- collections$member
clean_members <- sapply(as.character(named_rock_unit),scourgify_rock_unit_names,delete_rock_type=TRUE)
collections$member <- clean_members

# standardize member/formation ranks if possible
if (standardize_members)	{
	formations <- sort(unique(clean_formations))
	members <- sort(unique(clean_members))
	confusion <- sum(members %in% formations)
	if (confusion>0)	{
		member_or_formation <- members[(1:length(members))[members %in% formations]]
		for (c in 1:confusion)	{
			# Use latest opinion.  If latest opinion is "member," then reassign
			#	all collections to the latest formation/member combo
			# If the latest opinion is tied, then go with majority rule.  If that
			#	is tied, too, then just make the damned thing a formation....
			if (member_or_formation[c]!="")	{
				vote_formation <- (1:ttl_coll)[clean_formations %in% member_or_formation[c]]
				vote_member <- (1:ttl_coll)[clean_members %in% member_or_formation[c]]
				if (max(collections$ref_pubyr[vote_formation]) > max(collections$ref_pubyr[vote_member]))	{
					### elevate member to formation in appropriate collections
					collections$formation[vote_member] <- member_or_formation[c]
					collections$member[vote_member] <- ""
					}	else if (max(collections$ref_pubyr[vote_formation]) < max(collections$ref_pubyr[vote_member]))	{
#					for (cc in 1:length(vote_formation))	{
					## get the latest opinion, and assign the rock unit as a member to that formation
					latest_opinion <- vote_member[match(max(collections$ref_pubyr[vote_member]),collections$ref_pubyr[vote_member])]
					collections$formation[vote_formation] <- collections$formation[latest_opinion]
					collections$member[vote_formation] <- member_or_formation[c]
#						}
					} else if (length(vote_formation) < length(vote_member))	{
					## get the latest opinion, and assign the rock unit as a member to that formation
					latest_opinion <- vote_member[match(max(collections$ref_pubyr[vote_member]),collections$ref_pubyr[vote_member])]
					collections$formation[vote_formation] <- collections$formation[latest_opinion]
					collections$member[vote_formation] <- member_or_formation[c]
					} else	{
					collections$formation[vote_member] <- member_or_formation[c]
					collections$member[vote_member] <- ""
					}
				}
			}
		}
	}

if (save_files)	{
	if (output_type==".csv")	{
		sepr <- ","
		}	else	{
		sepr <- "\t"
		}
	if (onset!=end)	{
		timespan <- paste(onset,"-",end,sep="")
		}	else	timespan <- onset
	output <- paste(timespan,"_",taxa,"_Collections",output_type,sep="")
	write.table(collections,file=output,sep = sepr,row.names = FALSE)
	}
return(collections)
}

# this simply gets all of the occurrence, taxonomic and locality data for a taxon.
get_a_group <- function(taxon,file_format)	{
httpO <- paste("http://paleobiodb.org/data1.1/occs/list.txt?base_name=",taxon,"&show=time,loc,crmod&limit=all",sep="")
occurrences <- read.table(httpO, sep=',', header=T)
#finds <- dim(occurrences)[1]
# get occurrence data: if this can be modified to make it species-only, then do that.
occurrences$taxon_name <- paleodb_species_occurrence_name_cleaning(occurrences$taxon_name)

paleodb_species_occurrence_name_cleaning(occurrences$taxon_name)
occ_species <- rbind(subset(occurrences, occurrences$taxon_rank=="species"))
colnames(occ_species) <- colnames(occurrences)
  
httpC <- paste("http://paleobiodb.org/data1.1/colls/list.txt?base_name=",taxon,"&show=time,stratext,geo,lithext,paleoloc,crmod&limit=all",sep="")
collections <- read.table(httpC, sep=',', header=T)
#locales <- dim(collections)[1]
  
httpT <- paste("http://paleobiodb.org/data1.1/taxa/list.txt?name=",taxon,"&rel=all_children&show=attr&limit=99999&rank=genus,subgenus,species",sep="")
taxonomy <- read.table(httpT, sep=',', header=T)
interval_names <- unique(collections$early_interval)
stages <- length(interval_names)
strat_names <- pbdb_intervals(vocab="pbdb",limit="all") #read all intervals at once
for (s in 1:stages) {
	i <- match(interval_names[s],strat_names$interval_name)
	if (s==1) {
		strat <- strat_names[i,]
    } else {
    strat <- rbind(strat,strat_names[i,])  
    }
  }
  #mid <- (early_late[,1]+early_late[,2])/2
filename1 <- paste(taxon,"Records",sep="_")
filename1 <- paste(filename1,file_format,sep="")
  
filename2 <- paste(taxon,"Localities",sep="_")
filename2 <- paste(filename2,file_format,sep="")
  
filename3 <- paste(taxon,"Chronostratigraphic_Intervals",sep="_")
filename3 <- paste(filename3,file_format,sep="")
  
sortdata <- c("early_age","late_age")
strat <- strat[do.call("order",strat[sortdata]),]
  
if (file_format!=".csv")  {
	write.table(occ_species,filename1,sep="\t",eol="\n",row.names=FALSE, col.names=TRUE)
	write.table(collections,filename2,sep="\t",eol="\n",row.names=FALSE)
	write.table(strat,filename3,sep="\t",eol="\n",row.names=FALSE)
  }
if (file_format==".csv")  {
	write.table(occ_species,filename1,sep=",",eol="\n",row.names=FALSE, col.names=TRUE)
	write.table(collections,filename2,sep=",",eol="\n",row.names=FALSE)
	write.table(strat,filename3,sep=",",eol="\n",row.names=FALSE)
  }
}

# this gets occurrence, taxonomic and/or locality data for a taxon.  Enter 'TRUE' to get records, localities and/or taxonomy.
paleodb_data_download <- function(taxon,oldest,youngest,get_records,get_localities,get_abundances,get_taxonomy,get_ranks,file_format)	{
  
# get temporal information for search
strat_names <- pbdb_intervals(limit="all") #read all intervals at once
  
if(oldest %in% strat_names$nam) {
    start_ma <- strat_names$eag[match(oldest, strat_names$nam)]
  } else 
  { print("Error! ",oldest," is not a term in our chronostratigraphic lexicon!")}
  
if(youngest %in% strat_names$nam) {
    end_ma <- strat_names$lag[match(youngest, strat_names$nam)]
  } else { 
  print("Error! ",youngest," is not a term in our chronostratigraphic lexicon!")
  }
  
if (oldest!=youngest)	{
	time_range <- paste(oldest,"-",youngest,sep="")
  } else {
  time_range <- oldest
  }
  
if (get_records==TRUE)	{
	if (get_abundances=="TRUE") {
      httpO <- paste("http://paleobiodb.org/data1.1/occs/list.txt?base_name=",taxon,"&max_ma=",start_ma,"&min_ma=",end_ma,"&show=abund,time,loc,crmod&limit=all",sep="")
    } else {
      httpO <- paste("http://paleobiodb.org/data1.1/occs/list.txt?base_name=",taxon,"&max_ma=",start_ma,"&min_ma=",end_ma,"&show=time,loc,crmod&limit=all",sep="")  
    }
	occurrences <- read.table(httpO, sep=',', header=T)
#	dodod <- subset(occurrences,occurrences$matched_rank=="species")
#	dim(occurrences)
#	dim(dodod)
    #finds <- dim(occurrences)[1]
    # get occurrence data: if this can be modified to make it species-only, then do that.
    occurrences$taxon_name <- paleodb_species_occurrence_name_cleaning(occurrences$taxon_name)
    filename <- paste(taxon,time_range,"Records",sep="_")
    filename <- paste(filename,file_format,sep="")
    if (file_format!=".csv")	{
      write.table(occurrences,filename,sep="\t",eol="\n",row.names=FALSE)
    } else if (file_format==".csv")	{
      write.table(occurrences,filename,sep=",",eol="\n",row.names=FALSE)
    }
  }	
  
  if (get_localities==TRUE)	{
    httpC <- paste("http://paleobiodb.org/data1.1/colls/list.txt?base_name=",taxon,"&max_ma=",start_ma,"&min_ma=",end_ma,"&show=time,stratext,geo,lithext,paleoloc,crmod&limit=all",sep="")
    collections <- read.table(httpC, sep=',', header=T)
    collections$formation <- clean_rock_unit_names(collections$formation,delete_rock_type=TRUE)
    collections$member <- clean_rock_unit_names(collections$member,delete_rock_type=FALSE)
    collections$late_age <- clear_na_from_vector(collections$late_age, 0)
    collections$paleomodel <- clear_na_from_vector(collections$paleomodel, "")
    collections$paleomodel2 <- clear_na_from_vector(collections$paleomodel2, "")
    collections$paleomodel3 <- clear_na_from_vector(collections$paleomodel3, "")
    collections$paleomodel4 <- clear_na_from_vector(collections$paleomodel4, "")
    collections$paleolng2 <- clear_na_from_vector(collections$paleolng2, "")
    collections$paleolng3 <- clear_na_from_vector(collections$paleolng3, "")
    collections$paleolng4 <- clear_na_from_vector(collections$paleolng4, "")
    collections$paleolat2 <- clear_na_from_vector(collections$paleolat2, "")
    collections$paleolat3 <- clear_na_from_vector(collections$paleolat3, "")
    collections$paleolat4 <- clear_na_from_vector(collections$paleolat4, "")
    collections$geoplate2 <- clear_na_from_vector(collections$geoplate2, "")
    collections$regionalsection <- clear_na_from_vector(collections$regionalsection, "")
    collections$regionalbed <- clear_na_from_vector(collections$regionalbed, "")
    collections$regionalorder <- clear_na_from_vector(collections$regionalorder, "")
    collections$collection_subset <- clear_na_from_vector(collections$collection_subset, "")
    filename <- paste(taxon,time_range,"Collections",sep="_")
    filename <- paste(filename,file_format,sep="")
    if (file_format!=".csv")	{
      write.table(collections,filename,sep="\t",eol="\n",row.names=FALSE)
    } else if (file_format==".csv")	{
      write.table(collections,filename,sep=",",eol="\n",row.names=FALSE)
    }
  }
  
  if (get_taxonomy==TRUE)	{
    httpT <- paste("http://paleobiodb.org/data1.1/taxa/list.txt?name=",taxon,"&rel=all_children&show=attr&limit=99999&rank=genus,subgenus,species",sep="")
    taxonomy <- read.table(httpT, sep=',', header=T)
    filename <- paste(taxon,time_range,"Taxonomy",sep="_")
    filename <- paste(filename,file_format,sep="")
    if (file_format!=".csv")	{
      write.table(taxonomy,filename,sep="\t",eol="\n",row.names=FALSE)
    } else if (file_format==".csv")	{
      write.table(taxonomy,filename,sep=",",eol="\n",row.names=FALSE)
    }
  }
}

# this gets occurrence, taxonomic and/or locality data for a taxon.  Enter 'TRUE' to get records, localities and/or taxonomy.  Note that it gets only occurrences with specimen or individual occurrences
paleodb_abundance_data_download <- function(taxon,oldest,youngest,get_localities,file_format)	{
	
# get temporal information for search
strat_names <- pbdb_intervals(limit="all") #read all intervals at once
	
if(oldest %in% strat_names$nam) {
	start_ma <- strat_names$eag[match(oldest, strat_names$nam)]
	} else {
	print("Error! ",oldest," is not a term in our chronostratigraphic lexicon!")}
	if(youngest %in% strat_names$nam) {
		end_ma <- strat_names$lag[match(youngest, strat_names$nam)]
		} else { 
		print("Error! ",youngest," is not a term in our chronostratigraphic lexicon!")
		}
	
	if (oldest!=youngest)	{
		time_range <- paste(oldest,"-",youngest,sep="")
		} else {
		time_range <- oldest
		}
	
	httpO <- paste("http://paleobiodb.org/data1.1/occs/list.txt?base_name=",taxon,"&max_ma=",start_ma,"&min_ma=",end_ma,"&show=abund,time,loc,crmod&limit=all",sep="")
	occurrences <- read.table(httpO, sep=',', header=T)
	abundances <- subset(occurrences,occurrences$abund_unit=="specimens")
	abundances <- rbind(abundances,subset(occurrences,occurrences$abund_unit=="individuals"))
	
	abundances$taxon_name <- paleodb_species_occurrence_name_cleaning(abundances$taxon_name)
	filename <- paste(taxon,time_range,"Abundance_Records",sep="_")
	filename <- paste(filename,file_format,sep="")
	if (file_format!=".csv")	{
		write.table(occurrences,filename,sep="\t",eol="\n",row.names=FALSE)
	} else if (file_format==".csv")	{
		write.table(occurrences,filename,sep=",",eol="\n",row.names=FALSE)
	}
	
if (get_localities==TRUE)  {
	localities <- unique(abundances$collection_no)
	coll <- length(localities)
	httpC <- paste("http://paleobiodb.org/data1.1/colls/list.txt?id=",localities[1],"&show=loc,time",sep="")
	collections <- read.table(httpC, sep=',', header=T)
	for (c in 2:coll) {
		httpC <- paste("http://paleobiodb.org/data1.1/colls/list.txt?id=",localities[c],"&show=loc,time",sep="")
	collections <- rbind(collections,read.table(httpC, sep=',', header=T))
	}
		
	collections$formation <- clean_rock_unit_names(collections$formation,delete_rock_type=TRUE)
	collections$member <- clean_rock_unit_names(collections$member,delete_rock_type=FALSE)
	collections$late_age <- clear_na_from_vector(collections$late_age, 0)
	collections$paleomodel <- clear_na_from_vector(collections$paleomodel, "")
	collections$paleomodel2 <- clear_na_from_vector(collections$paleomodel2, "")
	collections$paleomodel3 <- clear_na_from_vector(collections$paleomodel3, "")
	collections$paleomodel4 <- clear_na_from_vector(collections$paleomodel4, "")
	collections$paleolng2 <- clear_na_from_vector(collections$paleolng2, "")
	collections$paleolng3 <- clear_na_from_vector(collections$paleolng3, "")
	collections$paleolng4 <- clear_na_from_vector(collections$paleolng4, "")
	collections$paleolat2 <- clear_na_from_vector(collections$paleolat2, "")
	collections$paleolat3 <- clear_na_from_vector(collections$paleolat3, "")
	collections$paleolat4 <- clear_na_from_vector(collections$paleolat4, "")
	collections$geoplate2 <- clear_na_from_vector(collections$geoplate2, "")
	collections$regionalsection <- clear_na_from_vector(collections$regionalsection, "")
	collections$regionalbed <- clear_na_from_vector(collections$regionalbed, "")
	collections$regionalorder <- clear_na_from_vector(collections$regionalorder, "")
	collections$collection_subset <- clear_na_from_vector(collections$collection_subset, "")
		filename <- paste(taxon,time_range,"Collections",sep="_")
		filename <- paste(filename,file_format,sep="")
		if (file_format!=".csv")	{
			write.table(collections,filename,sep="\t",eol="\n",row.names=FALSE)
		} else if (file_format==".csv")	{
			write.table(collections,filename,sep=",",eol="\n",row.names=FALSE)
		}
	}
}

# this gets occurrence, taxonomic and/or locality data for a taxon.  Enter 'TRUE' to get records, localities and/or taxonomy. This retrieves occurrences only from one country.
paleodb_data_download_by_country <- function(taxon,oldest,youngest,get_records,get_localities,get_abundances,get_taxonomy,country,file_format)	{
	
	# get temporal information for search
strat_names <- pbdb_intervals(limit="all") #read all intervals at once
	
if(oldest %in% strat_names$nam) {
	start_ma <- strat_names$eag[match(oldest, strat_names$nam)]
	} else {
	print("Error! ",oldest," is not a term in our chronostratigraphic lexicon!")
	}
	
if (youngest %in% strat_names$nam) {
	end_ma <- strat_names$lag[match(youngest, strat_names$nam)]
	} else { 
	print("Error! ",youngest," is not a term in our chronostratigraphic lexicon!")
	}
	
if (oldest!=youngest)	{
	time_range <- paste(oldest,"-",youngest,sep="")
	} else {
	time_range <- oldest
	}
	
if (get_records==TRUE)	{
	if (get_abundances=="TRUE") {
		httpO <- paste("http://paleobiodb.org/data1.1/occs/list.txt?base_name=",taxon,"&max_ma=",start_ma,"&min_ma=",end_ma,"&show=abund,time,loc,crmod&limit=all",sep="")
		} else {
		httpO <- paste("http://paleobiodb.org/data1.1/occs/list.txt?base_name=",taxon,"&max_ma=",start_ma,"&min_ma=",end_ma,"&show=time,loc,crmod&limit=all",sep="")
		}
	oc <- read.table(httpO, sep=',', header=T)
	occurrences <- subset(oc,oc$cc==country)
		#finds <- dim(occurrences)[1]
		# get occurrence data: if this can be modified to make it species-only, then do that.
	occurrences$taxon_name <- paleodb_species_occurrence_name_cleaning(occurrences$taxon_name)
	filename <- paste(taxon,time_range,"Records_from",country,sep="_")
	filename <- paste(filename,file_format,sep="")
	if (file_format!=".csv")	{
		write.table(occurrences,filename,sep="\t",eol="\n",row.names=FALSE)
		} else if (file_format==".csv")	{
		write.table(occurrences,filename,sep=",",eol="\n",row.names=FALSE)
		}
	}	
	
if (get_localities==TRUE)	{
	httpC <- paste("http://paleobiodb.org/data1.1/colls/list.txt?base_name=",taxon,"&max_ma=",start_ma,"&min_ma=",end_ma,"&show=time,stratext,geo,lithext,loc,paleoloc,crmod&limit=all",sep="")
	coll <- read.table(httpC, sep=',', header=T)
	collections <- subset(coll, coll$cc==country)
	collections$formation <- clean_rock_unit_names(collections$formation,delete_rock_type=TRUE)
	collections$member <- clean_rock_unit_names(collections$member,delete_rock_type=FALSE)
	collections$late_interval <- clear_na_from_vector(collections$late_interval,"")
	collections$paleomodel <- clear_na_from_vector(collections$paleomodel, "")
	collections$paleomodel2 <- clear_na_from_vector(collections$paleomodel2, "")
	collections$paleomodel3 <- clear_na_from_vector(collections$paleomodel3, "")
	collections$paleomodel4 <- clear_na_from_vector(collections$paleomodel4, "")
	collections$paleolng2 <- clear_na_from_vector(collections$paleolng2, "")
	collections$paleolng3 <- clear_na_from_vector(collections$paleolng3, "")
	collections$paleolng4 <- clear_na_from_vector(collections$paleolng4, "")
	collections$paleolat2 <- clear_na_from_vector(collections$paleolat2, "")
	collections$paleolat3 <- clear_na_from_vector(collections$paleolat3, "")
	collections$paleolat4 <- clear_na_from_vector(collections$paleolat4, "")
	collections$geoplate2 <- clear_na_from_vector(collections$geoplate2, "")
	collections$regionalsection <- clear_na_from_vector(collections$regionalsection, "")
	collections$regionalbed <- clear_na_from_vector(collections$regionalbed, "")
	collections$regionalorder <- clear_na_from_vector(collections$regionalorder, "")
	collections$collection_subset <- clear_na_from_vector(collections$collection_subset, "")
	filename <- paste(taxon,time_range,"Collections_from",country,sep="_")
	filename <- paste(filename,file_format,sep="")
	if (file_format!=".csv")	{
		write.table(collections,filename,sep="\t",eol="\n",row.names=FALSE)
		} else if (file_format==".csv")	{
		write.table(collections,filename,sep=",",eol="\n",row.names=FALSE)
		}
	}
	
if (get_taxonomy==TRUE)	{
	httpT <- paste("http://paleobiodb.org/data1.1/taxa/list.txt?name=",taxon,"&rel=all_children&show=attr&limit=99999&rank=genus,subgenus,species",sep="")
	taxonomy <- read.table(httpT, sep=',', header=T)
	filename <- paste(taxon,time_range,"Taxonomy",sep="_")
	filename <- paste(filename,file_format,sep="")
	if (file_format!=".csv")	{
		write.table(taxonomy,filename,sep="\t",eol="\n",row.names=FALSE)
		} else if (file_format==".csv")	{
		write.table(taxonomy,filename,sep=",",eol="\n",row.names=FALSE)
		}
	}
}

# this gets occurrence, taxonomic and/or locality data for a taxon.  Enter 'TRUE' to get records, localities and/or taxonomy. This retrieves occurrences only from one country.
paleodb_data_download_by_country <- function(taxon,oldest,youngest,get_records,get_localities,get_abundances,get_taxonomy,environment,file_format)	{

if (environment=="terrestrial") {
  environs <- c("terrestrial indet.", "fluvial indet.", "alluvial fan", "channel lag", "coarse channel fill", "fine channel fill", "channel", "wet floodplain", "dry floodplain", "&quot;floodplain&quot;", "crevasse splay", "levee", "mire/swamp", "fluvial-lacustrine indet.", "delta plain", "fluvial-deltaic indet.", "lacustrine - large", "lacustrine - small", "pond", "crater lake", "lacustrine delta plain", "lacustrine interdistributary bay", "lacustrine delta front", "lacustrine prodelta", "lacustrine deltaic indet.", "lacustrine indet.", "dune", "interdune", "loess", "eolian indet.", "cave", "fissure fill", "sinkhole", "karst indet.", "tar", "mire/swamp", "spring", "glacial")
  } else if (environment=="reef") {
  environs <- c("reef, buildup or bioherm","perireef or subreef","intrashelf/intraplatform reef","platform/shelf-margin reef","slope/ramp reef","basin reef")
  } else if (environment=="carbonate")  {
  environs <- c("carbonate indet.","peritidal","shallow subtidal indet.","open shallow subtidal","lagoonal/restricted shallow subtidal","sand shoal","reef, buildup or bioherm","perireef or subreef","intrashelf/intraplatform reef","platform/shelf-margin reef","slope/ramp reef","basin reef","deep subtidal ramp","deep subtidal shelf","deep subtidal indet.","offshore ramp","offshore shelf","offshore indet.","slope","basinal (carbonate)","basinal (siliceous)")  
  } else if (environment=="siliciclastic")  {
  environs <- c("marine indet.","marginal marine indet.","coastal indet.","estuary/bay","lagoonal","paralic indet.","delta plain","interdistributary bay","delta front","prodelta","deltaic indet.","foreshore","shoreface","transition zone/lower shoreface","offshore","submarine fan","basinal (siliciclastic)","deep-water indet.")
  } else if (environment=="marine") {
  env1 <- c("carbonate indet.","peritidal","shallow subtidal indet.","open shallow subtidal","lagoonal/restricted shallow subtidal","sand shoal","reef, buildup or bioherm","perireef or subreef","intrashelf/intraplatform reef","platform/shelf-margin reef","slope/ramp reef","basin reef","deep subtidal ramp","deep subtidal shelf","deep subtidal indet.","offshore ramp","offshore shelf","offshore indet.","slope","basinal (carbonate)","basinal (siliceous)")  
  env2 <- c("marine indet.","marginal marine indet.","coastal indet.","estuary/bay","lagoonal","paralic indet.","delta plain","interdistributary bay","delta front","prodelta","deltaic indet.","foreshore","shoreface","transition zone/lower shoreface","offshore","submarine fan","basinal (siliciclastic)","deep-water indet.")
  environs <- c(env1,env2)
  }

# get temporal information for search
strat_names <- pbdb_intervals(limit="all") #read all intervals at once

if(oldest %in% strat_names$nam) {
  start_ma <- strat_names$eag[match(oldest, strat_names$nam)]
  } else {
  print("Error! ",oldest," is not a term in our chronostratigraphic lexicon!")
  }

if (youngest %in% strat_names$nam) {
  end_ma <- strat_names$lag[match(youngest, strat_names$nam)]
  } else { 
  print("Error! ",youngest," is not a term in our chronostratigraphic lexicon!")
  }

if (oldest!=youngest)	{
  time_range <- paste(oldest,"-",youngest,sep="")
  } else {
  time_range <- oldest
  }

if (get_records==TRUE)	{
  if (get_abundances=="TRUE") {
    httpO <- paste("http://paleobiodb.org/data1.1/occs/list.txt?base_name=",taxon,"&max_ma=",start_ma,"&min_ma=",end_ma,"&show=abund,time,loc,crmod&limit=all",sep="")
    } else {
    httpO <- paste("http://paleobiodb.org/data1.1/occs/list.txt?base_name=",taxon,"&max_ma=",start_ma,"&min_ma=",end_ma,"&show=time,loc,geo,crmod&limit=all",sep="")
    }
  oc <- read.table(httpO, sep=',', header=T)
  oc$environment <- clear_na_from_vector(oc$environment, "")
  occurrences <- subset(oc,oc$environment==environs)
  occurrences$reid_no <- clear_na_from_vector(occurrences$reid_no,"")
  occurrences$superceded <- clear_na_from_vector(occurrences$superceded,"")
  #finds <- dim(occurrences)[1]
  # get occurrence data: if this can be modified to make it species-only, then do that.
  occurrences$taxon_name <- paleodb_species_occurrence_name_cleaning(occurrences$taxon_name)
  filename <- paste(taxon,time_range,"Records_from",environment,sep="_")
  filename <- paste(filename,file_format,sep="")
  if (file_format!=".csv")	{
    write.table(occurrences,filename,sep="\t",eol="\n",row.names=FALSE)
    } else if (file_format==".csv")	{
    write.table(occurrences,filename,sep=",",eol="\n",row.names=FALSE)
    }
  }


if (get_localities==TRUE)	{
  httpC <- paste("http://paleobiodb.org/data1.1/colls/list.txt?base_name=",taxon,"&max_ma=",start_ma,"&min_ma=",end_ma,"&show=time,stratext,geo,lithext,loc,paleoloc,crmod&limit=all",sep="")
  coll <- read.table(httpC, sep=',', header=T)
  collections <- subset(coll, coll$cc==country)
  collections$formation <- clean_rock_unit_names(collections$formation,delete_rock_type=TRUE)
  collections$member <- clean_rock_unit_names(collections$member,delete_rock_type=FALSE)
  collections$late_interval <- clear_na_from_vector(collections$late_interval,"")
  collections$paleomodel <- clear_na_from_vector(collections$paleomodel, "")
  collections$paleomodel2 <- clear_na_from_vector(collections$paleomodel2, "")
  collections$paleomodel3 <- clear_na_from_vector(collections$paleomodel3, "")
  collections$paleomodel4 <- clear_na_from_vector(collections$paleomodel4, "")
  collections$paleolng2 <- clear_na_from_vector(collections$paleolng2, "")
  collections$paleolng3 <- clear_na_from_vector(collections$paleolng3, "")
  collections$paleolng4 <- clear_na_from_vector(collections$paleolng4, "")
  collections$paleolat2 <- clear_na_from_vector(collections$paleolat2, "")
  collections$paleolat3 <- clear_na_from_vector(collections$paleolat3, "")
  collections$paleolat4 <- clear_na_from_vector(collections$paleolat4, "")
  collections$geoplate2 <- clear_na_from_vector(collections$geoplate2, "")
  collections$regionalsection <- clear_na_from_vector(collections$regionalsection, "")
  collections$regionalbed <- clear_na_from_vector(collections$regionalbed, "")
  collections$regionalorder <- clear_na_from_vector(collections$regionalorder, "")
  collections$collection_subset <- clear_na_from_vector(collections$collection_subset, "")
  filename <- paste(taxon,time_range,"Collections_from",country,sep="_")
  filename <- paste(filename,file_format,sep="")
  if (file_format!=".csv")	{
    write.table(collections,filename,sep="\t",eol="\n",row.names=FALSE)
    } else if (file_format==".csv")	{
    write.table(collections,filename,sep=",",eol="\n",row.names=FALSE)
    }
}

if (get_taxonomy==TRUE)	{
  httpT <- paste("http://paleobiodb.org/data1.1/taxa/list.txt?name=",taxon,"&rel=all_children&show=attr&limit=99999&rank=genus,subgenus,species",sep="")
  taxonomy <- read.table(httpT, sep=',', header=T)
  filename <- paste(taxon,time_range,"Taxonomy",sep="_")
  filename <- paste(filename,file_format,sep="")
  if (file_format!=".csv")	{
    write.table(taxonomy,filename,sep="\t",eol="\n",row.names=FALSE)
	  } else if (file_format==".csv")	{
    write.table(taxonomy,filename,sep=",",eol="\n",row.names=FALSE)
	  }
	}
#occs.marine <- occs.raw %>% filter(cx_int_no %in% stage$interval_no) %>%filter(environment %in% marine_env)
}

# this gets occurrence, taxonomic and/or locality data for a taxon.  Enter 'TRUE' to get records, localities and/or taxonomy.
get_species_finds <- function(taxon,oldest,youngest,file_format)	{
  
# get temporal information for search
strat_names <- pbdb_intervals(limit="all") #read all intervals at once
  
if(oldest %in% strat_names$nam) {
    start_ma <- strat_names$eag[match(oldest, strat_names$nam)]
  } else 
  { print("Error! ",oldest," is not a term in our chronostratigraphic lexicon!")}
  
if(youngest %in% strat_names$nam) {
    end_ma <- strat_names$lag[match(youngest, strat_names$nam)]
  } else { 
  print("Error! ",youngest," is not a term in our chronostratigraphic lexicon!")
  }
  
if (oldest!=youngest)	{
	time_range <- paste(oldest,"-",youngest,sep="")
  } else {
  time_range <- oldest
  }
  
httpO <- paste("http://paleobiodb.org/data1.1/occs/list.txt?base_name=",taxon,"&max_ma=",start_ma,"&min_ma=",end_ma,"&show=time,loc,crmod&limit=all",sep="")  
first_finds <- read.table(httpO, sep=',', header=T)
occurrences <- subset(first_finds,first_finds$matched_rank=="species")
occurrences$matched_name <- paleodb_species_occurrence_name_cleaning(occurrences$taxon_name)
species <- unique(sort(occurrences$matched_name))
nsp <- length(species)
finds <- vector(length=nsp)
for (i in 1:nsp)	{
	finds[i] <- dim(subset(occurrences,occurrences$matched_name==species[i]))[1]
	}

    filename <- paste(taxon,time_range,"Records",sep="_")
    filename <- paste(filename,file_format,sep="")
    if (file_format!=".csv")	{
      write.table(occurrences,filename,sep="\t",eol="\n",row.names=FALSE)
    } else if (file_format==".csv")	{
      write.table(occurrences,filename,sep=",",eol="\n",row.names=FALSE)
    }
  }	
  
get_taxonomic_ranks <- function()	{
taxonomic_rank <- "subspecies"
taxonomic_rank <- c(taxonomic_rank,"species")
taxonomic_rank <- c(taxonomic_rank,"subgenus")
taxonomic_rank <- c(taxonomic_rank,"genus")
taxonomic_rank <- c(taxonomic_rank,"subtribe")
taxonomic_rank <- c(taxonomic_rank,"tribe")
taxonomic_rank <- c(taxonomic_rank,"subfamily")
taxonomic_rank <- c(taxonomic_rank,"family")
taxonomic_rank <- c(taxonomic_rank,"superfamily")
taxonomic_rank <- c(taxonomic_rank,"infraorder")
taxonomic_rank <- c(taxonomic_rank,"suborder")
taxonomic_rank <- c(taxonomic_rank,"order")
taxonomic_rank <- c(taxonomic_rank,"superorder")
taxonomic_rank <- c(taxonomic_rank,"infraclass")
taxonomic_rank <- c(taxonomic_rank,"subclass")
taxonomic_rank <- c(taxonomic_rank,"class")
taxonomic_rank <- c(taxonomic_rank,"superclass")
taxonomic_rank <- c(taxonomic_rank,"subphylum")
taxonomic_rank <- c(taxonomic_rank,"phylum")
taxonomic_rank <- c(taxonomic_rank,"superphylum")
taxonomic_rank <- c(taxonomic_rank,"subkingdom")
taxonomic_rank <- c(taxonomic_rank,"kingdom")
taxonomic_rank <- c(taxonomic_rank,"superkingdom")

return(taxonomic_rank)
}

get_classification_tree <- function(taxon,high_rank,low_rank,file_format)	{
taxon_ranks <- get_taxonomic_ranks()
start <- match(high_rank,taxon_ranks)
end <- match(low_rank,taxon_ranks)
tranks <- taxon_ranks[start]
for (t in (start-1):end)	tranks <- paste(tranks,",",taxon_ranks[t],sep="")
http <- paste("http://paleobiodb.org/data1.2/taxa/list.txt?base_name=",taxon,sep="")
http <- paste(http,"&rank=",tranks,sep="")
if (oldest!="")		http <- paste(http,"&interval=",oldest,",",youngest)
http <- paste(http,"&show=attr,parent",sep="")
taxonomy <- read.table(http, sep=',', header=T)
genera <- subset(taxonomy,taxonomy$taxon_rank=="genus")
subgenera <- subset(taxonomy,taxonomy$taxon_rank=="subgenus")
#word(subgenera$taxon_name[1],2)
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
#	print(c(g,mtch,uniq_subgenus_genera[g]))
	}
if (prob>0)	{
	for (p in 1:prob)	{
		orig_mtch <- match(prob_gen[p],combined[,2])
		if (is.na(orig_mtch))	{
			xxx <- clear_na_from_vector(match(combined[,1],prob_gen[p]),0)
			ttl <- sum(xxx)
			for (i in 1:length(xxx))	xxx[i] <- i*(xxx[i])
			yyy <- subset(xxx,xxx>0)
			line1 <- paste(combined[yyy,2]," is assigned to: ",combined[yyy,1]," (",combined[yyy,2],")",sep="")
			line2 <- paste("   however, ",combined[yyy,1]," is not a currently accepted genus name or spelling.",sep="")
			print(line1)
			print(line2)
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
			xxx <- clear_na_from_vector(match(combined[,1],prob_gen[p]),0)
			ttl <- sum(xxx)
			for (i in 1:length(xxx))	xxx[i] <- i*(xxx[i])
			yyy <- subset(xxx,xxx>0)
			if (length(yyy)==1)	{
				dy <- yyy
				yyy <- vector(length=1)
				yyy[1] <- dy
				}
			for (i in 1:length(yyy))	{
				line2 <- paste("     ",combined[yyy[i],1]," (",combined[yyy[i],2],")",sep="")
				print(line2)
				}
			}
		}
	}	else print("All genera with subgenera are last ranked as genera.")
}

find_inconsistent_genus_and_subgenus_ranks <- function(taxon,file_format)	{
#taxon: a higher taxon that includes genera and subgenera
#file_format:
#	csv: comma-delimited
#	txt: tab-delimited text
#	tab: tab-delimited text
#	xls: tab-delimited text that usually will be opened directly by Excel
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
	outfile <- paste(taxon,"_genus_and_subgenus_conflicts.",file_format,sep="")
	if (file_format=="csv")	{
		flsp <- ","
		} else flsp <- "\t"
	write(output,outfile,)
	}	else print("All genera with subgenera are last ranked as genera.")
}

get_oldest_members_of_a_taxon_ideal <- function(taxon)	{
eras <- c("Cryogenian","Ediacaran","Cambrian","Ordovician","Silurian","Devonian","Carboniferous","Permian","Triassic","Jurassic","Cretaceous","Paleogene","Neogene")
finds <- 0
stg <- 1
while (finds==0 && stg<length(eras))	{
	httpO <- paste("http://paleobiodb.org/data1.2/occs/list.tsv?base_name=",taxon,"&interval=",eras[stg],",",eras[stg],sep="")
#	occurrences <- read.table(httpO, sep='\t', header=T)
#	httpO <- paste("http://paleobiodb.org/data1.2/occs/list.txt?base_name=",taxon,"&interval=",eras[stg],",",eras[stg],"&show=stratext",sep="")
	occurrences <- read.table(httpO, sep=',', header=T)
	if (occurrences[1,1]=="THIS REQUEST RETURNED NO RECORDS")	{
		stg <- stg+1
		}	else {
		finds <- dim(occurrences)[1]
		httpO <- paste("http://paleobiodb.org/data1.2/occs/list.txt?base_name=",taxon,"&interval=",eras[stg],",",eras[stg],"&show=stratext",sep="")
		occurrences <- read.table(httpO, sep=',', header=T)
#		for (i in 1:dim(occurrences)[1])	for (j in 1:dim(occurrences[2]))	if (is.na(data[i,j]))	occurrences[i,j] <- ""
		}
	}
elder <- order(occurrences$max_ma,decreasing=TRUE)
filename1 <- paste(taxon,"Oldest_Records.csv",sep="_")
write.table(occurrences[elder,],filename1,sep=",",eol="\n",row.names=FALSE, col.names=TRUE)
}

get_oldest_species_of_a_taxon <- function(taxon)	{
#eras <- c("Cryogenian","Ediacaran","Cambrian","Ordovician","Silurian","Devonian","Carboniferous","Permian","Triassic","Jurassic","Cretaceous","Paleogene","Neogene")
#finds <- 0
#stg <- 1
httpO <- paste("http://paleobiodb.org/data1.2/occs/list.txt?base_name=",taxon,"&taxon_reso=species",sep="")
occurrences <- read.table(httpO, sep=',', header=T)
maxma <- max(occurrences$max_ma)
elder <- subset(occurrences,occurrences$max_ma>=(maxma-(maxma/20)))
eldest <- order(elder$max_ma,decreasing=TRUE)
filename1 <- paste(taxon,"Oldest_Records.csv",sep="_")
write.table(elder[eldest,],filename1,sep=",",eol="\n",row.names=FALSE, col.names=TRUE)
}

get_oldest_members_of_a_taxon <- function(taxon)	{
httpO <- paste("http://paleobiodb.org/data1.2/occs/list.txt?base_name=",taxon,"&show=stratext",sep="")
occurrences <- read.table(httpO, sep=',', header=T)
maxma <- max(occurrences$max_ma)
elder <- subset(occurrences,occurrences$max_ma>=(maxma-(maxma/20)))
eldest <- order(elder$max_ma,decreasing=TRUE)
filename1 <- paste(taxon,"Oldest_Records.csv",sep="_")
write.table(elder[eldest,],filename1,sep=",",eol="\n",row.names=FALSE, col.names=TRUE)
}

get_oldest_members_for_all_subtaxa_within_taxon <- function(taxon,rank)	{
#http://paleobiodb.org/data1.2/taxa/list.txt?base_name=Gastropoda&rank=order
httpT <- paste("http://paleobiodb.org/data1.2/taxa/list.txt?base_name=",taxon,"&rank=",tolower(rank),"&show=subcounts",sep="")
taxonomy <- read.table(httpT, sep=',', header=T)
taxa_w_spc <- subset(taxonomy,taxonomy$n_species>0)
subtaxa <- unique(as.character(subset(taxa_w_spc,taxa_w_spc$n_occs>0)$accepted_name))
lapply(subtaxa,get_oldest_members_of_a_taxon)
}

accio_compendium_from_file_of_taxa <- function(taxonfilename,analysis_name,output_file=TRUE)	{
taxon_name <- read.table(file=taxonfilename,stringsAsFactors = TRUE)[,1]
taxon_list <- sapply(taxon_name,scourgify_taxon_names)
compendium <- accio_compendium_for_list_of_taxa(taxon_list,analysis_name,output_file)
return(compendium)
}

accio_compendium_for_list_of_taxa <- function(taxon_list,analysis_name,output_file=TRUE)	{
# remove any funny symbols from taxon names
taxon_name <- taxon_list
taxon_list <- sapply(taxon_name,scourgify_taxon_names)
ntaxa <- length(taxon_list)

compendium <- data.frame()
#entries <- array(0,ntaxa)
tx <- 1
while (tx <= ntaxa)	{
#for (tx in 1:ntaxa)	{
	print(paste(taxon_list[tx],", ",tx," of ",ntaxa,sep=""))
	http1 <- paste("http://paleobiodb.org/data1.2/taxa/list.csv?base_name=",taxon_list[tx],"&rank=genus,subgenus&private&show=attr,app",sep="")
	accio <- RCurl::getURL(http1)
	taxon_info <- data.frame(utils::read.csv(text = accio, header = TRUE, stringsAsFactors=TRUE))

	## If the taxon is absent, then something weird will happen: kill it with fire.
	if (ncol(taxon_info)<=2 || (taxon_info[1,1]=="THIS REQUEST RETURNED NO RECORDS" || taxon_info=="THIS.REQUEST.RETURNED.NO.RECORDS"))	{
		taxon_info <- compendium[1,]
		for(i in 1:length(taxon_info))	taxon_info[i] <-"?"
		tn <- match("taxon_name",colnames(compendium))
		taxon_info[tn] <- taxon_list[tx]
		no <- match("n_occs",colnames(compendium))
		taxon_info[no] <- 0
		}

	# if multiple taxa are returned, then winnow it down to one
	if (nrow(taxon_info)>1)	{
		tt <- match(taxon_list[tx],taxon_info$taxon_name)
		if (!is.na(tt))	{
			taxon_info <- taxon_info[tt,]
			} else	{
			taxon_info <- taxon_info[1,]
			tn <- match("taxon_name",colnames(compendium))
			taxon_info[tn] <- taxon_list[tx]
			}
		}

	## clear NAs
	if (is.na(taxon_info$difference))	taxon_info$difference <- ""
	taxon_info[(1:length(taxon_info))[is.na(taxon_info)]] <- ""
#	entries[tx] <- nrow(taxon_info)
#	if (entries[tx]>1)	print(paste(tx," has problems!!!"))
	# if the taxon is in the PaleoDB, then get the oldest & youngest occurrences
	oldest_rock <- youngest_rock <- "?"
	if (max(taxon_info$n_occs)>0)	{
		http2 <- paste("http://paleobiodb.org/data1.2/occs/list.csv?base_name=",taxon_list[tx],"&show=coll,strat,refattr&limit=all",sep = "")
		accio <- RCurl::getURL(http2)
		all_finds <- utils::read.csv(text = accio, header = TRUE, stringsAsFactors=TRUE)
		all_finds <- subset(all_finds,all_finds$identified_rank=="species")
		ttl_finds <- nrow(all_finds)
		if (ttl_finds>0)	{
			named_rock_unit <- all_finds$formation
			all_finds$formation <- sapply(named_rock_unit,scourgify_rock_unit_names,delete_rock_type=TRUE)
			named_rock_unit <- all_finds$member
			all_finds$member <- sapply(named_rock_unit,scourgify_rock_unit_names,delete_rock_type=TRUE)
			### get oldest rocks
			oldest <- (1:ttl_finds)[all_finds$max_ma==max(all_finds$max_ma)]
			oldest_seds <- c()
			for (rr in 1:length(oldest))	{
				ol <- oldest[rr]
				if (all_finds$member[ol]=="")	{
					zr <- as.character(all_finds$formation[ol])
					}	else	{
					zr <- paste(as.character(all_finds$formation[ol])," (",as.character(all_finds$member[ol]),")",sep="")
					}
				if (zr!="")	oldest_seds <- c(oldest_seds,zr)
				}
			if (length(oldest_seds)==0)	oldest_seds <- ""
			oldest_seds <- sort(unique(oldest_seds))
			oldest_rock <- c()
			if (length(oldest_seds)>0)	{
				rr <- 1
				while (rr<length(oldest_seds))	{
					oldest_rock <- paste(oldest_rock,oldest_seds[rr],"; ",sep="")
					rr <- rr+1
					}
				oldest_rock <- paste(oldest_rock,oldest_seds[rr],sep="")
				}	else	oldest_rock <- ""

			### get youngest rocks
			youngest <- (1:ttl_finds)[all_finds$min_ma==min(all_finds$min_ma)]
			youngest_seds <- c()
			for (rr in 1:length(youngest))	{
				yn <- youngest[rr]
				if (all_finds$member[yn]=="")	{
					zr <- as.character(all_finds$formation[yn])
					}	else	{
					zr <- paste(as.character(all_finds$formation[yn])," (",as.character(all_finds$member[yn]),")",sep="")
					}
				if (zr!="")	youngest_seds <- c(youngest_seds,zr)
				}
			if (length(youngest_seds)==0)	youngest_seds <- ""
			youngest_seds <- sort(unique(youngest_seds))
			youngest_rock <- c()
			if (length(youngest_seds)>0)	{
				rr <- 1
				while (rr<length(youngest_seds))	{
					youngest_rock <- paste(youngest_rock,youngest_seds[rr],"; ",sep="")
					rr <- rr+1
					}
				youngest_rock <- paste(youngest_rock,youngest_seds[rr],sep="")
				}	else	youngest_rock <- ""
			}
		} else	{
		a <- match("firstapp_max_ma",colnames(taxon_info))
		z <- match("late_interval",colnames(taxon_info))
		taxon_info[a:z] <- "?"
		}

	gee <- base::t(data.frame(c(oldest_rock,youngest_rock)))
	colnames(gee) <- c("oldest_rocks","youngest_rocks")
#	gee_science <- rbind(gee_science,c(oldest_rock,youngest_rock))
	compendium <- rbind(compendium, merge(taxon_info,gee))
	tx <- tx+1
	}

if (output_file)	{
	output_file_name <- paste(analysis_name,"_Compendium.xls",sep="")	
	write.table(compendium,output_file_name,row.names = FALSE,col.names = TRUE,sep="\t")
	}
return(compendium)
}

accio_occurrences_from_file_of_taxa <- function(taxonfilename,analysis_name,output_file=TRUE)	{
taxon_name <- read.table(file=taxonfilename,stringsAsFactors = TRUE,sep="\t")[,1]
taxon_list <- sapply(taxon_name,scourgify_taxon_names)
occompendium <- accio_occurrences_for_list_of_taxa(taxon_list,analysis_name,output_file)
return(occompendium)
}

accio_occurrences_for_list_of_taxa <- function(taxon_list,analysis_name,output_file=TRUE)	{
# remove any funny symbols from taxon names
taxon_name <- taxon_list
taxon_list <- sapply(taxon_name,scourgify_taxon_names)
ntaxa <- length(taxon_list)

occompendium <- data.frame(row.names=FALSE)	# occurrences
collpendium <- data.frame(row.names=FALSE)		# collections
ttl_finds <- data.frame(row.names=FALSE)
#entries <- array(0,ntaxa)
tx <- 1
while (tx <= ntaxa)	{
#	print(paste(taxon_list[tx],", ",tx," of ",ntaxa,sep=""))
	#paste("http://paleobiodb.org/data1.2/occs/list.csv?base_name=",taxon_list[tx],"&show=coll,strat,refattr&limit=all",sep = "")	accio <- RCurl::getURL(http1)
	taxon <- gsub(" ","%20",taxon_list[tx])
	httpT <- paste("http://paleobiodb.org/data1.2/taxa/list.csv?base_name=",taxon,"&show=attr,app",sep="")
	## If the taxon is absent, then something weird will happen: kill it with fire.
	accio <- RCurl::getURL(httpT)
	taxon_info <- data.frame(utils::read.csv(text = accio, header = TRUE, stringsAsFactors=TRUE))
	if (ncol(taxon_info)<=2 || (taxon_info[1,1]=="THIS REQUEST RETURNED NO RECORDS" || taxon_info=="THIS.REQUEST.RETURNED.NO.RECORDS"))	{
		ttl_finds <- rbind(ttl_finds,cbind(taxon_list[tx],0))
		}	else	{
		http <- paste("http://paleobiodb.org/data1.2/colls/list.csv?base_name=",taxon,"&show=full",sep="")
		accio <- RCurl::getURL(http)
		taxon_finds <- data.frame(utils::read.csv(text = accio, header = TRUE, stringsAsFactors=TRUE))
		if (nrow(taxon_finds)==1 && taxon_finds$collection_no=="THIS REQUEST RETURNED NO RECORDS")	{
			n_occs <- 0
			} else	{
			n_occs <- nrow(taxon_finds)
			if (n_occs==1)	{
				taxon_finds <- clear_na_from_vector(taxon_finds, "")
				} else	{
				taxon_finds <- clear_na_from_matrix(taxon_finds, "")
				}
			collpendium <- rbind(collpendium,taxon_finds)
			occompendium <- rbind(occompendium,cbind(taxon_finds$collection_no,rep(taxon_list[tx],n_occs)))
			}
		ttl_finds <- rbind(ttl_finds,cbind(taxon_list[tx],taxon_info$n_occs[1]))
		}
	tx <- tx+1
	}

# reduce collections to just unique ones & sort them
novel_coll <- sort(unique(collpendium$collection_no))
xx <- (1:nrow(collpendium))[match(novel_coll,collpendium$collection_no)]
collpendium <- collpendium[xx,]
# put column names on the occurrences & taxon finds
colnames(occompendium) <- c("collection_no","taxon")
colnames(ttl_finds) <- c("taxon","occurrences")

if (output_file)	{
	output_file1 <- paste(analysis_name,"_Collections.xls",sep="")
	output_file2 <- paste(analysis_name,"_Occurrences.xls",sep="")
	output_file3 <- paste(analysis_name,"_Finds_per_Taxon.xls",sep="")
	write.table(collpendium,output_file1,col.names=TRUE,row.names=FALSE,sep="\t")
	write.table(occompendium,output_file2,col.names=TRUE,row.names=FALSE,sep="\t")
	write.table(ttl_finds,output_file3,col.names=TRUE,row.names=FALSE,sep="\t")
	}
pendia <- list(occompendium,collpendium,ttl_finds)
return(pendia)
}

what_am_I_doing_with_my_life <- function(taxonfilename,analysis_name,output_file=TRUE)	{
taxon_name <- read.table(file=taxonfilename,stringsAsFactors = TRUE,sep="\t")[,1]
taxon_list <- sapply(taxon_name,scourgify_taxon_names)
what <- what_am_I_doing(taxon_list,analysis_name,output_file=TRUE)
}

what_am_I_doing <- function(taxonfilename,analysis_name,output_file=TRUE)	{
taxon_name <- taxon_list
taxon_list <- sapply(taxon_name,scourgify_taxon_names)
ntaxa <- length(taxon_list)

tactics <- data.frame(row.names=FALSE)	# occurrences
tx <- 1
httpD <- "http://paleobiodb.org/data1.2/taxa/list.csv?base_name=Lophospira&rank=genus&show=ecospace"
accio <- RCurl::getURL(httpD)
default <- data.frame(utils::read.csv(text = accio, header = TRUE, stringsAsFactors=TRUE))
dd <- match("taxon_environment",names(default))
for (d in dd:length(default))	default[d] <- "?"
default$orig_no <- default$taxon_no <- default$accepted_no <- default$parent_no <- default$reference_no<- default$n_occs <- 0
while (tx <= ntaxa)	{
#	print(paste(taxon_list[tx],", ",tx," of ",ntaxa,sep=""))
	#paste("http://paleobiodb.org/data1.2/occs/list.csv?base_name=",taxon_list[tx],"&show=coll,strat,refattr&limit=all",sep = "")	accio <- RCurl::getURL(http1)
	taxon <- gsub(" ","%20",taxon_list[tx])
	httpT <- paste("http://paleobiodb.org/data1.2/taxa/list.csv?base_name=",taxon,"&show=ecospace",sep="")
	## If the taxon is absent, then something weird will happen: kill it with fire.
	accio <- RCurl::getURL(httpT)
	taxon_info <- data.frame(utils::read.csv(text = accio, header = TRUE, stringsAsFactors=TRUE))
	if (ncol(taxon_info)<=2 || (taxon_info[1,1]=="THIS REQUEST RETURNED NO RECORDS" || taxon_info=="THIS.REQUEST.RETURNED.NO.RECORDS"))	{
		if (length(unique(strsplit(taxon_list[tx],split=" ")[[1]]))>1)	{
			higher_taxon <- unique(strsplit(taxon_list[tx],split=" ")[[1]])[1]
			httpHT <- paste("http://paleobiodb.org/data1.2/taxa/list.csv?base_name=",higher_taxon,"&rank=genus&show=ecospace",sep="")
			accio <- RCurl::getURL(httpHT)
			taxon_info <- data.frame(utils::read.csv(text = accio, header = TRUE, stringsAsFactors=TRUE))
			if (ncol(taxon_info)<=2 || (taxon_info[1,1]=="THIS REQUEST RETURNED NO RECORDS" || taxon_info=="THIS.REQUEST.RETURNED.NO.RECORDS"))	{
				taxon_info <- default
				taxon_info$taxon_name <- taxon_list[tx]
				taxon_info$accepted_name <- "Please Check!"
				} else	{
				taxon_info$taxon_name <- taxon_list[tx]
				}
			} else	{
			taxon_info <- default
			taxon_info$taxon_name <- taxon_list[tx]
			taxon_info$accepted_name <- "Please Check!"
			}
		}
	tactics <- rbind(tactics,taxon_info)
	tx <- tx+1
	}
if (output_file)	{
	output_file_name <- paste(analysis_name,"_Life_Habits.xls",sep="")
	write.table(tactics,output_file_name,col.names=TRUE,row.names=FALSE,sep="\t")
	}
return(tactics)
}

scourgify_rock_unit_names <- function(named_rock_unit,delete_rock_type=TRUE)	{
# named_rock_unit: string giving the name of a formation, member or group
# delete_rock_type: if true, the "Burgess Shale" becomes "Burgess" This is here because workers are
#	inconsistent about including rock-types in formation names
if (is.na(named_rock_unit))	named_rock_unit <- ""
named_rock_unit <- gsub("\\?","",named_rock_unit)
named_rock_unit <- gsub("\"", "",named_rock_unit)
named_rock_unit <- gsub("Ft. \\.","Fort ",named_rock_unit)
named_rock_unit <- gsub(" Fm\\.", "",named_rock_unit)
named_rock_unit <- gsub(" Fm", "",named_rock_unit)
named_rock_unit <- gsub(" Formation", "",named_rock_unit)
named_rock_unit <- gsub(" formation", "",named_rock_unit)
named_rock_unit <- gsub(" Member", "",named_rock_unit)
named_rock_unit <- gsub(" member", "",named_rock_unit)
named_rock_unit <- gsub(" Mbr ", "",named_rock_unit)
named_rock_unit <- gsub(" Mbr", "",named_rock_unit)
named_rock_unit <- gsub("Á","A",named_rock_unit)
named_rock_unit <- gsub("Ä","A",named_rock_unit)
named_rock_unit <- gsub("ä","a",named_rock_unit)
named_rock_unit <- gsub("á","a",named_rock_unit)
named_rock_unit <- gsub("ç","c",named_rock_unit)
named_rock_unit <- gsub("é","e",named_rock_unit)
named_rock_unit <- gsub("è","e",named_rock_unit)
named_rock_unit <- gsub("ñ","n",named_rock_unit)
named_rock_unit <- gsub("Ö","O",named_rock_unit)
named_rock_unit <- gsub("Ø","O",named_rock_unit)
named_rock_unit <- gsub("ø","o",named_rock_unit)
named_rock_unit <- gsub("ó","o",named_rock_unit)
named_rock_unit <- gsub("ö","o",named_rock_unit)
named_rock_unit <- gsub("õ","o",named_rock_unit)
named_rock_unit <- gsub("ü","u",named_rock_unit)
	#named_rock_unit <- gsub("’","\\'",,named_rock_unit)
named_rock_unit <- gsub("Ã„","A",named_rock_unit)
named_rock_unit <- gsub("Á","A",named_rock_unit)
named_rock_unit <- gsub("Ã¡","a",named_rock_unit)
named_rock_unit <- gsub("Ã¤","a",named_rock_unit)
named_rock_unit <- gsub("Ã¥","a",named_rock_unit)
named_rock_unit <- gsub("Ã§","c",named_rock_unit)
named_rock_unit <- gsub("Ã©","e",named_rock_unit)
named_rock_unit <- gsub("Ã¨","e",named_rock_unit)
named_rock_unit <- gsub("Ã±","n",named_rock_unit)
named_rock_unit <- gsub("Ã–","O",named_rock_unit)
named_rock_unit <- gsub("Ã¸","o",named_rock_unit)
named_rock_unit <- gsub("Ã¶","o",named_rock_unit)
named_rock_unit <- gsub("Ãµ","o",named_rock_unit)
named_rock_unit <- gsub("Ã´","o",named_rock_unit)
named_rock_unit <- gsub("Ã¼","u",named_rock_unit)
	#named_rock_unit <- gsub("â€™","\\'",,named_rock_unit)
	
if (delete_rock_type)	{
	named_rock_unit <- gsub(" Ls", "",named_rock_unit)
	named_rock_unit <- gsub(" Lst", "",named_rock_unit)
	named_rock_unit <- gsub(" Ls\\.", "",named_rock_unit)
	named_rock_unit <- gsub(" Lst\\.", "",named_rock_unit)
	named_rock_unit <- gsub(" Qzt\\.", "",named_rock_unit)
	named_rock_unit <- gsub(" Quartzite", "",named_rock_unit)
	named_rock_unit <- gsub(" limestone", "",named_rock_unit)
	named_rock_unit <- gsub(" Limestone", "",named_rock_unit)
	named_rock_unit <- gsub(" Limeston", "",named_rock_unit)
	named_rock_unit <- gsub(" Dolomite", "",named_rock_unit)
	named_rock_unit <- gsub(" Dolomites", "",named_rock_unit)
	named_rock_unit <- gsub(" Sandstone", "",named_rock_unit)
	named_rock_unit <- gsub(" Sandstones", "",named_rock_unit)
	named_rock_unit <- gsub(" Shale", "",named_rock_unit)
	named_rock_unit <- gsub(" Shales", "",named_rock_unit)
	named_rock_unit <- gsub(" Marl", "",named_rock_unit)
	named_rock_unit <- gsub(" Marls", "",named_rock_unit)
	}
return(named_rock_unit)
}

scourgify_taxon_names <- function(taxon_name)	{
# taxon_name: string giving species name
taxon_name <- gsub(" n\\. sp\\.","",taxon_name)
taxon_name <- gsub("n\\. gen\\. ","",taxon_name)
taxon_name <- gsub(" n\\. subgen\\.","",taxon_name)
taxon_name <- gsub(" cf\\.","",taxon_name)
taxon_name <- gsub(" aff\\.","",taxon_name)
taxon_name <- gsub(" ex gr\\.","",taxon_name)
taxon_name <- gsub(" informal","",taxon_name)
taxon_name <- gsub(" sensu lato","",taxon_name)
taxon_name <- gsub("\" ", "",taxon_name)
taxon_name <- gsub(" \"", "",taxon_name)
taxon_name <- gsub(" \\?" ,"",taxon_name)
taxon_name <- gsub("\\? " ,"",taxon_name)
taxon_name <- gsub("\\?" ,"",taxon_name)
taxon_name <- gsub("  " ," ",taxon_name)
return(taxon_name)
}