# scripts to download & sort PaleoDB data for the purposes of curation.
#	Concocted by Peter Wagner (pjwagner3@gmail.com)

paleodb_entered_species_for_curation <- function(taxa,onset="Cambrian",end="Holocene",oldest_entry,latest_entry=NA,save_files=TRUE,output_type=".txt") {
# Relevant functions
# function to get occurrence data for a group from a general interval
#	of time and entered between particular dates
# Will return lists, but also can output csv, txt, tab, etc., files directly
# Returns three basic results:
# 1. Genus species combinations with no authority data
# 2. Sources for occurrences of those species
# 3. Genus species combinations for species with authority data
if (is.na(latest_entry))	{
	now <- as.character(Sys.time())
	latest_entry <- strsplit(now," ")[[1]][1]
	}
taxa <- paste(taxa, collapse = ",")
http <- paste("http://paleobiodb.org/data1.2/occs/list.csv?base_name=",taxa,"&interval=",onset,",",end,"&occs_created_after=",oldest_entry,"&occs_created_before=",latest_entry,"&show=coords,paleoloc,refattr,ref&limit=all",sep = "")
fetch <- RCurl::getURL(http)
all_finds <- utils::read.csv(text = fetch, header = TRUE, stringsAsFactors=TRUE)
species_finds <- rbind(subset(all_finds,all_finds$identified_rank=="species"),subset(all_finds,all_finds$identified_rank=="subspecies"))
cleaned_names <- sapply(as.character(species_finds$identified_name),scourgify_occurrence_identifications)
species_finds$identified_name <- cleaned_names

# set aside finds for species with no taxonomic information
unentered_species_finds <- accio_records_of_unentered_species(species_finds)
unentered_species <- sort(unique(unentered_species_finds$identified_name))
us <- length(unentered_species)
# get the references generating the species lacking taxonomic information
unentered_species_sources <- c()
for (i in 1:us)	{
	unentered_species_sources <- 
	rbind(unentered_species_sources,accio_sources_for_taxon_finds(taxon=unentered_species[i],finds=unentered_species_finds))
	}
unentered_epithets <- sapply(unentered_species,accio_species_epithets)
untered_species_finds <- sapply(unentered_species,count_collections_with_taxon,finds=unentered_species_finds)
unentered_species <- cbind(unentered_species,unentered_epithets,untered_species_finds)
colnames(unentered_species) <- c("identified_name","epithet","occs")

ttl_finds <- dim(species_finds)[1]
xx <- (1:ttl_finds)[!species_finds$occurrence_no %in% unentered_species_finds$occurrence_no]
entered_species_finds <- species_finds[xx,]
redone <- entered_species_finds[order(entered_species_finds$accepted_name,entered_species_finds$identified_name),]
entered_species <- unique(cbind(redone$accepted_no,as.character(redone$accepted_name),redone$identified_no,redone$identified_name))
colnames(entered_species) <- c("accepted_no","accepted_name","identified_no","identified_name")
taxon_attr <- sapply(as.numeric(entered_species[,3]),accio_taxon_authors)
entered_species <- cbind(entered_species,taxon_attr)

if (save_files)	{
	if (output_type==".csv")	{
		sepr <- ","
		}	else	{
		sepr <- "\t"
		}
	if (onset!=end)	{
		timespan <- paste(onset,"-",end,sep="")
		}	else	timespan <- onset
	output1 <- paste(timespan,"_",taxa,"_Unentered_Species",output_type,sep="")
	write.table(unentered_species,file=output1,sep = sepr,row.names = FALSE)

	output2 <- paste(timespan,"_",taxa,"_Sources_for_Unentered_Species",output_type,sep="")
	write.table(unentered_species_sources,file=output2,sep = sepr,row.names = FALSE)

	output3 <- paste(timespan,"_",taxa,"_Entered_Species",output_type,sep="")
	write.table(entered_species,file=output3,sep = sepr,row.names = FALSE)

	output4 <- paste(timespan,"_",taxa,"_Occurrences",output_type,sep="")
	write.table(species_finds,file=output4,sep = sepr,row.names = FALSE)
	}

output <- list(unentered_species,unentered_species_sources,entered_species)
names(output) <- c("Unentered_Species","Unentered_Species_Sources","Entered_Species")
return(output)
}

accio_records_of_unentered_species <- function(species_finds)	{
# function to separate records of species that do not have authority data
types_of_diffs <- as.character(unique(species_finds$difference))
tods <- length(types_of_diffs)
unentered <- c()
for (t in 1:tods)	{
	issues <- simplify2array(strsplit(types_of_diffs[t],", "))
	if (!is.na(match("species not entered",issues)))
		unentered <- c(unentered,t)
	}
if (length(unentered)>0)	{
	unentered_species_finds <- subset(species_finds,species_finds$difference==types_of_diffs[unentered[1]])
	if (length(unentered)>1)	{
		for (u in 2:length(unentered))	{
			unentered_species_finds <- rbind(unentered_species_finds,subset(species_finds,species_finds$difference==types_of_diffs[unentered[u]]))	
			}
		}
	}
return(unentered_species_finds)
}

accio_taxon_authors <- function(taxon_no)	{
# function to get authors of taxa with authority data
http <- paste("http://paleobiodb.org/data1.2/taxa/single.csv?id=txn:",taxon_no,"&show=attr",sep="")
GotURL <- RCurl::getURL(http)
taxonomy <- utils::read.csv(text = GotURL, header = TRUE, stringsAsFactors=TRUE)
return(as.character(taxonomy$taxon_attr))
}

accio_collections_with_taxon <- function(taxon,finds)	{
# function to get collections containing a particular taxon
ttl_finds <- dim(finds)[1]
return((1:ttl_finds)[finds$identified_name %in% taxon])
}

count_collections_with_taxon <- function(taxon,finds)	{
# function to get collections containing a particular taxon
return(length((1:dim(finds)[1])[finds$identified_name %in% taxon]))
}

accio_sources_for_taxon_finds <- function(taxon,finds)	{
# function to get original sources for occurrences of a particular taxon
collections <- accio_collections_with_taxon(taxon,finds)
sources <- unique(finds$reference_no[collections])
examples <- match(sources,finds$reference_no)
output <- matrix(0,length(sources),3)
colnames(output) <- c("identified_name","reference_no","primary_reference")
for (i in 1:length(sources))	output[i,] <- c(taxon,sources[i],as.character(finds$primary_reference[examples[i]]))
return(output)
}

accio_species_epithets <- function(binomen)	{
# function to extract the species epithet of a binomial combination
parts <- simplify2array(strsplit(binomen," "))
return(parts[length(parts)])
}

scourgify_occurrence_identifications <- function(taxon_name)	{
# function to remove question marks, cf.s, affs., etc.
text <- substring(taxon_name, seq(1, nchar(taxon_name), 1), seq(1, nchar(taxon_name), 1))
ttlch <- nchar(taxon_name)
if (text[1]=="\"")	{
	taxon_name <- pracma::strcat(text[2:ttlch])
	text <- text[2:ttlch]
	}
ttlch <- nchar(taxon_name)
if (text[ttlch]=="\"")
	taxon_name <- pracma::strcat(text[1:(ttlch-1)])
taxon_name <- gsub(" n\\. sp\\.","",taxon_name)
taxon_name <- gsub("n\\. gen\\. ","",taxon_name)
taxon_name <- gsub(" n\\. subgen\\.","",taxon_name)
taxon_name <- gsub(" cf\\.","",taxon_name)
taxon_name <- gsub("cf\\. ","",taxon_name)
taxon_name <- gsub(" aff\\.","",taxon_name)
taxon_name <- gsub("aff\\. ","",taxon_name)
taxon_name <- gsub(" ex gr\\.","",taxon_name)
taxon_name <- gsub(" informal","",taxon_name)
taxon_name <- gsub(" sensu lato","",taxon_name)
taxon_name <- gsub("\" "," ",taxon_name)
taxon_name <- gsub(" \""," ",taxon_name)
taxon_name <- gsub("\" ","",taxon_name)
taxon_name <- gsub(" \\?" ,"",taxon_name)
taxon_name <- gsub("\\? " ,"",taxon_name)
taxon_name <- gsub("  " ," ",taxon_name)
taxon_name <- gsub("  " ," ",taxon_name)
taxon_name <- gsub("  " ," ",taxon_name)
return(as.character(taxon_name))
}

paleodb_rock_units_for_curation <- function(onset="Cambrian", end="Holocene", taxa=c("Animalia","Protista","Plantae"),system="All",delete_rock_type=TRUE,save_files=TRUE,output_type=".txt")	{
# function for general curation of rock units in PaleoDB.
# Routine returns (and will output as files) three summaries:
#	1.  All stages (local & international) to which a rock unit is assigned
#	2.  All zones (local & international) to which a rock unit is assigned
#	3.  All combinations for rock units alternativel classified as formations
#		or members, or different formation assignments for members
# Input: onset & end for oldest & latest rocks
# taxa: leave this set in order to get basically everything
# system: "all" "marine" or "terrestrial
# delete rock type: if TRUE< then "BLANK Limestone" or "BLAH Sandstone" has
#	lithology removed from name; this make Antelope Valley and Antelope Valley
#	Limestone the same rock unit
# save_files: if true, then output is saved by function with name based on
#	time span & system.
# output_type: ending (e.g., .txt or .csv) for file.
taxa <- paste(taxa, collapse = ",")
http <- paste("http://paleobiodb.org/data1.2/colls/list.csv?base_name=",taxa,"&interval=",onset,",",end,"&show=strat,stratext,refattr,ref,geo",sep="")
fetch <- RCurl::getURL(http)
collections <- utils::read.csv(text = fetch, header = TRUE, stringsAsFactors=TRUE)
ttl_coll <- dim(collections)[1]

if (system!="all" && system!="All")	{
	general_environments <- get_general_environments()
	environments <- c()
	if (!is.na(match("marine",system)) || !is.na(match("Marine",system)))	{
		environments <- c(environments,general_environments$marine)
		}	else if (!is.na(match("terrestrial",system)) || !is.na(match("Terrestrial",system)))
		environments <- general_environments$terrestrial
	relevant <- (1:ttl_coll)[collections$environment %in% environments]
	}	else	{
	relevant <- (1:ttl_coll)
	}
relevant_collections <- collections[relevant,]

# clean up rock unit names
clean_groups <- sapply(as.character(relevant_collections$stratgroup),scourgify_rock_unit_names,delete_rock_type)
relevant_collections$stratgroup <- clean_groups
clean_formations <- sapply(as.character(relevant_collections$formation),scourgify_rock_unit_names,delete_rock_type)
relevant_collections$formation <- clean_formations
clean_members <- sapply(as.character(relevant_collections$member),scourgify_rock_unit_names,delete_rock_type)
relevant_collections$member <- clean_members

rel_coll <- dim(relevant_collections)[1]
formation_members <- c()
formation_members_sep <- matrix(0,rel_coll,2)
for (c in 1:rel_coll)	{
	formation_members_sep[c,] <- c(as.character(relevant_collections$formation[c]),as.character(relevant_collections$member[c]))
	if (relevant_collections$formation[c]!="" && relevant_collections$member[c]!="")	{
		formation_members <- c(formation_members,as.character(paste(relevant_collections$formation[c]," [",relevant_collections$member[c],"]",sep="")))
		} else if (relevant_collections$formation[c]!="")	{
		formation_members <- c(formation_members,as.character(relevant_collections$formation[c]))
		} else if (relevant_collections$member[c]!="")	{
		formation_members <- c(formation_members,as.character(paste("[",relevant_collections$member[c],"]",sep="")))
		}	else if (relevant_collections$stratgroup[c]!="")	{
		ad_hoc <- paste(as.character(relevant_collections$stratgroup[c])," Group Unnamed Formation",sep="")
		formation_members <- c(formation_members,ad_hoc)
		formation_members_sep[c,1] <- as.character(ad_hoc)
		}	else	{
		formation_members <- c(formation_members,"")
		}
#	c <- c+1
	}
relevant_collections <- cbind(relevant_collections,formation_members)
rock_order <- order(formation_members)
formation_members <- sort(unique(formation_members))
formation_members_sep <- unique(formation_members_sep[rock_order,])
rock_units <- length(formation_members)

rock_unit_zones <- c()
rock_unit_stages <- c()
for (r in 1:rock_units)	{
	if (formation_members[r]!="" && formation_members[r]!=" ")	{
		## for some reason, subset command is malfunctioning here....
		cases <- (1:rel_coll)[relevant_collections$formation_members %in% formation_members[r]]
		this_rock <- relevant_collections[cases,]
		zones <- as.character(unique(this_rock$zone[this_rock$zone!=""]))
		if (length(zones)>0)	{
			zn <- length(zones)
			zones <- sort(zones)
			rock_unit_zones <- rbind(rock_unit_zones,cbind(rep(formation_members_sep[r,1],zn),rep(formation_members_sep[r,2],zn),zones))
			}
		stages <- sort(unique(c(as.character(this_rock$early_interval[this_rock$early_interval!=""]),as.character(this_rock$late_interval[this_rock$late_interval!=""]))))
		st <- length(stages)
		if (st>0)	{
			rock_unit_stages <- rbind(rock_unit_stages,cbind(rep(formation_members_sep[r,1],st),rep(formation_members_sep[r,2],st),stages))
			}
		}
#	r <- r+1
	}
colnames(rock_unit_stages) <- c("Formation","Member","Stage")
colnames(rock_unit_zones) <- c("Formation","Member","Zone")
### prepare all combinations of formations & members, zones & strat_units

#accio_latest_rock_unit_assignment(formations=relevant_collections$formation,members=relevant_collections$member,ref_pubyr=relevant_collections$ref_pubyr,reference_no=relevant_collections$reference_no)
confused_rocks <- confunded_rock_assignments(formations=relevant_collections$formation,members=relevant_collections$member,reference_no=relevant_collections$reference_no,ref_pubyr=relevant_collections$ref_pubyr,primary_reference=relevant_collections$primary_reference)

if (save_files)	{
	if (output_type==".csv")	{
		sepr <- ","
		}	else	{
		sepr <- "\t"
		}
	if (onset!=end)	{
		timespan <- paste(onset,"-",end,sep="")
		}	else	timespan <- onset
	output1 <- paste(timespan,"_Stage_Assignments_for_",system,"_Strata",output_type,sep="")
	write.table(rock_unit_stages,file=output1,sep = sepr,row.names = FALSE,col.names=TRUE)

	output2 <- paste(timespan,"_Zone_Assignments_for_",system,"_Strata",output_type,sep="")
	write.table(rock_unit_zones,file=output2,sep = sepr,row.names = FALSE,col.names=TRUE)

	output3 <- paste(timespan,"_",system,"_Strata_Entered_as_Formations_and_Members",output_type,sep="")
	write.table(confused_rocks,file=output3,sep = sepr,row.names = FALSE,col.names=TRUE)
	
	output4 <- paste(timespan,"_",system,"_Relevant_Collections",output_type,sep="")
	write.table(relevant_collections,file=output4,sep = sepr,row.names = FALSE,col.names=TRUE)
	}
output <- list(output1,output2,output3)
names(output) <- c("Stage_Assignments","Zone Assignments","Inconsistent_Ranks")
#print("That's all folks....")
return(output)
}

get_general_environments <- function()	{
# function to list particular environments that are either marine or terrestrial
marine <- c("marine indet.","Carbonate marine","carbonate indet.","peritidal","shallow subtidal indet.","open shallow subtidal","lagoonal/restricted shallow subtidal","sand shoal","reef, buildup or bioherm","perireef or subreef","intrashelf/intraplatform reef","platform/shelf-margin reef","slope/ramp reef","basin reef","deep subtidal ramp","deep subtidal shelf","deep subtidal indet.","offshore ramp","offshore shelf","offshore indet.","slope","basinal (carbonate)","basinal (siliceous)","Siliciclastic marine","marginal marine indet.","coastal indet.","estuary/bay","lagoonal","paralic indet.","delta plain","interdistributary bay","delta front","prodelta","deltaic indet.","foreshore","shoreface","transition zone/lower shoreface","offshore","coastal indet.","submarine fan","basinal (siliciclastic)","basinal (siliceous)","basinal (carbonate)","deep-water indet.")
terrestrial <- c("terrestrial indet.","fluvial indet.","alluvial fan","channel lag","coarse channel fill","fine channel fill","\"channel\"","channel","wet floodplain","dry floodplain","\"floodplain\"","floodplain","crevasse splay","levee","mire/swamp","fluvial-lacustrine indet.","delta plain","fluvial-deltaic indet.","lacustrine - large","lacustrine - small","pond","crater lake","lacustrine delta plain","lacustrine interdistributary bay","lacustrine delta front","lacustrine prodelta","lacustrine deltaic indet.","lacustrine indet.")
output <- list(marine,terrestrial)
names(output) <- c("marine","terrestrial")
return(output)
}

scourgify_rock_unit_names <- function(named_rock_unit,delete_rock_type=TRUE)	{
# This removes quotes, accented letters (or their corruptions), lithology addenda
# 	etc. from the rock unit name.
ttlch <- nchar(named_rock_unit)
if (ttlch>0 && named_rock_unit != " ")	{
	text <- substring(named_rock_unit, seq(1, nchar(named_rock_unit), 1), seq(1, nchar(named_rock_unit), 1))
	while (text[1]==" " && ttlch>1)	{
		named_rock_unit <- pracma::strcat(text[2:ttlch])
		text <- text[2:ttlch]
		ttlch <- ttlch-1
		}
	if (text[1]=="\"")	{
		named_rock_unit <- pracma::strcat(text[2:ttlch])
		text <- text[2:ttlch]
		}
	ttlch <- nchar(named_rock_unit)
	if (text[ttlch]=="\"")
		named_rock_unit <- pracma::strcat(text[1:(ttlch-1)])
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
	named_rock_unit <- gsub("Ã¥","a",named_rock_unit)
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
		named_rock_unit <- gsub(" Sandstone", "",named_rock_unit)
		named_rock_unit <- gsub(" Shale", "",named_rock_unit)
		}
	return(named_rock_unit)
	}	else	{
	named_rock_unit <- ""
	return(named_rock_unit)
	}
}

accio_latest_rock_unit_assignment <- function(formations,members,ref_pubyr,reference_no)	{
# function to get the latest formation/member combination for rocks with
#	multiple such combinations and/or rankings in the PaleoDB.
member_or_formation <- members[(1:length(members))[members %in% formations]]
member_or_formation <- member_or_formation[member_or_formation!=""]
ttl_coll <- length(members)
last_authors <- vector(length=ttl_coll)
if (length(member_or_formation)>0)	{
	for (c in 1:length(member_or_formation))	{
		# Use latest opinion.  If latest opinion is "member," then reassign
		#	all collections to the latest formation/member combo
		# If the latest opinion is tied, then go with majority rule.  If that
		#	is tied, too, then just make the damned thing a formation....
		vote_formation <- (1:ttl_coll)[formations %in% member_or_formation[c]]
		vote_member <- (1:ttl_coll)[members %in% member_or_formation[c]]
		if (length(vote_member)>0 && length(vote_formation)>0)	{
			if (max(ref_pubyr[vote_formation]) > max(ref_pubyr[vote_member]) || (max(ref_pubyr[vote_formation]) == max(ref_pubyr[vote_member]) && length(vote_formation) > length(vote_member)))	{
				### elevate member to formation in appropriate collections
				formations[vote_member] <- member_or_formation[c]
				members[vote_member] <- ""
				latest_opinion <- vote_formation[match(max(ref_pubyr[vote_formation]),ref_pubyr[vote_formation])]
				last_authors[c(vote_formation,vote_member)] <- reference_no[latest_opinion]
				}	else if ((max(ref_pubyr[vote_formation]) < max(ref_pubyr[vote_member])) || (max(ref_pubyr[vote_formation]) == max(ref_pubyr[vote_member]) && length(vote_formation) < length(vote_member)))	{
				latest_opinion <- vote_member[match(max(ref_pubyr[vote_member]),ref_pubyr[vote_member])]
				formations[vote_formation] <- formations[latest_opinion]
				members[vote_formation] <- member_or_formation[c]
				last_authors[c(vote_formation,vote_member)] <- reference_no[latest_opinion]
				} else if (sum(vote_formation==vote_member)==length(vote_member))	{
				members[vote_member] <- ""
				} else	{
				latest_opinion <- vote_member[match(max(ref_pubyr[vote_member]),ref_pubyr[vote_member])]
				formations[vote_formation] <- formations[latest_opinion]
				members[vote_formation] <- member_or_formation[c]
				last_authors[c(vote_formation,vote_member)] <- reference_no[latest_opinion]
				}
			}
#		c <- c+1
		}
	}
output <- cbind(formations,members)
colnames(output) <- c("Formation_Latest","Member_Latest")
return(output)
}

confunded_rock_assignments <- function(formations,members,reference_no,ref_pubyr,primary_reference)	{
# function to separate out rock units that seem to have both formation and
#	member status at different localities.  This returns the different references
#	expressing or using the opinions in order for someone to decide what is
#	the best combination to use now.  
member_or_formation <- members[(1:length(members))[members %in% formations]]
member_or_formation <- member_or_formation[member_or_formation!=""]
ttl_coll <- length(members)
confunded <- c()
if (length(member_or_formation)>0)	{
	confused <- sort(unique(member_or_formation))
	confused <- confused[nchar(confused)>1]
	for (c in 1:length(confused))	{
		# Use latest opinion.  If latest opinion is "member," then reassign
		#	all collections to the latest formation/member combo
		# If the latest opinion is tied, then go with majority rule.  If that
		#	is tied, too, then just make the damned thing a formation....
		vote_formation <- (1:ttl_coll)[formations %in% confused[c]]
		vote_member <- (1:ttl_coll)[members %in% confused[c]]
		rocks <- c(vote_formation,vote_member)
		all_refs <- reference_no[rocks]
		uniq_refs <-unique(all_refs)
		relv_colls <- rocks[match(uniq_refs,all_refs)]
		confunded <- rbind(confunded,cbind(formations[relv_colls],members[relv_colls],ref_pubyr[relv_colls],reference_no[relv_colls],as.character(primary_reference[relv_colls])))
		}
	}
colnames(confunded) <- c("formation","member","ref_pubyr","reference_no","primary_reference")
return(confunded)
}