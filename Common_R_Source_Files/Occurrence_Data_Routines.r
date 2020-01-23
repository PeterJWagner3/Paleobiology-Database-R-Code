#source('~/Documents/R_Projects/Common_R_Source_Files/Wagner_kluges.r')  #
#dev.off()
# get taxon data
#species_info <- read.table(file="Species_Data_More.tab", header=FALSE, stringsAsFactors=FALSE, sep="\t")
#occurrences <- read.table(file="Occurrences_More.tab", header=FALSE, stringsAsFactors=FALSE, sep="\t")
library(stringr)

##### OCCURRENCE & COLLECTION DATA CLEANING #####
evanesco_indeterminate_species <- function(paleodb_finds)	{
taxon_name <- sort(unique(paleodb_finds$accepted_name));
species_epithets <- sapply(taxon_name,diffindo_species_epithets);
indet_species <- c("sp.");
indet_species <- c(indet_species,paste("sp.",LETTERS));
indet_species <- c(indet_species,paste("sp.",letters));
for (i in 1:100)	indet_species <- c(indet_species,paste("sp.",i));
for (i in 1:100)	indet_species <- c(indet_species,paste("nov.",i));
for (i in 1:100)	indet_species <- c(indet_species,paste("sp. nov.",i));
for (i in 1:100)	indet_species <- c(indet_species,paste("indet.",i));
species_epithets <- sapply(species_epithets,accio_embedded_informal_names);
taxon_names <- taxon_name[!species_epithets %in% indet_species];
echino_species <- taxon_names <- unique(taxon_names);
taxon_names <- sapply(echino_species,echinoscrub);
taxon_names <- taxon_names[taxon_names!=""];
paleodb_finds <- paleodb_finds[paleodb_finds$accepted_name %in% taxon_names,];
return(paleodb_finds);
}

accio_embedded_informal_names <- function(species_epithets)	{
for (se in 1:length(species_epithets))	{
	if (length(strsplit(species_epithets[se],split=" ")[[1]])>1)	{
		if (sum(strsplit(species_epithets[se],split=" ")[[1]] %in% c(0:9,letters,LETTERS,"genus","sp.","indet."))>0)
			species_epithets[se] <- "sp. 1";
		}
	}
return(species_epithets)
}

echinoscrub <- function(echino_species)	{
echinobabble <- c("columnals","debris","stem","stems","holdfast","holdfasts","miscellanea","miscellaneus","ossicle","ossicles","plate","plates");
j <- simplify2array(strsplit(echino_species," ")[[1]]);
if (sum(tolower(j) %in% echinobabble)>0)	{
	return("");
	} else	{
	return(echino_species);
	}
}

# count occurrences by subintervals, with finds attributable only to general intervals spread among subintervals
tally_occurrences_by_subinterval <- function(coll_no,early_interval,late_interval,hierarchical_chronostrat,lump_cooccr=T,constrain=F,temporal_precision=0.1)	{
# coll_no: colletion number
# early_interval: earliest possible chronostratigraphic interval for corresponding coll_no
# late_interval: lateest possible chronostratigraphic interval for corresponding coll_no
# hierarchical_chronostrat: table denoting which chronostratigraphic units are subunits of others
# lump_cooccr: it true, the co-occurrences in the same collection are lumped into 1 find
# constrain: ???	

if (constrain)	{
	fa_latest <- min(unique(hierarchical_chronostrat$bin_last[match(late_interval,hierarchical_chronostrat$interval)]));
	la_earliest <- max(unique(hierarchical_chronostrat$bin_first[match(early_interval,hierarchical_chronostrat$interval)]));
	} else	{
	fa_latest <- min(hierarchical_chronostrat$bin_first);
	la_earliest <- max(hierarchical_chronostrat$bin_last);
	}
	
ttl_finds <- length(coll_no);
bins_early <- hierarchical_chronostrat$bin_first[match(early_interval,hierarchical_chronostrat$interval)];
bins_late <- hierarchical_chronostrat$bin_last[match(late_interval,hierarchical_chronostrat$interval)];
bins_early[bins_early<fa_latest] <- fa_latest;
bins_late[bins_late>la_earliest] <- la_earliest;
span <- 1+abs(bins_late-bins_early);
bin_finds <- vector(length=max(hierarchical_chronostrat$bin_last));
for (f in 1:ttl_finds)	{
	bin_finds[bins_early[f]:bins_late[f]] <- bin_finds[bins_early[f]:bins_late[f]]+1/span[f];
	}
return(bin_finds);
}

# count rock-units occupied by subintervals, with possible fractional counts
tally_collections_occupied_by_subinterval <- function(taxon_collections,hierarchical_chronostrat,constrain=F,temporal_precision=0.1)	{
# coll_no: collection number
# early_interval: earliest possible chronostratigraphic interval for corresponding coll_no
# late_interval: lateest possible chronostratigraphic interval for corresponding coll_no
# hierarchical_chronostrat: table denoting which chronostratigraphic units are subunits of others
# constrain: return only relevant time intervals; otherwise, return results for entire time scale.
finest_chronostrat <- subset(hierarchical_chronostrat,hierarchical_chronostrat$bin_first==hierarchical_chronostrat$bin_last);
finest_chronostrat <- finest_chronostrat[unique(match(finest_chronostrat$bin_first,finest_chronostrat$bin_first)),];
if (constrain)	{
	early_interval <- as.character(taxon_collections$interval_lb);
	late_interval <- as.character(taxon_collections$interval_ub);
	fa_latest <- min(unique(hierarchical_chronostrat$bin_last[match(late_interval,hierarchical_chronostrat$interval)]));
	la_earliest <- max(unique(hierarchical_chronostrat$bin_first[match(early_interval,hierarchical_chronostrat$interval)]));
	} else	{
	fa_latest <- min(hierarchical_chronostrat$bin_first);
	la_earliest <- max(hierarchical_chronostrat$bin_last);
	}
collection_finds <- data.frame(array(0,dim=c(nrow(taxon_collections),1+(la_earliest-fa_latest))));
colnames(collection_finds) <- finest_chronostrat$interval;
for (cn in 1:nrow(taxon_collections))	{
	relevant_collections <- taxon_collections[cn,];
	collection_finds[cn,] <- count_units_per_bin_fuzzily(relevant_collections,finest_chronostrat = finest_chronostrat,temporal_precision = temporal_precision);
	}
return(colSums(collection_finds));
}

# count rock-units occupied by subintervals, with possible fractional counts
#tally_rock_units_occupied_by_subinterval <- function(taxon_collections,hierarchical_chronostrat,constrain=F,temporal_precision=0.1)	{
tally_rock_units_occupied_by_subinterval <- function(taxon_collections,finest_chronostrat,constrain=F,temporal_precision=0.1)	{
# coll_no: collection number
# early_interval: earliest possible chronostratigraphic interval for corresponding coll_no
# late_interval: lateest possible chronostratigraphic interval for corresponding coll_no
# hierarchical_chronostrat: table denoting which chronostratigraphic units are subunits of others
# constrain: return only relevant time intervals; otherwise, return results for entire time scale.
#finest_chronostrat <- subset(hierarchical_chronostrat,hierarchical_chronostrat$bin_first==hierarchical_chronostrat$bin_last);
#finest_chronostrat <- finest_chronostrat[unique(match(finest_chronostrat$bin_first,finest_chronostrat$bin_first)),];
if (constrain)	{
	early_interval <- as.character(taxon_collections$interval_lb);
	late_interval <- as.character(taxon_collections$interval_ub);
	fa_latest <- min(unique(finest_chronostrat$bin_last[match(late_interval,finest_chronostrat$interval)]));
	la_earliest <- max(unique(finest_chronostrat$bin_first[match(early_interval,finest_chronostrat$interval)]));
#	fa_latest <- min(unique(hierarchical_chronostrat$bin_last[match(late_interval,hierarchical_chronostrat$interval)]));
#	la_earliest <- max(unique(hierarchical_chronostrat$bin_first[match(early_interval,hierarchical_chronostrat$interval)]));
	} else	{
	fa_latest <- min(finest_chronostrat$bin_first);
	la_earliest <- max(finest_chronostrat$bin_last);
#	fa_latest <- min(hierarchical_chronostrat$bin_first);
#	la_earliest <- max(hierarchical_chronostrat$bin_last);
	}
taxon_collections <- name_unnamed_rock_units(paleodb_collections=taxon_collections,finest_chronostrat);
unique_rocks <- sort(unique(taxon_collections$rock_no_sr[taxon_collections$rock_no_sr>0]));
rock_finds <- data.frame(array(0,dim=c(length(unique_rocks),1+(la_earliest-fa_latest))));
colnames(rock_finds) <- finest_chronostrat$interval;
# get rock-units
rn <- 0;
while (rn < length(unique_rocks))	{
	rn <- rn+1;
	relevant_collections <- subset(taxon_collections,taxon_collections$rock_no_sr==unique_rocks[rn]);
	rock_finds[rn,] <- count_units_per_bin_fuzzily(relevant_collections,finest_chronostrat = finest_chronostrat,temporal_precision = temporal_precision);
	}
return(colSums(rock_finds))
}

count_units_per_bin_fuzzily <- function(relevant_collections,finest_chronostrat,temporal_precision=0.1)	{
#	prob_find_this_rock_bin <- array(0,dim=c(nrow(this_rock),length(finest_chronostrat$interval)));
prob_find <- array(0,dim=c(length(finest_chronostrat$interval)));
prob_find_this_set <- c();
round_level <- ceiling(-log10(temporal_precision));
r_c <- 0;
while (r_c < nrow(relevant_collections))	{
	r_c <- r_c+1;
	prob_find_this_case <- array(0,dim=c(length(finest_chronostrat$interval)));
	aa <- min(finest_chronostrat$bin_first[match(relevant_collections$interval_lb[r_c],finest_chronostrat$interval)]);
	zz <- max(finest_chronostrat$bin_last[match(relevant_collections$interval_ub[r_c],finest_chronostrat$interval)]);
	if (aa==zz)	{
		prob_find_this_case[aa] <- 1;
		} else	{
		if (relevant_collections$ma_lb[r_c]!=relevant_collections$ma_ub[r_c])	{
			ma_range <- abs(relevant_collections$ma_lb[r_c]-relevant_collections$ma_ub[r_c]);
			range_start <- relevant_collections$ma_lb[r_c];
			range_end <- relevant_collections$ma_ub[r_c];
			} else	{
			ma_range <- abs(mean(finest_chronostrat$ma_lb[aa],finest_chronostrat$ma_ub[aa])-mean(finest_chronostrat$ma_lb[zz],finest_chronostrat$ma_ub[zz]));
			range_start <- mean(finest_chronostrat$ma_lb[aa],finest_chronostrat$ma_ub[aa]);
			range_end <- mean(finest_chronostrat$ma_lb[zz],finest_chronostrat$ma_ub[zz]);
			}
		ma_range_finds <- array(1/(ma_range/temporal_precision),dim=c(ma_range/temporal_precision));
		ma_range_subbin_finds <- temporal_precision/ma_range;
		ma_range_find_ages <- round(seq(range_start-temporal_precision,range_end,by=-temporal_precision),round_level);
		for (bn in aa:zz)	{
			prob_bin <- ma_range_subbin_finds*sum(ma_range_find_ages[ma_range_find_ages<round(finest_chronostrat$ma_lb[bn],round_level)] %in% ma_range_find_ages[ma_range_find_ages>=round(finest_chronostrat$ma_ub[bn],round_level)]);
		#		rock_finds[rn,bn] <- max(rock_finds[rn,bn],prob_bin);
			prob_find_this_case[bn] <- round(prob_bin,round_level);
			}
		}
	prob_find_this_set <- rbind(prob_find_this_set,prob_find_this_case);
#	prob_find_this_set;
	}
if (length(prob_find_this_set)==nrow(finest_chronostrat))	{
	prob_find <- 1-exp(log(1-prob_find_this_set));
	}	else	{
	prob_find <- 1-exp(colSums(log(1-prob_find_this_set)));
	}
return(prob_find);
}

# count occurrences by subintervals, with finds attributable only to general intervals spread among subintervals
tally_definite_occurrences_by_subinterval <- function(coll_no,early_interval,late_interval,hierarchical_chronostrat,lump_cooccr=T)	{
# coll_no: colletion number
# early_interval: earliest possible chronostratigraphic interval for corresponding coll_no
# late_interval: lateest possible chronostratigraphic interval for corresponding coll_no
# hierarchical_chronostrat: table denoting which chronostratigraphic units are subunits of others
# lump_cooccr: it true, the co-occurrences in the same collection are lumped into 1 find
if (lump_cooccr)	{
	orig_coll_no <- coll_no;
	coll_no <- unique(coll_no);
	early_interval <- early_interval[match(coll_no,orig_coll_no)];
	late_interval <- late_interval[match(coll_no,orig_coll_no)];
	}
ttl_finds <- length(coll_no);
bins_early <- hierarchical_chronostrat$bin_first[match(early_interval,hierarchical_chronostrat$interval)];
bins_late <- hierarchical_chronostrat$bin_last[match(late_interval,hierarchical_chronostrat$interval)];
span <- 1+abs(bins_late-bins_early);
bin_finds <- vector(length=max(hierarchical_chronostrat$bin_last));
for (f in 1:ttl_finds)
	if (bins_early[f]==bins_late[f])
		bin_finds[bins_early[f]] <- bin_finds[bins_early[f]]+1;
return(bin_finds);
}

# count rock-units occupied by subintervals, with possible fractional counts
tally_rock_units_definitely_occupied_by_subinterval <- function(taxon_collections,hierarchical_chronostrat,constrain=F)	{
# coll_no: collection number
# early_interval: earliest possible chronostratigraphic interval for corresponding coll_no
# late_interval: lateest possible chronostratigraphic interval for corresponding coll_no
# hierarchical_chronostrat: table denoting which chronostratigraphic units are subunits of others
# constrain: return only relevant time intervals; otherwise, return results for entire time scale.
finest_chronostrat <- subset(hierarchical_chronostrat,hierarchical_chronostrat$bin_first==hierarchical_chronostrat$bin_last);
if (constrain)	{
	early_interval <- as.character(taxon_collections$interval_lb);
	late_interval <- as.character(taxon_collections$interval_ub);
	fa_latest <- min(unique(hierarchical_chronostrat$bin_last[match(late_interval,hierarchical_chronostrat$interval)]));
	la_earliest <- max(unique(hierarchical_chronostrat$bin_first[match(early_interval,hierarchical_chronostrat$interval)]));
	} else	{
	fa_latest <- min(hierarchical_chronostrat$bin_first);
	la_earliest <- max(hierarchical_chronostrat$bin_last);
	}
unique_rocks <- sort(unique(taxon_collections$rock_no_sr[taxon_collections$rock_no_sr>0]));
rock_finds <- data.frame(array(0,dim=c(length(unique_rocks),1+(la_earliest-fa_latest))));
rownames(rock_finds) <- unique_rocks;
colnames(rock_finds) <- finest_chronostrat$interval;
# get rock-units
# this needs to be inside time loop!
#tx_colls <- nrow(taxon_finds);
set_collections <- subset(taxon_collections,taxon_collections$interval_lb==taxon_collections$interval_ub);
unset_collections <- subset(taxon_collections,taxon_collections$interval_lb!=taxon_collections$interval_ub);
for (rn in 1:length(unique_rocks))	{
	this_rock <- subset(set_collections,set_collections$rock_no_sr==unique_rocks[rn]);
	if (nrow(this_rock)>0)	{
		this_rock_bins <- unique(hierarchical_chronostrat$bin_first[match(this_rock$interval_lb,hierarchical_chronostrat$interval)])
		rock_finds[rn,this_rock_bins] <- 1;
		}
	}
return(colSums(rock_finds[fa_latest:la_earliest]));
}

# count occurrences by subintervals, with finds attributable only to general intervals spread among subintervals
tally_presence_by_subinterval <- function(coll_no,early_interval,late_interval,hierarchical_chronostrat,lump_cooccr=T)	{
# coll_no: colletion number
# early_interval: earliest possible chronostratigraphic interval for corresponding coll_no
# late_interval: lateest possible chronostratigraphic interval for corresponding coll_no
# hierarchical_chronostrat: table denoting which chronostratigraphic units are subunits of others
# lump_cooccr: it true, the co-occurrences in the same collection are lumped into 1 find
if (lump_cooccr)	{
	orig_coll_no <- coll_no;
	coll_no <- unique(coll_no);
	early_interval <- early_interval[match(coll_no,orig_coll_no)];
	late_interval <- late_interval[match(coll_no,orig_coll_no)];
	}
ttl_finds <- length(coll_no);
bins_early <- hierarchical_chronostrat$bin_first[match(early_interval,hierarchical_chronostrat$interval)];
bins_late <- hierarchical_chronostrat$bin_last[match(late_interval,hierarchical_chronostrat$interval)];
span <- 1+abs(bins_late-bins_early);
bin_finds <- vector(length=max(hierarchical_chronostrat$bin_last));
for (f in 1:ttl_finds)
	if (bins_early[f]==bins_late[f])
		bin_finds[bins_early[f]] <- 1;
return(bin_finds);
}

# count occurrences by subintervals, with finds attributable only to general intervals spread among subintervals
tally_presence_by_subinterval_minimum <- function(early_interval,late_interval,hierarchical_chronostrat,lump_cooccr=T,constrain=F)	{
# coll_no: colletion number
# early_interval: earliest possible chronostratigraphic interval for corresponding coll_no
# late_interval: lateest possible chronostratigraphic interval for corresponding coll_no
# hierarchical_chronostrat: table denoting which chronostratigraphic units are subunits of others
# lump_cooccr: it true, the co-occurrences in the same collection are lumped into 1 find

fa_latest <- min(unique(hierarchical_chronostrat$bin_last[match(late_interval,hierarchical_chronostrat$interval)]));
la_earliest <- max(unique(hierarchical_chronostrat$bin_first[match(early_interval,hierarchical_chronostrat$interval)]));

ttl_finds <- length(coll_no);
bins_early <- hierarchical_chronostrat$bin_first[match(early_interval,hierarchical_chronostrat$interval)];
bins_late <- hierarchical_chronostrat$bin_last[match(late_interval,hierarchical_chronostrat$interval)];
bins_early[bins_early<fa_latest] <- fa_latest;
bins_late[bins_late>la_earliest] <- la_earliest;
span <- 1+abs(bins_late-bins_early);
bin_finds <- vector(length=max(hierarchical_chronostrat$bin_last));
for (f in 1:ttl_finds)	{
	bin_finds[bins_early[f]:bins_late[f]] <- bin_finds[bins_early[f]:bins_late[f]]+1/span[f];
	}
return(bin_finds);
}

evanesco_duplicate_occurrences <- function(occurrences,species_info)	{
ttl_finds_cnd <- ttl_finds <- dim(occurrences)[1]
species <- sort(unique(species_info$OTU_Sr))
otu <- length(species)
deletions <- 0
for (s in 1:otu)	{
	spc <- species[s]
	species_finds <- subset(occurrences,occurrences$OTU_Sr==spc)
	duplicates <- length(species_finds$Locality)-length(unique(species_finds$Locality))
	if (duplicates>0)	{
		seniors <- match(unique(species_finds$Locality),species_finds$Locality)
		row_nos <- as.numeric(rownames(species_finds))
		suspects <- length(row_nos)
		for (i in suspects:1)	{
			if (is.na(match(i,seniors)))	{
#				rn <- row_nos[i]
#				dummy <- rbind(occurrences[1:(rn-1),],occurrences[(rn+1):ttl_finds,])
				if (deletions==0)	{
					rows_to_cut <- row_nos[i]
					}	else {
					rows_to_cut <- c(rows_to_cut,row_nos[i])
					}
				deletions <- deletions+1
				}
			}
		}
	}

rows_to_cut <- sort(rows_to_cut,decreasing=FALSE)
rows_to_keep <- (1:ttl_finds)[!(1:ttl_finds) %in% rows_to_cut]
return(occurrences[rows_to_keep,])
}

evanesco_duplicate_occurrences_paleodb <- function(occurrences)	{
ttl_finds_cnd <- ttl_finds <- nrow(occurrences);
species <- sort(unique(occurrences$accepted_name));
otu <- length(species)

reduced_finds <- data.frame(collection_no=as.numeric(occurrences$collection_no),
							accepted_name=as.character(occurrences$accepted_name));
reduced_finds <- unique(reduced_finds);
keepers <- as.numeric(rownames(reduced_finds));
return(occurrences[keepers,]);
}

# assume that co-occuring congeneric singleton species are a single species.
lump_cooccuring_singletons <- function(occurrences,species_info)	{
species <- sort(unique(species_info$OTU_Sr))
genera <- sort(unique(species_info$Genus))
ttl_genera <- length(genera)

colnames(occurrences) <- c("Locality", "OTU", "OTU_Sr");
otus <- length(species)
localities <- unique(occurrences$Locality)
ttl_localities <- length(localities)

species_summary <- matrix(0,otus,3)
colnames(species_summary) <- c("Spc No", "Gen No", "Finds")
senior_species_info <- subset(species_info,species_info$OTU==species_info$OTU_Sr)
for (s in 1:otus)	{
	species_summary[s,1] <- senior_species_info$OTU_Sr[s]
	species_summary[s,2] <- match(as.character(senior_species_info$Genus[s]),genera)
	species_summary[s,3] <- length(unique(subset(occurrences,occurrences$OTU_Sr==species_summary[s,1])[,1]))
	}

singletons <- subset(species_summary,species_summary[,3]==1)
ttl_singletons <- dim(singletons)[1]

# add genus numbers to occurrence information
genus_nos <- sort(unique(species_summary[,2]))
ttl_finds <- dim(occurrences)[1]
GTU <- vector(length=ttl_finds)
for (i in 1:ttl_finds)	{
	s <- match(occurrences[i,3],species_summary[,1])
	GTU[i] <- species_summary[s,2]
	}
occurrences <- cbind(occurrences,GTU)
#write.table(occurrences,"Look for Problems.xls",row.names=FALSE,sep="\t")

genus_singletons_per_locality <- matrix(0,ttl_genera,ttl_localities)
for (s in 1:ttl_singletons)	{
	spc <- singletons[s,1]
	gen <- match(singletons[s,2],genus_nos)
	oc <- match(spc,occurrences$OTU_Sr)
	lc <- match(occurrences$Locality[oc],localities)
	genus_singletons_per_locality[gen,lc] <- genus_singletons_per_locality[gen,lc]+1
	}

genus_species_per_locality <- matrix(0,ttl_genera,ttl_localities)
genus_finds_all <- genus_finds_unq <- ttl_genera
for (g in 1:ttl_genera)	{
	gen <- genus_nos[g]
	genus_finds <- subset(occurrences,occurrences$GTU==gen)
	genus_finds_all[g] <- length(genus_finds$Locality)
	if (genus_finds_all[g]>1)	{
		unique_locs <- unique(genus_finds$Locality)
		genus_finds_unq[g] <- length(unique(unique_locs))
		for (i in 1:genus_finds_unq[g])	{
			loc <- match(unique_locs[i],localities)
			genus_species_per_locality[g,loc] <- sum(genus_finds$Locality==unique_locs[i])
			}
		}	else {
		genus_finds_unq[g] <- genus_finds_all[g]
		loc <- match(genus_finds[1],localities)
		genus_species_per_locality[g,loc] <- 1
		}
	}
gen_ord <- order(singletons[,2])
singletons <- singletons[gen_ord,]
done <- vector(length=dim(singletons)[1])
singletons <- cbind(singletons,done)
genera_w_singletons <- unique(singletons[,2])
genus_singletons <- vector(length=length(genera_w_singletons))
names(genus_singletons) <- genera[genera_w_singletons]
for (g in 1:length(genera_w_singletons))	{
	genus_singletons[g] <- sum(singletons[,2]==genera_w_singletons[g])
	}
prob_genera <- sort(genus_singletons,decreasing=TRUE)
lumped <- vector(length=length(prob_genera))
names(lumped) <- names(prob_genera)
poss_cases <- sum(prob_genera>1)
for (g in 1:poss_cases)	{
	prob_child <- names(prob_genera)[g]
	gen_no <- match(prob_child,genera)
	relevant_finds <- subset(occurrences,occurrences$GTU==gen_no)
	singles <- subset(singletons,singletons[,2]==gen_no)
	lumped[g] <- 0
	
	for (s in 1:prob_genera[g])	{
		while (singles[s,4]==1 && s<prob_genera[g])	s <- s+1
		loner <- singles[s,1]
		oc_no <- match(loner,relevant_finds$OTU_Sr)
		loc <- relevant_finds$Locality[oc_no]
		loner_congen <- subset(relevant_finds,relevant_finds$Locality==loc)
		if (length(loner_congen$OTU_Sr)>1)	{
			# put on # finds to identify other singletons
			test <- cbind(loner_congen$OTU_Sr,species_summary[match(loner_congen$OTU_Sr,species_summary[,1]),3])
			others <- subset(test,test[,2]==1)
			if (length(others[,1])>1)	{
				# "synonymize" species
				change <- match(others[,1],occurrences$OTU_Sr)
				occurrences$OTU_Sr[change] <- loner
				# check other singletons off as done
				redone <- match(others[,1],singles[,1])
				singles[redone,4] <- 1
				# add to sum of lumped
				lumped[g] <- lumped[g]+(length(redone)-1)
				}
			}
		}
	}

# need to clear occurrences here.  Delete duplicate finds for senior species
cleaned_occurrences <- evanesco_duplicate_occurrences(occurrences,species_info)
return(cleaned_occurrences[,1:3])
}

lump_singletons_w_multifind_congenerics <- function(occurrences,species_info)	{
colnames(species_info) <- c("OTU","OTU_Sr","Genus","species","Family")
species <- sort(unique(species_info$OTU_Sr))
genera <- sort(unique(species_info$Genus))
ttl_genera <- length(genera)

colnames(occurrences) <- c("Locality", "OTU", "OTU_Sr")
otus <- length(species)
localities <- unique(occurrences$Locality)
ttl_localities <- length(localities)

species_summary <- matrix(0,otus,3)
colnames(species_summary) <- c("Spc No", "Gen No", "Finds")
senior_species_info <- subset(species_info,species_info$OTU==species_info$OTU_Sr)
for (s in 1:otus)	{
	species_summary[s,1] <- senior_species_info$OTU_Sr[s]
	species_summary[s,2] <- match(as.character(senior_species_info$Genus[s]),genera)
	species_summary[s,3] <- length(unique(subset(occurrences,occurrences$OTU_Sr==species_summary[s,1])[,1]))
	}

singletons <- subset(species_summary,species_summary[,3]==1)
ttl_singletons <- dim(singletons)[1]

# add genus numbers to occurrence information
genus_nos <- sort(unique(species_summary[,2]))
ttl_finds <- dim(occurrences)[1]
GTU <- vector(length=ttl_finds)
for (i in 1:ttl_finds)	{
	s <- match(occurrences[i,3],species_summary[,1])
	GTU[i] <- species_summary[s,2]
	}
occurrences <- cbind(occurrences,GTU)
#write.table(occurrences,"Look for Problems.xls",row.names=FALSE,sep="\t")

genus_singletons_per_locality <- matrix(0,ttl_genera,ttl_localities)
for (s in 1:ttl_singletons)	{
	spc <- singletons[s,1]
	gen <- match(singletons[s,2],genus_nos)
	oc <- match(spc,occurrences$OTU_Sr)
	lc <- match(occurrences$Locality[oc],localities)
	genus_singletons_per_locality[gen,lc] <- genus_singletons_per_locality[gen,lc]+1
	}

genus_species_per_locality <- matrix(0,ttl_genera,ttl_localities)
genus_finds_all <- genus_finds_unq <- ttl_genera
for (g in 1:ttl_genera)	{
	gen <- genus_nos[g]
	genus_finds <- subset(occurrences,occurrences$GTU==gen)
	genus_finds_all[g] <- length(genus_finds$Locality)
	if (genus_finds_all[g]>1)	{
		unique_locs <- unique(genus_finds$Locality)
		genus_finds_unq[g] <- length(unique(unique_locs))
		for (i in 1:genus_finds_unq[g])	{
			loc <- match(unique_locs[i],localities)
			genus_species_per_locality[g,loc] <- sum(genus_finds$Locality==unique_locs[i])
			}
		}	else {
		genus_finds_unq[g] <- genus_finds_all[g]
		loc <- match(genus_finds[1],localities)
		genus_species_per_locality[g,loc] <- 1
		}
	}
gen_ord <- order(singletons[,2])
singletons <- singletons[gen_ord,]
done <- vector(length=dim(singletons)[1])
singletons <- cbind(singletons,done)
genera_w_singletons <- unique(singletons[,2])
genus_singletons <- vector(length=length(genera_w_singletons))
names(genus_singletons) <- genera[genera_w_singletons]
for (g in 1:length(genera_w_singletons))	{
	genus_singletons[g] <- sum(singletons[,2]==genera_w_singletons[g])
	}
prob_genera <- sort(genus_singletons,decreasing=TRUE)
lumped <- vector(length=length(prob_genera))
names(lumped) <- names(prob_genera)
poss_cases <- sum(prob_genera>1)
for (g in 1:poss_cases)	{
	prob_child <- names(prob_genera)[g]
	gen_no <- match(prob_child,genera)
	relevant_finds <- subset(occurrences,occurrences$GTU==gen_no)
	singles <- subset(singletons,singletons[,2]==gen_no)
	lumped[g] <- 0
	
	for (s in 1:prob_genera[g])	{
		while (singles[s,4]==1 && s<prob_genera[g])	s <- s+1
		loner <- singles[s,1]
		oc_no <- match(loner,relevant_finds$OTU_Sr)
		loc <- relevant_finds$Locality[oc_no]
		loner_congen <- subset(relevant_finds,relevant_finds$Locality==loc)
		if (length(loner_congen$OTU_Sr)>1)	{
			# put on # finds to identify other singletons
			test <- cbind(loner_congen$OTU_Sr,species_summary[match(loner_congen$OTU_Sr,species_summary[,1]),3])
			others <- subset(test,test[,2]==1)
			if (length(others[,1])>1)	{
				# "synonymize" species
				change <- match(others[,1],occurrences$OTU_Sr)
				occurrences$OTU_Sr[change] <- loner
				# check other singletons off as done
				redone <- match(others[,1],singles[,1])
				singles[redone,4] <- 1
				# add to sum of lumped
				lumped[g] <- lumped[g]+(length(redone)-1)
				}
			}
		}
	}

# need to clear occurrences here.  Delete duplicate finds for senior species
cleaned_occurrences <- evanesco_duplicate_occurrences(occurrences,species_info)
return(cleaned_occurrences[,1:3])
}

put_rock_units_into_regime_partitions <- function(code)	{
rocks <- length(code)
regime <- vector(length=rocks)
for (r in 1:rocks)	{
	if (word(code[r],-1)=="Toquima-Tablehead")	{
		regime[r] <- word(code[r],-1)
		}	else if (word(code[r],1)=="Old-World")	{
		regime[r] <- "Gondwana"
		}	else {
		regime[r] <- word(code[r],1)
		}
	}
return(regime)
}

partition_localities_by_sampling_regimes <- function(locality_data,rock_unit,regime)	{
rocks <- length(rock_unit[,1])
locality_regimes <- vector(length=length(locality_data$Locality))
for (r in 1:rocks)	{
	# skip rocks with no collections
	while (sum(locality_data$Formation_Sr==as.numeric(rock_unit[r,1]))==0 && r<rocks)	r <- r+1
	form_coll <- subset(locality_data,locality_data$Formation_Sr==as.numeric(rock_unit[r,1]))
	# those locality_data localities that are in the localities assigned to formation scored
	locality_regimes[locality_data$Locality %in% form_coll$Locality] <- regime[r]
#	print(c(r,length(locality_regimes)))
	}
return(locality_regimes)
}

rock_units_per_stage_per_regime <- function(locality_data,rock_unit)	{
code <- paste(formation_info$Province,formation_info$Subprovince,sep=" ")
regime <- put_rock_units_into_regime_partitions(code)
locality_regimes <- partition_localities_by_sampling_regimes(locality_data,rock_unit,regime)
resolved_localities <- cbind(locality_data,locality_regimes)
resolved_localities <- subset(resolved_localities,locality_data$inexact==0)
mnbin <- min(as.numeric(resolved_localities$bin_lb))
mxbin <- max(as.numeric(resolved_localities$bin_lb))
realms <- sort(unique(regime))
partitions <- length(realms)
res_coll <- dim(resolved_localities)[1]
rock_units_per_stage_per_realm <- matrix(0,mxbin,partitions)
colnames(rock_units_per_stage_per_realm) <- realms
rownames(rock_units_per_stage_per_realm) <- rock_unit[match(sort(unique(as.numeric(rock_unit[,8]))),as.numeric(rock_unit[,8])),6]

#dim(subset(resolved_localities,resolved_localities$locality_regimes=="Laurentia"))
for (b in mnbin:mxbin)	{
	bin_localities <- subset(resolved_localities,resolved_localities$bin_lb==b)
	bin_regimes <- sort(unique(bin_localities$locality_regimes))
	for (p in 1:length(bin_regimes))	{
		xx <- subset(bin_localities,bin_localities$locality_regimes==as.character(bin_regimes[p]))
		pp <- match(bin_regimes[p],realms)
		rock_units_per_stage_per_realm[b,pp] <- length(unique(xx$Formation_Sr))
		}
	}
return(rock_units_per_stage_per_realm)
}

localities_per_stage_per_regime <- function(locality_data,rock_unit)	{
code <- paste(formation_info$Province,formation_info$Subprovince,sep=" ")
regime <- put_rock_units_into_regime_partitions(code)
locality_regimes <- partition_localities_by_sampling_regimes(locality_data,rock_unit,regime)
resolved_localities <- cbind(locality_data,locality_regimes)
resolved_localities <- subset(resolved_localities,locality_data$inexact==0)
mnbin <- min(as.numeric(resolved_localities$bin_lb))
mxbin <- max(as.numeric(resolved_localities$bin_lb))
realms <- sort(unique(regime))
partitions <- length(realms)
res_coll <- dim(resolved_localities)[1]
collections_per_stage_per_realm <- matrix(0,mxbin,partitions)
colnames(collections_per_stage_per_realm) <- realms
rownames(collections_per_stage_per_realm) <- rock_unit[match(sort(unique(as.numeric(rock_unit[,8]))),as.numeric(rock_unit[,8])),6]

dim(subset(resolved_localities,resolved_localities$locality_regimes=="Laurentia"))
for (b in mnbin:mxbin)	{
	bin_localities <- subset(resolved_localities,resolved_localities$bin_lb==b)
	bin_regimes <- sort(unique(bin_localities$locality_regimes))
	for (p in 1:length(bin_regimes))	{
		xx <- subset(bin_localities,bin_localities$locality_regimes==as.character(bin_regimes[p]))
		pp <- match(bin_regimes[p],realms)
		collections_per_stage_per_realm[b,pp] <- length(xx$bin_lb)
		}
	}
return(collections_per_stage_per_realm)
}

revelio_controlled_localities <- function(species,control_category,control_groups,occurrences,locality_info)	{
controls <- length(control_groups)
# sum binary vectors for whether species belong to one of the control groups
# NOTE: there must be some way to do this in one command....
#sum(control_category %in% control_groups)
controlled <- (control_category==control_groups[1])
for (cg in 2:controls)	controlled <- controlled+(control_category==control_groups[cg])
controlled_find <- vector(length=length(occurrences$OTU_Sr))
for (sp in 1:length(species))
	controlled_find[occurrences$OTU_Sr %in% species[sp]] <- controlled[sp]
# create vector stating whether occurrence fits controls
# match(occurrences$OTU_Sr,species) converts initial species number to ranked species number
keeper <- controlled[match(occurrences$OTU_Sr,species)]
# list localities and whether species is in control groups
link <- cbind(occurrences$Locality,keeper)
# reduce to just each collection with 1+ control taxa
useful <- unique(subset(link,link[,2]==1))
# separate out the remainder
useless <- cbind(locality_info$Locality[!locality_info$Locality %in% useful],integer(length(length(locality_info$Locality)-length(useful[,1]))))
all <- rbind(unique(useful),unique(useless))
dd <- order(all[,1])
# return unique localities
return(unique(all[dd,])[,2])
}

revelio_controlled_rock_units <- function(species,control_category,control_groups,occurrences,locality_info)	{
controls <- length(control_groups)
# sum binary vectors for whether species belong to one of the control groups
# NOTE: there must be some way to do this in one command....
controlled <- (control_category==control_groups[1])
for (cg in 2:controls)	controlled <- controlled+(control_category==control_groups[cg])
controlled_find <- vector(length=length(occurrences$OTU_Sr))
for (sp in 1:length(species))		controlled_find[occurrences$OTU_Sr %in% species[sp]] <- controlled[sp]
# create vector stating whether occurrence fits controls
# match(occurrences$OTU_Sr,species) converts initial species number to ranked species number
keeper <- controlled[match(occurrences$OTU_Sr,species)]
# list localities and whether species is in control groups
link <- cbind(occurrences$Locality,keeper)
# reduce to just each collection with 1+ control taxa
useful <- unique(subset(link,link[,2]==1))
# separate out the remainder
useless <- cbind(locality_info$Locality[!locality_info$Locality %in% useful],integer(length(length(locality_info$Locality)-length(useful[,1]))))
all <- rbind(unique(useful),unique(useless))
dd <- order(all[,1])
# return unique localities
return(unique(all[dd,])[,2])
}

count_rock_units_per_bin_simple <- function(rock_unit,bin)	{
return(sum((as.numeric(set_rocks[,8])==bin)*(as.numeric(set_rocks[,9])==bin)))
}

count_rock_units_per_bin <- function(bin,locality_info)	{
bin_locales <- subset(locality_info,locality_info$bin_lb==bin)
bin_locales <- subset(bin_locales,bin_locales$bin_ub==bin)
return(length(unique(bin_locales$Formation_Sr)))
}

count_collections_per_bin <- function(bin,locality_info)	{
return(sum((locality_info$bin_lb==bin)*(locality_info$bin_ub==bin)))
}

count_occurrences_per_bin <- function(bin,occurrences,locality_info)	{
bin_locales <- subset(locality_info,locality_info$bin_lb==bin)
bin_locales <- subset(bin_locales,bin_locales$bin_ub==bin)
return(length(occurrences$Locality[occurrences$Locality %in% bin_locales$Locality]))
}

count_occurrences_per_taxon <- function(taxon,finds)	{
# finds should just be a list of (say) species numbers that go to locality numbers
# 	the localtiy numbers themselves are not needed here
return(sum(finds==taxon))
}

count_localities_per_taxon <- function(taxon,occr_data,loc_col,taxon_col)	{
# finds should just be a list of (say) species numbers that go to locality numbers
#	the localtiy numbers themselves are not needed here
rel_data <- subset(occr_data,occr_data[,taxon_col]==taxon)
return(length(unique(rel_data[,loc_col])))
}

count_subtaxa_per_taxon <- function(taxon,subtaxa)	{
# subtaxa should just be a list of (say) genus numbers that go to species numbers
#	the species numbers themselves are not needed here
return(sum(subtaxa==taxon))
}

#reducto_lists(occurrence_data=bin_finds_spc,criterion=daughters[1],column=match("SPC",colnames(bin_finds_spc)))
reducto_lists <- function(occurrence_data,criterion,column)	{
print(paste("Doing ",criterion,sep=""))
return(subset(occurrence_data,occurrence_data[,column]==criterion))
}

accio_rock_units_for_a_taxon_per_bin <- function(taxon,occurrences,locality_info)	{
taxon_finds <- subset(occurrences,occurrences$OTU_Sr==taxon)
taxon_localities <- locality_info[(locality_info$Locality %in% taxon_finds$Locality),]
bin_nos <- (1:max(locality_info$bin_ub))
return(simplify2array(lapply(bin_nos,count_rock_units_per_bin,locality_info=taxon_localities)))
}

accio_collections_per_taxon_per_bin <- function(taxon,occurrences,locality_info)	{
taxon_finds <- subset(occurrences,occurrences$OTU_Sr==taxon)
taxon_localities <- locality_info[(locality_info$Locality %in% taxon_finds$Locality),]
bin_nos <- (1:max(locality_info$bin_ub))
return(simplify2array(lapply(bin_nos,count_rock_units_per_bin,locality_info=taxon_localities)))
}

accio_taxon_ranges_from_occurrence_data <- function(taxon,occurrences,locality_info)	{
taxon_finds <- subset(occurrences,occurrences$OTU_Sr==taxon)
taxon_localities <- locality_info[(locality_info$Locality %in% taxon_finds$Locality),]
if (dim(taxon_localities)[1]==1)	{
	output_range <- c(1,0.1)
	names(output_range) <- c("Bins","Ma")
	return(c(1,0.1))
	} else	{
	rng_bin <- 1+(max(taxon_localities$bin_ub)-min(taxon_localities$bin_lb))
	rng_myr <- max(abs(taxon_localities$ma_lb))-min(abs(taxon_localities$ma_ub))
	output_range <- c(rng_bin,rng_myr)
	names(output_range) <- c("Bins","Ma")
	return(c(rng_bin,rng_myr))
	}
}

## get first and last appearance data
accio_taxon_fa_and_la_data <- function(taxon,occurrences,locality_info)	{
if (sum(occurrences$OTU_Sr==taxon)>0)	{
	taxon_finds <- subset(occurrences,occurrences$OTU_Sr==taxon)
	taxon_localities <- locality_info[(locality_info$Locality %in% taxon_finds$Locality),]
	fa_bin <- min(taxon_localities$bin_lb)
	la_bin <- max(taxon_localities$bin_ub)
	fa_myr <- max(abs(taxon_localities$ma_lb))
	la_myr <- min(abs(taxon_localities$ma_ub))
	output_fa_la <- c(fa_bin,la_bin,fa_myr,la_myr)
	}	else {
	output_fa_la <- c(0,0,0,0)
	}
names(output_fa_la) <- c("FA_Bin","LA_Bin","FA_Ma","LA_Ma")
return(output_fa_la)
}

### ROUTINE TO GET FINDS FOR ONE TAXON FROM STANDARD OCCURRENCE FILE 2017-05-04
accio_occurrences_per_bin_per_taxon_from_occurence_file <- function(file_name,header=TRUE,bin_col=1,loc_col=2,taxon_col,lump=TRUE)	{
# file_name: occurrence data from a flat file with columns giving bin #, locality # & taxon #
# bin_col: column giving the bin number
# loc_col: column giving the locality number
# taxon_col: column giving the taxon number
find_data <- read.table(file=file_name, header, stringsAsFactors=FALSE)
taxa <- sort(unique(find_data[,taxon_col]))
#ntaxa <- max(find_data[,taxon_col])	# get maximum taxon number
ntaxa <- length(taxa)	# get maximum taxon number
mxbin <- max(find_data[,bin_col])
mnbin <- min(find_data[,bin_col])
taxon_bin_find_matrix <- matrix(0,mxbin,ntaxa)
tbfm <- sapply((mnbin:mxbin),accio_all_taxon_occurrences_per_bin,find_data,taxa,bin_col,loc_col,taxon_col,lump)
colnames(tbfm) <- (mnbin:mxbin)
rownames(tbfm) <- taxa
taxon_bin_find_matrix <- base::t(tbfm)
#for (b in mnbin:bins)	{
#	bin_data <- subset(find_data,find_data[,bin_col]==b)
#	bin_taxa <- sort(unique(bin_data[,taxon_col]))
#	bin_finds <- sapply(bin_taxa,accio_individual_taxon_occurrences_per_bin,bin_data,loc_col,taxon_col,lump)
#	taxon_bin_find_matrix[b,(1:ntaxa)[taxa %in% bin_taxa]] <- bin_finds
#	}
return(taxon_bin_find_matrix)
}

### ROUTINE TO GET FINDS FROM STANDARD OCCURRENCE FILE MATRIX 2017-05-04
accio_all_taxon_occurrences_per_bin <- function(b,find_data,taxa,bin_col=1,loc_col=2,taxon_col,lump=TRUE)	{
ntaxa <- length(taxa)
bin_data <- subset(find_data,find_data[,bin_col]==b)
bin_taxa <- sort(unique(bin_data[,taxon_col]))
bin_finds <- sapply(bin_taxa,accio_individual_taxon_occurrences_per_bin,bin_data,loc_col,taxon_col,lump)
bf <- rep(0,ntaxa)
bf[(1:ntaxa)[taxa %in% bin_taxa]] <- bin_finds
return(bf)
}

### ROUTINE TO GET FINDS FOR ONE TAXON FROM STANDARD OCCURRENCE FILE MATRIX 2017-05-04
accio_individual_taxon_occurrences_per_bin <- function(taxon,bin_data,loc_col=2,taxon_col,lump=TRUE)	{
taxon_data <- subset(bin_data,bin_data[,taxon_col]==taxon)
if(lump)	{
	return(length(unique(taxon_data[,loc_col])))
	} else	{
	return(nrow(taxon_data))
	}
}

### ROUTINE TO GET BIN RANGES FROM FIND MATRIX 2017-05-04
accio_strat_ranges_from_bin_by_taxa_matrix <- function(find_matrix,bins,taxa)	{
#find_matrix: bin x taxon matrix giving # occurrences per taxon per bin
#bins: array giving bin numbers corresponding to rows
#taxa: array giving taxon numbers corresponding to columns
ntaxa <- length(taxa)
nbins <- length(bins)
tx <- taxa
#ranges <- sapply(tx=taxa,extract_taxon_range_from_bin_by_taxa_matrix,find_matrix,bins)
ranges <- c()
for (tx in 1:ntaxa)	{
	ranges <- rbind(ranges,extract_taxon_range_from_bin_by_taxa_matrix(tx,find_matrix,bins))
#	taxon <- taxa[tx]
#	taxon_finds <- find_matrix[,tx]
	}
rownames(ranges) <- taxa
colnames(ranges) <- c("FA","LA")
return(ranges)
}

### ROUTINE TO GET BIN RANGES FOR A TAXON FROM A VECTOR OF 0's for absence, >0 for present
extract_taxon_range_from_bin_by_taxa_matrix	<- function(tx,find_matrix,bins)	{
return(c(min(bins[find_matrix[,tx]>0]),max(bins[find_matrix[,tx]>0])))
}

# ROUTINE TO COUNT GAPS WITHIN BIN RANGES FOR A TAXON 2015-05-04
count_gaps_in_ranges_from_bin_by_taxa_matrix <- function(tx,find_matrix,bins)	{
return(1+(max((1:length(bins))[find_matrix[,tx]>0])-min((1:length(bins))[find_matrix[,tx]>0]))-sum(find_matrix[,tx]>0))
}

# 2017-05-04: Routine to break up taxa with overly long gaps into distinct taxa
breakup_gaps_within_ranges <- function(find_matrix,bins,taxa,max_gap=3)	{
ntaxa <- length(taxa)
nbins <- length(bins)
find_matrix_alt <- find_matrix
added_taxa <- c()
added_finds <- c()
for (tx in 1:ntaxa)	{
	gp <- count_gaps_in_ranges_from_bin_by_taxa_matrix(tx,find_matrix,bins)
	if (gp >= max_gap)	{
		tx_bins <- bins[find_matrix[,tx]>0]
		ind_gaps <- tx_bins[2:length(tx_bins)]-(1+tx_bins[1:(length(tx_bins)-1)])
		poly_gaps <- sum(ind_gaps>max_gap)
		if (poly_gaps>0)	{
			dummy_finds <- c()
			for (nb in 1:(poly_gaps+1))	dummy_finds <- cbind(dummy_finds,find_matrix[,tx])
			gg <- fst <- 1
			for (xg in 1:poly_gaps)	{
				while (ind_gaps[gg]<max_gap)	gg <- gg+1
				lst <- match(tx_bins[gg],bins)
				# clear "early" part
				for (bb in (lst+1):nbins)	dummy_finds[bb,xg] <- 0
				for (bb in fst:lst)	for (xxg in xg:poly_gaps)	dummy_finds[bb,xxg+1] <- 0
				fst <- lst+1
				gg <- gg+1
				}
			colnames(dummy_finds) <- rep(taxa[tx],ncol(dummy_finds))
			find_matrix_alt[,tx] <- dummy_finds[,1]
#			find_matrix_alt <- cbind(find_matrix_alt,dummy_finds[,2:(poly_gaps+1)])
			added_finds <- cbind(added_finds,dummy_finds[,2:(poly_gaps+1)])
#			added_taxa <- c(added_taxa,rep(taxa[tx],poly_gaps))
			for (pg in 1:poly_gaps)	added_taxa <- c(added_taxa,taxa[tx]+pg/(poly_gaps+1))
			}
		}
	}
#dim(find_matrix)
#dim(find_matrix_alt)
a <- 1+ncol(find_matrix)
b <- ncol(find_matrix)+ncol(added_finds)
colnames(added_finds) <- added_taxa
find_matrix_alt <- cbind(find_matrix_alt,added_finds)
#for (c in 1:length(added_taxa))	names(find_matrix_alt[,a+c]) <- added_taxa[c]
#colnames(find_matrix_alt[,a:b]) <- added_taxa
return(find_matrix_alt)
}

# 2017-05-04: count unique localities from occurrence file
accio_localities_per_bin_per_taxon_from_occurence_file <- function(file_name,header=TRUE,bin_col=1,loc_col=2)	{
find_data <- read.table(file=file_name, header, stringsAsFactors=FALSE)
bins <- sort(unique(find_data[,bin_col]))
nbins <- length(bins)
bcoll <- c()
for (b in 1:nbins)	{
	bcoll <- c(bcoll,length(unique(subset(find_data,find_data[,bin_col]==bins[b])[,loc_col])))
	}
return(bcoll)
}

# get samples per taxon
accio_samples_per_taxon <- function(taxa,sample_units)	{
unique_taxon <- sort(unique(taxa))
samples <- length(sample_units)
taxon_samples <- c()
for (u in 1:length(unique_taxon))	{
	tfinds <- (1:samples)[taxa==unique_taxon[u]]
	taxon_samples <- c(taxon_samples,length(unique(sample_units[tfinds])))
	}
return(cbind(unique_taxon,taxon_samples))
}

summed_common_species_occurrence_rate <- function(occupancy,sites)	{
return(sum(-1*log(1-(occupancy/sites))))
}

scor_lambda <- function(pij)	{
return(-1*log(1 - pij))
}

	#### JACK LIVES ####
sepkoskify_one_taxon_from_occurrence_data <- function(notu,taxon_no,interval_no,sample_no,sample_age_lb,sample_age_ub)  {
fossil_record <- data.frame(taxon_no,interval_no,sample_no,sample_age_lb,sample_age_ub);
taxon_record <- unique(subset(fossil_record,fossil_record$taxon_no==notu));
if (nrow(taxon_record)>0)	{
	ma_max <- max(taxon_record$sample_age_lb);
	ma_min <- min(taxon_record$sample_age_ub);
	bin_lb <- min(taxon_record$stage_no);
	bin_ub <- max(taxon_record$stage_no);
	} else	{
	# fill vectors with nothing if no finds.
	ma_max <- ma_min <- bin_lb <- bin_ub <- 0;
	}
return(c(bin_lb,bin_ub,ma_max,ma_min));
}

#taxon="Onniella bancrofti"
sepkoskify_paleodb_data_one_taxon <- function(taxon,pbdb_finds,transpose=FALSE)	{
taxon_record <- subset(pbdb_finds,pbdb_finds$accepted_name==taxon);
if (nrow(taxon_record)==0)	{
	# fill vectors with nothing if no finds.
	taxon_found <- which(pbdb_finds==taxon,arr.ind = T);
	if (length(taxon_found)>0)	{
		# taxon found under another name!
		taxon_columns <- match(c("phylum","class","order","family","genus","subgenus"),colnames(pbdb_finds));
		if (sum(taxon_columns %in% taxon_found[1,2])>0)	{
			coll_w_taxon <- which(pbdb_finds==taxon,arr.ind = T)[,1];
			taxon_record <- pbdb_finds[coll_w_taxon,];
			} else if (sum(unique(colnames(paleodb_finds)[taxon_found[,2]]) %in% "identified_name")>0)	{
			coll_w_taxon <- which(pbdb_finds==taxon,arr.ind = T)[,1];
			taxon_record <- pbdb_finds[coll_w_taxon,];
			}
		}
	}

if (nrow(taxon_record)==0)	{
	# fill vectors with nothing if no finds.
	ma_max <- ma_min <- bin_lb <- bin_ub <- 0;
	} else {
	taxon_record_set <- subset(taxon_record,taxon_record$bin_lb==taxon_record$bin_ub);
	taxon_record_fuzzy <- subset(taxon_record,taxon_record$bin_lb!=taxon_record$bin_ub);
	if (nrow(taxon_record_set)>0)	{
		ma_max <- max(taxon_record_set$ma_lb);
		ma_min <- min(taxon_record_set$ma_ub);
		bin_lb <- min(taxon_record_set$bin_lb);
		bin_ub <- max(taxon_record_set$bin_ub);
		if (nrow(taxon_record_fuzzy)>0)	{
			if (min(taxon_record_fuzzy$bin_ub) < bin_lb)	{
				bin_lb <- min(taxon_record_fuzzy$bin_ub);
				if (!is.na(match(bin_lb,taxon_record$bin_lb)))
					taxon_record$ma_lb[match(bin_lb,taxon_record$bin_lb)];
				} else if (max(taxon_record_fuzzy$bin_lb) > bin_ub)	{
				bin_ub <- max(taxon_record_fuzzy$bin_lb);
				if (!is.na(match(bin_ub,taxon_record$bin_ub)))
					taxon_record$ma_ub[match(bin_ub,taxon_record$bin_ub)];
				}
			}
		} else if (nrow(taxon_record_fuzzy)>0)	{
		ma_max <- max(taxon_record$ma_lb[taxon_record$ma_lb>=max(taxon_record$ma_ub)]);
		ma_min <- min(taxon_record$ma_ub[taxon_record$ma_ub<=min(taxon_record$ma_lb)]);
		bin_lb <- min(taxon_record$bin_lb[taxon_record$bin_lb<=min(taxon_record$bin_ub)]);
		bin_ub <- max(taxon_record$bin_ub[taxon_record$bin_ub>=max(taxon_record$bin_lb)]);
		}
	}
if (transpose)	{
	output <- base::t(data.frame(taxon=as.character(taxon),bin_lb=as.numeric(bin_lb),bin_ub=as.numeric(bin_ub),ma_max=as.numeric(ma_max),ma_min=as.numeric(ma_min)));
	colnames(output) <- taxon;
	} else	{
	output <- data.frame(taxon=as.character(taxon),bin_lb=as.numeric(bin_lb),bin_ub=as.numeric(bin_ub),ma_max=as.numeric(ma_max),ma_min=as.numeric(ma_min));
	rownames(output) <- taxon;
	}
return(output);
}
	
sepkoskify_paleodb_data <- function(pbdb_finds,taxon_names,interval_names="")  {
#taxon <- taxon_names;
tranpose <- FALSE;
compendium <- c();
#compend <- sapply(taxon,sepkoskify_paleodb_data_one_taxon,pbdb_finds,tranpose);
for (i in 1:length(taxon_names))
	compendium <- rbind(compendium,sepkoskify_paleodb_data_one_taxon(taxon=taxon_names[i],pbdb_finds))
#	print(c(i,nrow(compendium)));
#compend <- data.frame(base::t(sapply(notu,sepkoskify_one_taxon_from_occurrence_data,taxon_no,stage_no,sample_no,sample_age_lb,sample_age_ub)));
return(compendium);
}

sepkoskify_occurrence_data <- function(taxon_no,stage_no,sample_no,sample_age_lb,sample_age_ub,sampled_taxa,interval_names="",condense=TRUE)  {
# taxon_no: vector giving taxon number from PaleoDB or elsewhere for occurrences
# stage_no: vector given the rank of the stratigraphic unit, with 1 the oldest
# sample_no: the collection number from PaleoDB or elswhere; can also be rock unit number
# sample_age_lb: lower bound on sample age in millions of years, with 400 older than 300
# sample_age_ub: upper bound on sample age in millions of years, with 400 older than 300
# interval_names: names of the stratigraphic units, with interval_names[1] being the name going to stage_no[n]=1
# condense: if TRUE, then only information for taxa actually present is returned;
#	if FALSE, then absent numbers (e.g., 2 if 1 and 3 present) are inserted with no information
#	FALSE is useful if you are building the compendium over increasingly refined stratigraphic data
fossil_record <- data.frame(taxon_no=as.numeric(taxon_no),stage_no=as.numeric(stage_no),sample_no=as.numeric(sample_no),sample_age_lb=as.numeric(sample_age_lb),sample_age_ub=as.numeric(sample_age_ub));
taxonomic_dictionary <- unique(data.frame(taxon_no=as.numeric(taxon_no),sampled_taxa=as.character(sampled_taxa)));
taxonomic_dictionary <- taxonomic_dictionary[order(taxonomic_dictionary$taxon_no),]
ntaxa <- length(taxonomic_dictionary$sampled_taxa);
notu <- taxonomic_dictionary$taxon_no;
compend <- base::t(sapply(notu,sepkoskify_one_taxon_from_occurrence_data,taxon_no,stage_no,sample_no,sample_age_lb,sample_age_ub));
colnames(compend) <- c("bin_lb","bin_ub","ma_max","ma_min");
if (length(interval_names)>0)	{
	compendium <- data.frame(bin_lb=as.numeric(frak_it(compend[,1])),
							 bin_ub=as.numeric(frak_it(compend[,2])),
							 ma_max=as.numeric(frak_it(compend[,3])),
							 ma_min=as.numeric(frak_it(compend[,4])),
							 stage_lb=as.character(as.character(interval_names[frak_it(compend[,1])])),
							 stage_ub=as.character(as.character(interval_names[frak_it(compend[,2])])),
							 stringsAsFactors = FALSE);
	} else	{
	compendium <- data.frame(bin_lb=as.numeric(frak_it(compend[,1])),
							 bin_ub=as.numeric(frak_it(compend[,2])),
							 ma_max=as.numeric(frak_it(compend[,3])),
							 ma_min=as.numeric(frak_it(compend[,4])));
	}
rownames(compendium) <- taxonomic_dictionary$sampled_taxa;
if (!condense)	{
	absent_taxa <- (1:max(notu))[!(1:max(notu)) %in% notu];
	ataxa <- length(absent_taxa);
	all_taxa <- c(notu,absent_taxa);
	if (length(interval_names)>0)	{
		dummy <- data.frame(
		bin_lb=as.numeric(rep(0,ataxa)),
		bin_ub=as.numeric(rep(0,ataxa)),
		ma_max=as.numeric(rep(0,ataxa)),
		ma_min=as.numeric(rep(0,ataxa)),
		stage_lb=as.character(rep("",ataxa)),
		stage_ub=as.character(rep("",ataxa)),
		stringsAsFactors = FALSE);
		} else	{
		dummy <- data.frame(
		bin_lb=as.numeric(rep(0,ataxa)),
		bin_ub=as.numeric(rep(0,ataxa)),
		ma_max=as.numeric(rep(0,ataxa)),
		ma_min=as.numeric(rep(0,ataxa)));
		}
	rownames(dummy) <- absent_taxa;
	compendium <- rbind(compendium,dummy);
	compendium <- compendium[order(all_taxa),];
	}
return(compendium);
}

# score how many sum of range extensions required to make
procrustean_binning <- function(range,poss_range,debug=FALSE)	{
# range: taxon range in bin numbers (with 1 older than 2)
# poss_range: lower and upper bounds on collection age in bin numbers
if (debug)	print(range);
obs_bins <- range[1]:range[2];
poss_bins <- poss_range[1]:poss_range[2];
#ok_bins <- poss_bins[poss_bins %in% obs_bins];
#if (length(ok_bins)==0)
#	ok_bins <- poss_bins;
bin_score <- c();
for (pb in 1:length(poss_bins))	{
	bin_score <- c(bin_score,-min(abs(poss_bins[pb]-obs_bins)));
	}
return(bin_score);
}

# score how many sum of range extensions required to make
procrustean_binning_one_taxon <- function(range,poss_range,debug=FALSE)	{
# range: taxon range in bin numbers (with 1 older than 2)
# poss_range: lower and upper bounds on collection age in bin numbers
if (debug)	print(range);
obs_bins <- range[1]:range[2];
poss_bins <- poss_range[1]:poss_range[2];
#ok_bins <- poss_bins[poss_bins %in% obs_bins];
#if (length(ok_bins)==0)
#	ok_bins <- poss_bins;
bin_score <- c();
for (pb in 1:length(poss_bins))	{
	bin_score <- c(bin_score,min(abs(poss_bins[pb]-obs_bins)));
	}
return(bin_score);
}

transmogrify_accepted_species_name <- function(identified_name,accepted_genus)	{
#transmogrify_accepted_species_name <- function(name_to_fix)	{
#print(name_to_fix[1:2])
#accepted_genus <- name_to_fix[1];
species_epithet <- diffindo_species_epithets(identified_name);
return(paste(accepted_genus,species_epithet))
}

accio_minimum_strat_ranges <- function(taxa,paleodb_finds,paleodb_collections,hierarchical_chronostrat)	{
# taxa: vector of taxon names
# paleodb_finds: data.frame of Paleobiology Database finds
# paleodb_collections: data.frame of collections from the PaleoDB corresponding to finds & linked by collection_no
# hierarchical_chronostrat: chronostratigraphic scheme with information aboust subintervals
# taxon_rank: ranke (genus, species, etc.) at which analysis is conducted
# lump_cooccr: if TRUE, the only one occurrence per colleciton
# constrain: if TRUE, then do not extend possible range of find counts beyond the latest possible FA & earliest possible LA
notu <- length(taxa);
n_intervals <- max(hierarchical_chronostrat$bin_last);
if (taxon_rank=="species")
	taxon_rank <- "accepted_name";
taxon_col <- match(taxon_rank,colnames(paleodb_finds));

occurrence_in_bin <- array(0,dim=c(notu,n_intervals));
colnames(occurrence_in_bin) <- hierarchical_chronostrat$interval[hierarchical_chronostrat$bin_first==hierarchical_chronostrat$bin_last];
rownames(occurrence_in_bin) <- taxa;

for (tx in 1:notu)	{
	# get PaleoDB collection numbers
#	print(taxa[tx]);
	taxon_colls <- as.numeric(paleodb_finds$collection_no[paleodb_finds[,taxon_col]==taxa[tx]]);
	# get which collection in this data set the PaleoDB colleciton is
	coll_no <- match(taxon_colls,paleodb_collections$collection_no);
	early_interval <- as.character(refined_collections$interval_lb[coll_no]);
	late_interval <- as.character(refined_collections$interval_ub[coll_no]);
	f_p_b <- finds_per_bin[tx,] <- tally_occurrences_by_subinterval(coll_no,early_interval,late_interval,hierarchical_chronostrat,lump_cooccr,constrain);
	definite_finds_per_bin[tx,] <- tally_definite_occurrences_by_subinterval(coll_no,early_interval,late_interval,hierarchical_chronostrat);
	occurrence_in_bin[tx,] <- tally_presence_by_subinterval(coll_no,early_interval,late_interval,hierarchical_chronostrat);
	# latest possible first appearance would be the earliest late interval
	}
output <- list(finds_per_bin,definite_finds_per_bin,occurrence_in_bin);
names(output) <- c("finds_per_bin","definite_finds_per_bin","occurrence_in_bin");
return(output);
}

tally_minimum_stratigraphic_ranges <- function(ma_lb,ma_ub,interval_lb,interval_ub,hierarchical_chronostrat)	{
fa_latest <- min(unique(hierarchical_chronostrat$bin_last[match(interval_ub,hierarchical_chronostrat$interval)]));
la_earliest <- max(unique(hierarchical_chronostrat$bin_first[match(interval_lb,hierarchical_chronostrat$interval)]));

#ttl_finds <- length(coll_no);
bins_early <- hierarchical_chronostrat$bin_first[match(interval_lb,hierarchical_chronostrat$interval)];
bins_late <- hierarchical_chronostrat$bin_last[match(interval_ub,hierarchical_chronostrat$interval)];
bins_early[bins_early<fa_latest] <- fa_latest;
bins_late[bins_late>la_earliest] <- la_earliest;

definite_bins <- unique(bins_early[bins_early==bins_late]);
output <- data.frame(ma_fa_mn=0,ma_la_mx=0,interval_lb=as.character(""),interval_ub=as.character(""),stringsAsFactors = F);
if (length(definite_bins)>0)	{
	finest_chronostrat <- hierarchical_chronostrat[hierarchical_chronostrat$bin_first==hierarchical_chronostrat$bin_last,];
	if (sum(sort(ma_ub,decreasing = T) > sort(ma_lb))>0)	{
		output$ma_fa_mn <- as.numeric(max(ma_ub[ma_ub>min(ma_lb)]));
		output$ma_la_mx <- as.numeric(min(ma_lb[ma_lb<max(ma_ub)]));
#		output[1,] <- c(max(ma_ub[ma_ub>min(ma_lb)]),min(ma_lb[ma_lb<max(ma_ub)]),finest_chronostrat$interval[min(definite_bins)],finest_chronostrat$interval[max(definite_bins)]);
		} else	{
		output$ma_fa_mn <- as.numeric(min(ma_lb));
		output$ma_la_mx <- as.numeric(max(ma_ub));
#		output[1,] <- c(min(ma_lb),max(ma_ub),finest_chronostrat$interval[min(definite_bins)],finest_chronostrat$interval[max(definite_bins)])
		}
	output$interval_lb <- as.character(finest_chronostrat$interval[min(definite_bins)]);
	output$interval_ub <- as.character(finest_chronostrat$interval[max(definite_bins)]);
	}
return(output);
}

tally_fuzzy_stratigraphic_ranges <- function(ma_lb,ma_ub,interval_lb,interval_ub,hierarchical_chronostrat)	{
fa_latest <- min(unique(hierarchical_chronostrat$bin_last[match(interval_ub,hierarchical_chronostrat$interval)]));
la_earliest <- max(unique(hierarchical_chronostrat$bin_first[match(interval_lb,hierarchical_chronostrat$interval)]));

#ttl_finds <- length(coll_no);
bins_early <- hierarchical_chronostrat$bin_first[match(interval_lb,hierarchical_chronostrat$interval)];
bins_late <- hierarchical_chronostrat$bin_last[match(interval_ub,hierarchical_chronostrat$interval)];
bins_early[bins_early<fa_latest] <- fa_latest;
bins_late[bins_late>la_earliest] <- la_earliest;

definite_bins <- unique(bins_early[bins_early==bins_late]);
output <- data.frame(ma_fa_mx=0,ma_fa_mn=0,ma_la_mx=0,ma_la_mn=0,interval_lb=as.character(""),interval_ub=as.character(""),stringsAsFactors = F);
if (length(definite_bins)>0)	{
	finest_chronostrat <- hierarchical_chronostrat[hierarchical_chronostrat$bin_first==hierarchical_chronostrat$bin_last,];
	# extremes are easy
	output$ma_fa_mx <- as.numeric(max(ma_lb));						# earliest possible FA
	output$ma_la_mn <- as.numeric(min(ma_ub));						# latest possible LA
	# latest FA must precede or coincide with latest possible LA & vice versa
	output$ma_fa_mn <- as.numeric(max(ma_ub));						# latest possible FA
	output$ma_la_mx <- as.numeric(min(ma_lb));						# earliest possible LA
	output$interval_lb <- as.character(finest_chronostrat$interval[match(min(definite_bins),finest_chronostrat$bin_first)]);
	output$interval_ub <- as.character(finest_chronostrat$interval[match(max(definite_bins),finest_chronostrat$bin_first)]);
	}
return(output);
}

	#### BIOSTRATIGRAPHY ROUTINES ####
# routine to opimize the stage of uncertain collections based on fossil assemblages.
# paleodb_finds = rbind(control_occurrences,zone_occurrences)
# paleodb_collections_old = paleodb_collections
# paleodb_collections = refined_collections
# routine to opimize the stage of uncertain collections based on fossil assemblages.
optimo_paleodb_collection_and_occurrence_stratigraphy <- function(paleodb_finds,paleodb_collections,hierarchical_chronostrat,zone_database,update_search=T)	{
# rescore collections if there is any lumping of reported stages into useful stages
ncolls <- nrow(paleodb_collections);
nstages <- max(hierarchical_chronostrat$bin_last);

# delete redundant occurrences of species in localities.  (This usually reflects two co-occuring
# 	species being synonymized).
paleodb_finds <- evanesco_duplicate_occurrences_paleodb(occurrences=paleodb_finds);
	
# make sure that key collections fields are characters, not factors
paleodb_collections$early_interval <- as.character(paleodb_collections$early_interval);
paleodb_collections$late_interval <- as.character(paleodb_collections$late_interval);
paleodb_collections$interval_lb <- as.character(paleodb_collections$interval_lb);
paleodb_collections$interval_ub <- as.character(paleodb_collections$interval_ub);

# update names of un-entered species
noccr <- nrow(paleodb_finds);
to_fix <- (1:noccr)[paleodb_finds$accepted_rank %in% c("genus","subgenus")];
accepted_genus <- paleodb_finds$genus[to_fix];
identified_name <- paleodb_finds$identified_name[to_fix];
fixed_names <- c();
for (fn in 1:length(to_fix))
	fixed_names <- c(fixed_names,transmogrify_accepted_species_name(identified_name[fn],accepted_genus[fn]))
paleodb_finds$accepted_name[to_fix] <- fixed_names;

# cull out "sp.
paleodb_finds <- evanesco_indeterminate_species(paleodb_finds);
taxon_names <- sort(unique(paleodb_finds$accepted_name));

# update overall info
noccr <- nrow(paleodb_finds);
ntaxa <- length(taxon_names);
taxa <- paleodb_finds$accepted_name;

# number taxa
taxon_no <- match(taxa,taxon_names);
paleodb_finds <- cbind(paleodb_finds,taxon_no);

# add collection data to occurrence data
coll_to_find_key <- match(paleodb_finds$collection_no,paleodb_collections$collection_no);
paleodb_finds$ma_lb <- paleodb_collections$ma_lb[coll_to_find_key];
paleodb_finds$ma_ub <- paleodb_collections$ma_ub[coll_to_find_key];
paleodb_finds$interval_lb <- as.character(paleodb_collections$interval_lb[coll_to_find_key]);
paleodb_finds$interval_ub <- as.character(paleodb_collections$interval_ub[coll_to_find_key]);
#paleodb_finds <- cbind(paleodb_finds,ma_lb,ma_ub,interval_lb,interval_ub);
#paleodb_finds[105,]
#paleodb_finds[109,]
finest_chronostrat <- hierarchical_chronostrat[hierarchical_chronostrat$bin_first==hierarchical_chronostrat$bin_last,];
finest_chronostrat <- finest_chronostrat[match(finest_chronostrat$bin_first,finest_chronostrat$bin_first),];

for (tx in 1:ntaxa)	{
	taxon_finds <- subset(paleodb_finds,paleodb_finds$accepted_name==taxon_names[tx]);
	if (sum(taxon_finds$flags %in% c("uncertain species","uncertain genus, uncertain species")) < nrow(taxon_finds))
		taxon_finds <- taxon_finds[!taxon_finds$flags %in% c("uncertain species","uncertain genus, uncertain species"),];
	if (tx==1)	{
		minimum_range_data <- tally_fuzzy_stratigraphic_ranges(ma_lb=as.numeric(taxon_finds$ma_lb),ma_ub=as.numeric(taxon_finds$ma_ub),interval_lb=as.character(taxon_finds$interval_lb),interval_ub=as.character(taxon_finds$interval_ub),hierarchical_chronostrat);
		} else {
		minimum_range_data <- rbind(minimum_range_data,tally_fuzzy_stratigraphic_ranges(ma_lb=as.numeric(taxon_finds$ma_lb),ma_ub=as.numeric(taxon_finds$ma_ub),interval_lb=as.character(taxon_finds$interval_lb),interval_ub=as.character(taxon_finds$interval_ub),hierarchical_chronostrat));
		}
	}
rownames(minimum_range_data) <- taxon_names;
#which(is.na(minimum_range_data),arr.ind = T)
# redo things for zone taxa: their latest FA cannot be after the zone starts and earliest LA cannot be before zone ends.
zone_taxa <- taxon_names[taxon_names %in% zone_database$zone];
for (zt in 1:length(zone_taxa))	{
	zn <- match(zone_taxa[zt],taxon_names);
#	if (minimum_range_data$ma_fa_mn[zn]==0)
#		print(zt)
	this_zone <- subset(zone_database,zone_database$zone==zone_taxa[zt]);
	if (minimum_range_data$ma_fa_mn[zn]==0)	{
		minimum_range_data$ma_fa_mx[zn] <- minimum_range_data$ma_fa_mn[zn] <- max(this_zone$ma_lb);
		minimum_range_data$ma_la_mx[zn] <- minimum_range_data$ma_la_mn[zn] <- min(this_zone$ma_ub);
		minimum_range_data$interval_lb[zn] <- this_zone$interval_lb[match(max(this_zone$ma_lb),this_zone$ma_lb)];
		minimum_range_data$interval_ub[zn] <- this_zone$interval_ub[match(min(this_zone$ma_ub),this_zone$ma_ub)];
		} else	{
		if (minimum_range_data$ma_fa_mx[zn]<max(this_zone$ma_lb))	{
			minimum_range_data$ma_fa_mx[zn] <- max(this_zone$ma_lb);
			minimum_range_data$interval_lb[zn] <- finest_chronostrat$interval[max(1,sum(max(this_zone$ma_lb)<=finest_chronostrat$ma_lb))];
			}
		if (minimum_range_data$ma_fa_mn[zn]<max(this_zone$ma_lb)){
			minimum_range_data$ma_fa_mn[zn] <- max(this_zone$ma_lb);
			}
		if (minimum_range_data$ma_la_mx[zn]>min(this_zone$ma_ub)){
			minimum_range_data$ma_la_mx[zn] <- min(this_zone$ma_ub);
			minimum_range_data$interval_ub[zn] <- finest_chronostrat$interval[max(1,sum(max(this_zone$ma_ub)<=finest_chronostrat$ma_ub))];
			}
		if (minimum_range_data$ma_la_mn[zn]>min(this_zone$ma_ub))	{
			minimum_range_data$ma_la_mn[zn] <- min(this_zone$ma_ub);
			}
		}
	}
#ttt <- (1:ntaxa)[minimum_range_data$ma_fa_mx<minimum_range_data$ma_la_mx];
#which(is.na(minimum_range_data),arr.ind = T)
bin_lb <- hierarchical_chronostrat$bin_first[match(paleodb_collections$interval_lb,hierarchical_chronostrat$interval)];
bin_ub <- hierarchical_chronostrat$bin_last[match(paleodb_collections$interval_ub,hierarchical_chronostrat$interval)];
# kluge: eliminate the need for this;
xx <- (1:ncolls)[bin_lb>bin_ub];
if (length(xx)>0)	{
	dummy <- bin_lb[xx];
	bin_lb[xx] <- bin_ub[xx];
	bin_ub[xx] <- dummy;
	paleodb_collections$interval_lb[xx] <- finest_chronostrat$interval[bin_lb[xx]];
	paleodb_collections$interval_ub[xx] <- finest_chronostrat$interval[bin_ub[xx]];
	}
#paleodb_collections$early_interval[xx]
#paleodb_collections$late_interval[xx]
#paleodb_collections$interval_lb <- finest_chronostrat$interval[bin_lb];
#paleodb_collections$interval_ub <- finest_chronostrat$interval[bin_ub];

#unfixed_collections <- paleodb_collections$collection_no[as.character(paleodb_collections$interval_lb)!=as.character(paleodb_collections$interval_ub)];
bin_fuzz <- bin_ub-bin_lb;
unfixed_collections <- paleodb_collections$collection_no[bin_fuzz>0];

#paste(unfixed_collections,collapse = ",")
orig_unfixed <- unfixed <- length(unfixed_collections);
improved <- fixed <- ncolls - unfixed;
min_gap_to_set <- 0;
min_taxa_to_set <- 3;
attempt <- 1;
reboots <- 0;
# paste(unfixed_collections,collapse=",");
fuzz <- 0:max(bin_ub-bin_lb);
#hist(bin_ub-bin_lb)
fcolls <- hist(bin_fuzz,breaks=c(min(fuzz)-1,fuzz),plot=F)$counts;
progress <- cbind(fuzz,fcolls);
if (update_search)
	print(progress);
while (improved > 0)	{
	newly_fixed <- new_and_improved <- 0;
	uc <- 0;	# uc <- match(problems[12],unfixed_collections)
 	while (uc < unfixed)	{
		uc <- uc + 1;
		coll_no <- match(unfixed_collections[uc],paleodb_collections$collection_no);
		index_species <- paleodb_finds$accepted_name[paleodb_finds$collection_no==unfixed_collections[uc]];
		index_ranges_all <- minimum_range_data[match(index_species,rownames(minimum_range_data)),];
		index_ranges <- subset(index_ranges_all,index_ranges_all$ma_fa_mn!=0);
		if (nrow(index_ranges)>0)	{
			ranges <- cbind(hierarchical_chronostrat$bin_first[match(index_ranges$interval_lb,hierarchical_chronostrat$interval)],hierarchical_chronostrat$bin_last[match(index_ranges$interval_ub,hierarchical_chronostrat$interval)]);
			# something is resetting bin_lb & bin_ub to NA; might be ma_lb &/or ma_ub doing it
			poss_range <- bin_lb[coll_no]:bin_ub[coll_no];
			
			bin_gaps <- array(0,dim=c(length(poss_range)));
			for (rr in 1:nrow(ranges))
				bin_gaps <- bin_gaps+procrustean_binning_one_taxon(range=ranges[rr,],poss_range=c(min(poss_range),max(poss_range)),debug=F);

			# case where some of the previously assigned range induces no gaps in known species
			if ((min(bin_gaps)<=min_gap_to_set && length(poss_range)>sum(bin_gaps<=min_gap_to_set)) && nrow(index_ranges)>=min_taxa_to_set)	{
				bin_lb[coll_no] <- min(poss_range[bin_gaps==min(bin_gaps)]);
				bin_ub[coll_no] <- max(poss_range[bin_gaps==min(bin_gaps)]);
				new_and_improved <- new_and_improved+1;
				paleodb_collections$interval_lb[coll_no] <- finest_chronostrat$interval[match(bin_lb[coll_no],finest_chronostrat$bin_first)];
				paleodb_collections$interval_ub[coll_no] <- finest_chronostrat$interval[match(bin_ub[coll_no],finest_chronostrat$bin_last)];
				paleodb_collections$ma_lb[coll_no] <- min(paleodb_collections$ma_lb[coll_no],finest_chronostrat$ma_lb[match(bin_lb[coll_no],finest_chronostrat$bin_first)]);
				paleodb_collections$ma_ub[coll_no] <- max(paleodb_collections$ma_ub[coll_no],finest_chronostrat$ma_ub[match(bin_ub[coll_no],finest_chronostrat$bin_last)]);
				# case where we've narrowed it down to one bin
				if (sum(bin_gaps==min_gap_to_set)==1)	{
					# the earliest possible first occurrences must be at least as old as the oldest possible first occurrences
					#	for those taxa first known from this interval
					sub_index_ranges <- subset(index_ranges,index_ranges$interval_lb==paleodb_collections$interval_lb[coll_no]);
					if (nrow(sub_index_ranges)>0)	{
						sub_index_ranges$ma_fa_mx[sub_index_ranges$ma_fa_mx<max(sub_index_ranges$ma_fa_mn)] <- max(sub_index_ranges$ma_fa_mn);
						index_ranges[match(rownames(sub_index_ranges),rownames(index_ranges)),] <- sub_index_ranges
						}
					
					# the latest possible last occurrences must be at least as young as the youngest possible last occurrences
					#	for those taxa last known from this interval
					sub_index_ranges <- subset(index_ranges,index_ranges$interval_ub==paleodb_collections$interval_ub[coll_no]);
					if (nrow(sub_index_ranges)>0)	{
						sub_index_ranges$ma_la_mn[sub_index_ranges$ma_la_mn>min(sub_index_ranges$ma_la_mx)] <- min(sub_index_ranges$ma_la_mx);
						index_ranges[match(rownames(sub_index_ranges),rownames(index_ranges)),] <- sub_index_ranges
						}
					indexed_species <- rownames(index_ranges);
					i_s <- match(indexed_species,rownames(minimum_range_data));
					minimum_range_data[i_s,] <- index_ranges;
					
					if (nrow(index_ranges_all) > nrow(index_ranges))	{
					# update species that do not have constrained finds yet
						unindexed_species <- rownames(index_ranges_all)[index_ranges_all$ma_fa_mn==0];
						u_s <- match(unindexed_species,rownames(minimum_range_data));
						sub_index_ranges <- subset(index_ranges,index_ranges$interval_lb==paleodb_collections$interval_lb[coll_no]);
						if (nrow(sub_index_ranges)>0)	{
							minimum_range_data$ma_fa_mx[u_s] <- min(sub_index_ranges$ma_fa_mx);
							minimum_range_data$ma_fa_mn[u_s] <- max(sub_index_ranges$ma_fa_mn);
							} else	{
							minimum_range_data$ma_fa_mx[u_s] <- paleodb_collections$ma_lb[coll_no];
							minimum_range_data$ma_fa_mn[u_s] <- paleodb_collections$ma_ub[coll_no];
							}
						sub_index_ranges <- subset(index_ranges,index_ranges$interval_ub==paleodb_collections$interval_ub[coll_no]);
						if (nrow(sub_index_ranges)>0)	{
							minimum_range_data$ma_la_mx[u_s] <- min(sub_index_ranges$ma_la_mx);
							minimum_range_data$ma_la_mn[u_s] <- max(sub_index_ranges$ma_la_mn);
							} else	{
							minimum_range_data$ma_la_mx[u_s] <- paleodb_collections$ma_lb[coll_no];
							minimum_range_data$ma_la_mn[u_s] <- paleodb_collections$ma_ub[coll_no];
							}
						minimum_range_data$interval_lb[u_s] <- paleodb_collections$interval_lb[coll_no];
						minimum_range_data$interval_ub[u_s] <- paleodb_collections$interval_ub[coll_no];
						}
					newly_fixed <- newly_fixed+1;
					}
				}
 			}
 		if (sum(bin_lb==Inf)>0 || sum(is.na(bin_lb))>0)
 			print(uc);
#		print(cbind(uc,minimum_range_data[2309,]));
#		print(c(paleodb_collections$ma_lb[195],paleodb_collections$ma_ub[195]));
		}
#	(1:ncolls)[is.na(as.character(paleodb_collections$interval_lb)!=as.character(paleodb_collections$interval_ub))];
	unfixed_collections <- paleodb_collections$collection_no[as.character(paleodb_collections$interval_lb)!=as.character(paleodb_collections$interval_ub)];
	unfixed <- length(unfixed_collections);
	improved <- new_and_improved;
	fixed <- newly_fixed;
	unfixed <- length(unfixed_collections);
	if (improved==0 && min_taxa_to_set>0)	{
		reboots <- reboots+1;
		if ((reboots %% 2)==0)	{
			min_taxa_to_set <- min_taxa_to_set-1;
			} else	{
			min_gap_to_set <- min_gap_to_set+1;
			}
		if (min_taxa_to_set>0)
			improved <-max(1,improved);
		}
	attempt <- attempt+1;
	fcolls <- hist(bin_ub-bin_lb,breaks=c(min(fuzz)-1,fuzz),plot=F)$counts;
#	fcolls <- hist(bin_ub-bin_lb,breaks=c(fuzz,max(fuzz)+1),plot=F)$counts;
	progress <- cbind(progress,fcolls);
	if (update_search)
		print(progress);
	}
return(paleodb_collections);
}

optimo_collection_and_occurrence_stratigraphy <- function(find_data,coll_data,chronostrat_units,lump_stages="")	{
# rescore collections if there is any lumping of reported stages into useful stages
ncolls <- nrow(coll_data);
if (lump_stages=="")
	lump_stages <- data.frame(reported=as.character(chronostrat_units$Stage),used=as.character(chronostrat_units$Stage),stringsAsFactors = FALSE);
if (sum(as.character(lump_stages$reported)!=as.character(lump_stages$used))>0)	{
	csu <- nrow(chronostrat_units);
	orig_stages <- nrow(lump_stages);
	Stage <- as.character(unique(lump_stages$used));
	used_stages <- length(Stage);
	Ma_LB <- Ma_UB <- vector(length=used_stages);
	
	ided_stages <- nrow(lump_stages);
	stages_to_lump <- (1:orig_stages)[as.character(lump_stages$reported)!=as.character(lump_stages$used)];
	Ma_LB[(1:used_stages)[Stage %in% chronostrat_units$Stage]] <- chronostrat_units$Ma_LB[(1:csu)[chronostrat_units$Stage %in% Stage]];
	Ma_UB[(1:used_stages)[Stage %in% chronostrat_units$Stage]] <- chronostrat_units$Ma_UB[(1:csu)[chronostrat_units$Stage %in% Stage]];

	composites <- (1:used_stages)[Ma_LB==0];

	for (i in 1:length(composites))	{
		comps <- (1:nrow(lump_stages))[lump_stages$used %in% Stage[composites[i]]];
		Ma_LB[composites[i]] <- max(chronostrat_units$Ma_LB[comps]);
		Ma_UB[composites[i]] <- min(chronostrat_units$Ma_UB[comps]);
		}
#	cbind(Stage,Ma_LB,Ma_UB)
	
	chronostrat_units <- data.frame(Stage,Ma_LB,Ma_UB);
	
	for (stl in 1:length(stages_to_lump))	{
		reported <- as.character(lump_stages$reported[stages_to_lump[stl]]);
		colls_to_lump_lb <- (1:ncolls)[coll_data$stage_lb %in% reported];
		colls_to_lump_ub <- (1:ncolls)[coll_data$stage_ub %in% reported];
		used <- as.character(lump_stages$used[stages_to_lump[stl]]);
		coll_data$stage_lb[colls_to_lump_lb] <- used;
		coll_data$stage_ub[colls_to_lump_ub] <- used;
		sum(is.na(coll_data$stage_lb))
		}
	}
nstages <- nrow(chronostrat_units);
	
coll_to_find_key <- match(find_data$collection_no,coll_data$collection_no);
ma_lb <- coll_data$ma_lb[coll_to_find_key];
ma_ub <- coll_data$ma_ub[coll_to_find_key];
stage_lb <- coll_data$stage_lb[coll_to_find_key];
stage_ub <- coll_data$stage_ub[coll_to_find_key];
find_data <- cbind(find_data,ma_lb,ma_ub,stage_lb,stage_ub);

noccr <- nrow(find_data);
taxa <- find_data$accepted_name;
taxa[find_data$accepted_rank=="genus"] <- find_data$identified_name[find_data$accepted_rank=="genus"];
taxa[find_data$accepted_rank=="subgenus"] <- find_data$identified_name[find_data$accepted_rank=="subgenus"];
find_data$accepted_name <- taxa;
taxon_names <- sort(unique(taxa));
ntaxa <- length(taxon_names);
taxon_no <- match(taxa,taxon_names);
find_data <- cbind(find_data,taxon_no);

bin_lb <- match(find_data$stage_lb,chronostrat_units$Stage);
bin_ub <- match(find_data$stage_ub,chronostrat_units$Stage);
fuzz <- bin_ub-bin_lb;
duds1a <- (1:noccr)[fuzz<0];
duds1b <- (1:noccr)[is.na(fuzz)];
while (length(duds1a)>0)	{
	stage_dummy <- find_data$stage_lb[duds1a];
	ma_lb_dummy <- find_data$ma_lb[duds1a];
	bin_lb_dummy <- bin_lb[duds1a];
	find_data$stage_lb[duds1a] <- find_data$stage_ub[duds1a];
	find_data$ma_lb[duds1a] <- find_data$ma_ub[duds1a];
	bin_lb[duds1a] <- bin_ub[duds1a];
	find_data$stage_ub[duds1a] <- stage_dummy;
	find_data$ma_ub[duds1a] <- ma_lb_dummy;
	bin_ub[duds1a] <- bin_lb_dummy;
	fuzz <- bin_ub-bin_lb;
	duds1a <- (1:noccr)[fuzz<0];
	}
find_data <- cbind(find_data,bin_lb,bin_ub,fuzz);

bin_lb <- match(coll_data$stage_lb,chronostrat_units$Stage);
bin_ub <- match(coll_data$stage_ub,chronostrat_units$Stage);
fuzz <- bin_ub-bin_lb;
duds2a <- (1:ncolls)[fuzz<0];
duds2b <- (1:ncolls)[is.na(fuzz)];
#cbind(coll_data$collection_no[duds2a],coll_data$formation[duds2a],coll_data$member[duds2a],coll_data$stage_lb[duds2a],coll_data$stage_ub[duds2a]);
#coll_data$collection_no[];
#write.csv(coll_data[duds2a,],paste(directory,"Problem_Collections.csv",sep=""),row.names=FALSE);
coll_data <- cbind(coll_data,bin_lb,bin_ub,fuzz);

#stage_coll <- array(0,nstages)
iteration <- improvement <- 1;
obs_fuzz <- sort(unique(fuzz));	# number of different ranges of uncertainty in collection age
obs_fuzz <- obs_fuzz[obs_fuzz>=0]
max_fuzz <- max(obs_fuzz);

inexactness <- array(0,dim=c(1+max(obs_fuzz),1));	# track number of collections with x bin uncertainty
rownames(inexactness) <- ((1:nrow(inexactness))-1);

for (dd in 1:(max_fuzz+1))	{
	d <- obs_fuzz[dd];
	inexactness[d+1,1] <- sum(coll_data$fuzz==obs_fuzz[dd]);
	}
occr_by_bin <- array(0,dim=c(ntaxa,nstages));
rock_unit_info <- cbind(coll_data$rock_no_sr,coll_data$rock_unit_senior);

ttl_taxon_occr <- array(0,ntaxa)
colnames(occr_by_bin) <- chronostrat_units$Stage;
rownames(occr_by_bin) <- taxon_names;

# date uncertain collections to minimize range extensions.
#	Repeat this until there is no improvement
fixed_coll <- subset(coll_data,coll_data$fuzz==0);
fixed_finds <- subset(find_data,find_data$fuzz==0);
taxon_no <- fixed_finds$taxon_no;
# redo taxon numbers
#taxon_no <- match(fixed_finds$taxon_no,sort(unique(fixed_finds$taxon_no)));
#fixed_finds$taxon_no <- taxon_no;	# added to prevent matrix overwrites!
stage_no <- fixed_finds$bin_lb;
sample_no <- fixed_finds$collection_no;
interval_names <- as.character(chronostrat_units$Stage);
#fixed_taxon_sampling_info <- numerare_concha_ex_tempore(taxon_no,stage_no,sample_no,taxon_names,interval_names);
occr_by_bin <- numerare_concha_ex_tempore(taxon_no,stage_no,sample_no,taxon_names,interval_names);
min_stage <- min(stage_no);
ttl_taxon_occr <- rowSums(occr_by_bin);
sample_age_lb <- fixed_finds$ma_lb;
sample_age_ub <- fixed_finds$ma_ub;
sampled_taxa <- fixed_finds$accepted_name;
#cbind(taxon_no, stage_no, sample_no, sample_age_lb, sample_age_ub)
taxon_stage_ranges <- sepkoskify_occurrence_data(taxon_no,stage_no,sample_no,sample_age_lb,sample_age_ub,sampled_taxa,interval_names,condense=FALSE);
if(nrow(taxon_stage_ranges) < length(taxon_names))	{
	ataxa <- length(taxon_names) - nrow(taxon_stage_ranges);
	dummy <- c(rep(0,4),rep("",2));
	for (i in 1:ataxa)
		taxon_stage_ranges <- rbind(taxon_stage_ranges,dummy);
#	ataxa <- nrow(taxon_stage_ranges) - length(taxon_names);
#	dummy <- c(rep(0,4),rep("",2))
#	taxon_stage_ranges <- rbind(taxon_stage_ranges,dummy);
	}
rownames(taxon_stage_ranges) <- taxon_names;
#taxon_stage_ranges$bin_lb <- frak_it(taxon_stage_ranges$bin_lb);
#sum(coll_data$fuzz<0);
duds_pre <- (1:ncolls)[coll_data$fuzz<0];
#coll_data$collection_no[duds_pre]
#write.csv(coll_data[duds_pre,],file="WTF3.csv",row.names=FALSE)
# first
while (improvement>0)	{
	### make note of which collections go from fuzz>0 to fuzz = 0
	###		separate out the occurrrences for those collections
	###		separate out taxa from those localities and update their range & sampling data.
	#### separate out unfixed localities: focus on these.
#	unfixed_coll <- subset(coll_data,coll_data$fuzz>0);
#	uncoll <- nrow(unfixed_coll);
#	unfixed_finds <- subset(find_data,find_data$fuzz>0);
#	unoccr <- nrow(unfixed_finds);

	print(paste("iteration",iteration,improvement));
	print(inexactness[,1:iteration]);

	if (iteration==1)	{
#		stage_coll <- array(0,dim=c(nstages,1))
		strt <- 1
		}	else	{
#		stage_coll <- cbind(stage_coll,stage_coll[,iteration-1])
		strt <- 2
		}
	inexactness <- cbind(inexactness,vector(length=(1+max_fuzz)));
#	for (dd in strt:length(obs_fuzz))	{
	dd <- 1;
	fixed_coll_info <- c();
	improved_coll_info <- c();
	while (dd < length(obs_fuzz))	{
		dd <- dd+1;		# skip obs_fuzz[1] = 0;
		d <- obs_fuzz[dd];
		# separate out both PaleoDB collection numbers & matrix row numbers
		relv_coll_no <- coll_data$collection_no[coll_data$fuzz==d];
		relv_coll <- (1:ncolls)[coll_data$fuzz==d];
		rc <- 0;
		while (rc < length(relv_coll))	{
			rc <- rc+1;
			relv_occr <- subset(find_data,find_data$collection_no==relv_coll_no[rc]);
			relv_occr <- relv_occr[(1:nrow(relv_occr))[match(unique(relv_occr$accepted_name),relv_occr$accepted_name)],];
			# fuzzy gives stage numbers of possible rock & collection age
			fuzzy <- c(frak_it(coll_data$bin_lb[relv_coll[rc]]),frak_it(coll_data$bin_ub[relv_coll[rc]]));
			fuzzy_range <- fuzzy[1]:fuzzy[2];
			assemblage <- relv_occr$taxon_no;
			index_ranges <- taxon_stage_ranges[assemblage,];
			unknowables <- subset(index_ranges,index_ranges$bin_lb==0);
			index_ranges <- subset(index_ranges,index_ranges$bin_lb>0);
			bin_imp <- bin_fx <- -1;
			if (nrow(index_ranges)>0)	{
				restriction <- c(max(frak_it(index_ranges$bin_lb)),min(frak_it(index_ranges$bin_ub)));
				restriction_range <- sort(restriction[1]:restriction[2]);
				# go through various cases of how well fossils restrict age of collection
				# redo this.  If fuzzy overlaps only partly with restriction
				#	reduce fuzzy that much.
				if (restriction[1]==restriction[2] && length(fuzzy_range[fuzzy_range %in% restriction_range])==1)	{
					### fixed
					bin_fx <- restriction[2];
#					fixed_coll_nos <- c(fixed_coll_nos,relv_coll_no[rc])
#					} else if (restriction[1] > restriction[2] && iteration > 2)	{
					} else if (restriction[1] > restriction[2])	{
					# latest first appearance is after earliest last appearance
					restriction <- sort(restriction);
					fuzzy_range <- fuzzy_range[fuzzy_range %in% restriction_range];
					if (length(fuzzy_range)==1)	{
						restriction[1] <- restriction[2] <- bin_fx <- fuzzy_range
						fuzzy <- c(min(fuzzy_range),max(fuzzy_range));
						} else if (length(fuzzy_range)==0)	{
						if (restriction[1]>fuzzy[2])	{
							bin_fx <- restriction[1] <- restriction[2] <- fuzzy[1] <- fuzzy[2];
							} else	{
							bin_fx <- restriction[1] <- restriction[2] <- fuzzy[2] <- fuzzy[1];
							}
						fuzzy_range <- fuzzy[1]:fuzzy[2];
						} else	{
						fuzzy[1] <- min(fuzzy_range);
						fuzzy[2] <- max(fuzzy_range);
						range <- cbind(as.numeric(index_ranges$bin_lb),as.numeric(index_ranges$bin_ub));
#					sapply(range,procrustean_binning,poss_range=sort(restriction),debug=TRUE);
						fits <- c();
						restriction <- sort(restriction);
						zzz <- rep(0,length(fuzzy_range));
						if (iteration > 3)	{
							for (rr in 1:nrow(range))	{
#								fits <- rbind(fits,procrustean_binning(range=range[rr,],poss_range=sort(restriction)));
								fits <- rbind(fits,procrustean_binning(range=range[rr,],poss_range=fuzzy));
								} # fits gives summed range extensions (in bins) we need to put the collection in this bin.
							xxx <- nrow(fits)+colSums(fits);
							if (min(xxx)<0)	xxx <- xxx-min(xxx);
							zzz <- xxx/max(1,sum(xxx));
							}
						if (max(zzz)>0)	{
							cutoff <- ((sum(xxx)/2)-(1/iteration))/sum(xxx);
							} else	{
							cutoff <- 0.0;
							}
						# shave off bad fits at extremes of possible ages
						while (iteration > 3 && length(zzz)>1 && min(zzz) < cutoff && (match(min(zzz),zzz)==1 || match(min(zzz),zzz)==length(zzz)))	{
							if (match(min(zzz),zzz)==1)	{
#								fuzzy[1] <- fuzzy[1]+1;
								fuzzy_range <- fuzzy_range[2:length(fuzzy_range)];
								zzz <- zzz[2:length(zzz)];
								} else {
#								fuzzy[2] <- fuzzy[2]-1;
								fuzzy_range <- fuzzy_range[1:(length(fuzzy_range)-1)];
								zzz <- zzz[1:((length(zzz)-1))];
								}
#							print(fuzzy);
							zzz <- zzz/sum(zzz);
							}
#						if ((length(zzz)==1 && iteration > 3) || length(fuzzy_range==1))	{
						if (length(fuzzy_range)==1)	{
							# if we have gone through 3+ iterations, then start assigning collections to best-fit bin
#							bin_fx <- fuzzy[1];
							bin_fx <- fuzzy[1] <- fuzzy[2] <- fuzzy_range[1];
							} else if (length(zzz)<(1+abs(fuzzy[2]-fuzzy[1])))	{
							bin_imp <- 1;
 							fuzzy <- restriction <- c(min(fuzzy_range),max(fuzzy_range));
							}
#						print(c(d,rc,relv_coll_no[rc]));
						}
					# hullo gives "weight" of each possible bin
					} else if (length((fuzzy[1]:fuzzy[2])[fuzzy[1]:fuzzy[2] %in% restriction[1]:restriction[2]])==0) {
					# if possible ages do not overlap with known ages, then assign the closest possible age
						if (fuzzy[2] < min(restriction))	{
					# if all species are younger, then assume youngest possible age
						bin_fx <- fuzzy[2];
#					fixed_coll_nos <- c(fixed_coll_nos,relv_coll_no[rc])
						} else if (fuzzy[1] > max(restriction))	{
					# if all species are older, then assume oldest possible age
						bin_fx <- fuzzy[1];
						}
#					fixed_coll_nos <- c(fixed_coll_nos,relv_coll_no[rc])
					} else if (restriction[1] < restriction[2]) {
					### see if taxa at least reduce "fuzziness" of collections at all
					if (fuzzy[1] < restriction[1])	{
					# earliest first appearance is after lower bound of uncertainty
						fuzzy[1] <- restriction[1];
						bin_imp <- 1;
						}
					if (fuzzy[2] > restriction[2])	{
					# latest last appearance is before upper bound of uncertainty
						fuzzy[2] <- restriction[2];
						bin_imp <- 1;
						}
					if (fuzzy[1]==fuzzy[2])	{
						bin_fx <- fuzzy[1];
						bin_imp <- 1;
						}
					} else if ((fuzzy[2]>=restriction[2] && fuzzy[1] <  restriction[1]) ||
							   (fuzzy[2]> restriction[2] && fuzzy[1] <= restriction[1])) {
					# case where taxa are known from only portion of possible rock range
					print(paste("this actually happened for dd =",dd,"& rc =",rc));
					if (restriction[1]==restriction[2])	{
						bin_fx <- restriction[1];
						bin_imp <- 1;
						} else	{
						fuzzy <- restriction;
						bin_imp <- 1;
						}
					}
				}
#			print(c(bin_fx,bin_imp));
#			taxon_stage_ranges[assemblage,]
			if (bin_fx != -1)	{
				# emend collection with bin_fx being the accepted bin
				fc <- match(relv_coll_no[rc],coll_data$collection_no);	# fixed collection
				fixed_coll_info <- rbind(fixed_coll_info,c(fc,bin_fx));
				coll_data$bin_lb[fc] <- coll_data$bin_ub[fc] <- bin_fx;
				coll_data$stage_lb[fc] <- coll_data$stage_ub[fc] <- as.character(chronostrat_units$Stage[bin_fx]);
				if (coll_data$ma_lb[fc] > frak_it(chronostrat_units$Ma_LB[bin_fx]))
					coll_data$ma_lb[fc] <- frak_it(chronostrat_units$Ma_LB[bin_fx]);
				if (coll_data$ma_ub[fc] < frak_it(chronostrat_units$Ma_UB[bin_fx]))
					coll_data$ma_ub[fc] <- frak_it(chronostrat_units$Ma_UB[bin_fx]);
				coll_data$fuzz[fc] <- 0;

				# emend strat ranges & find data
				scored_bin <- 1+bin_fx-min_stage;	# adjust for first bin being >1
				if (scored_bin <= ncol(occr_by_bin) && scored_bin>0)	{
					occr_by_bin[assemblage,scored_bin] <- occr_by_bin[assemblage,scored_bin]+1;
					ttl_taxon_occr[assemblage] <- ttl_taxon_occr[assemblage]+1;
					taxon_stage_ranges$bin_lb[assemblage[taxon_stage_ranges$bin_lb[assemblage] > bin_fx]] <- bin_fx;
					taxon_stage_ranges$bin_lb[assemblage[taxon_stage_ranges$bin_lb[assemblage] == 0]] <- bin_fx;
					taxon_stage_ranges$bin_ub[assemblage[taxon_stage_ranges$bin_ub[assemblage] < bin_fx]] <- bin_fx;
					newbies <- assemblage[taxon_stage_ranges$stage_lb[assemblage]==""];
					taxon_stage_ranges$stage_lb[assemblage]<- as.character(chronostrat_units$Stage[frak_it(taxon_stage_ranges$bin_lb[assemblage])]);
					taxon_stage_ranges$stage_ub[assemblage] <- as.character(chronostrat_units$Stage[frak_it(taxon_stage_ranges$bin_ub[assemblage])]);
					taxon_stage_ranges$ma_max[assemblage[taxon_stage_ranges$ma_max[assemblage] < coll_data$ma_lb[fc]]] <- frak_it(coll_data$ma_lb[fc]);
					taxon_stage_ranges$ma_min[assemblage[taxon_stage_ranges$ma_min[assemblage] > coll_data$ma_ub[fc]]] <- frak_it(coll_data$ma_ub[fc]);
					taxon_stage_ranges$ma_min[assemblage[taxon_stage_ranges$ma_min[assemblage] == 0]] <- frak_it(coll_data$ma_ub[fc]);
					coll_finds <- (1:noccr)[find_data$collection_no %in% relv_coll_no[rc]];
					find_data$bin_lb[coll_finds] <- frak_it(coll_data$bin_lb[fc]);
					find_data$bin_ub[coll_finds] <- frak_it(coll_data$bin_ub[fc]);
					find_data$ma_lb[coll_finds] <- frak_it(coll_data$ma_lb[fc]);
					find_data$ma_ub[coll_finds] <- frak_it(coll_data$ma_ub[fc]);
					find_data$stage_lb[coll_finds] <- coll_data$stage_lb[fc];
					find_data$stage_ub[coll_finds] <- coll_data$stage_ub[fc];
					find_data$fuzz[coll_finds] <- coll_data$fuzz[fc];
					}

#				taxon_stage_ranges[assemblage,]
				} else if (bin_imp==1)	{
#				print(c(dd,rc,relv_coll_no[rc]))
#				ic <- match(relv_coll_no[rc],coll_data$collection_no);	# fixed collection
				ic <- relv_coll[rc];
				improved_coll_info <- rbind(improved_coll_info,c(ic,fuzzy));
#				coll_data[ic,];
				# restrict the collection data
				if (coll_data$bin_lb[ic] < restriction[1])	{
					# collection is younger than currently allowed
					coll_data$bin_lb[ic] <- restriction[1];
					coll_data$stage_lb[ic] <- as.character(chronostrat_units$Stage[restriction[1]]);
#					if (coll_data$ma_lb[ic] > min(as.numeric(index_ranges$ma_max)));
#						coll_data$ma_lb[ic] <- min(as.numeric(index_ranges$ma_max));
					if (coll_data$ma_lb[ic] > chronostrat_units$Ma_LB[restriction[1]])
						coll_data$ma_lb[ic] <- chronostrat_units$Ma_LB[restriction[1]];
					}
				if (coll_data$bin_ub[ic] > restriction[2])	{
					# collection is older than currently allowed
					coll_data$bin_ub[ic] <- restriction[2];
					coll_data$stage_ub[ic] <- as.character(chronostrat_units$Stage[restriction[2]]);
#					if (coll_data$ma_ub[ic] > max(as.numeric(index_ranges$ma_min)))
#						coll_data$ma_ub[ic] <- max(index_ranges$ma_min);
					if (coll_data$ma_ub[ic] < chronostrat_units$Ma_UB[restriction[2]])
						coll_data$ma_ub[ic] <- chronostrat_units$Ma_UB[restriction[2]];
					}
				coll_data$fuzz[ic] <- coll_data$bin_ub[ic]-coll_data$bin_lb[ic];
				
				coll_finds <- (1:noccr)[find_data$collection_no %in% relv_coll_no[rc]];
				find_data$bin_lb[coll_finds] <- coll_data$bin_lb[ic];
				find_data$bin_ub[coll_finds] <- coll_data$bin_ub[ic];
				find_data$ma_lb[coll_finds] <- as.numeric(coll_data$ma_lb[ic]);
				find_data$ma_ub[coll_finds] <- as.numeric(coll_data$ma_ub[ic]);
				find_data$stage_lb[coll_finds] <- coll_data$stage_lb[ic];
				find_data$stage_ub[coll_finds] <- coll_data$stage_ub[ic];
				find_data$fuzz[coll_finds] <- coll_data$fuzz[ic];
				}
#			print(rc)
#			taxon_stage_ranges[assemblage,]
			}	# end case where we can reduce the viable bins
		} # end range of inexactness
	
#	fixed_coll <- subset(coll_data,coll_data$fuzz==0);
	obs_fuzz <- sort(unique(coll_data$fuzz[coll_data$fuzz>=0]));

	for (ddd in 1:(max_fuzz+1))	{
		inexactness[ddd,iteration+1] <- sum(coll_data$fuzz==(ddd-1));
		}
	improvement <- sum(abs(inexactness[,iteration]-inexactness[,iteration+1]));
	iteration <- iteration + 1;
#	test_taxon <- match("Leptaena (Ygdrasilomena) roomusoksi",rownames(taxon_stage_ranges))
#	print(taxon_stage_ranges[test_taxon,]);
	}
duds_post <- (1:ncolls)[coll_data$fuzz<0];
#sum(coll_data$fuzz<0);
#colSums(inexactness)
rock_sampling_info <- numerare_concha_ex_tempore(taxon_no,stage_no,sample_no,taxon_names="",interval_names="")
fixed_coll_data <- subset(coll_data,coll_data$fuzz==0);
fixed_coll_data <- subset(fixed_coll_data,fixed_coll_data$rock_no_sr>0);
#rock_ranges <- sepkoskify_occurrence_data_allow_unsampled(taxon_no=fixed_coll_data$rock_no_sr,stage_no=fixed_coll_data$bin_lb,sample_no=fixed_coll_data$collection_no,sample_age_lb = fixed_coll_data$ma_lb,sample_age_ub=fixed_coll_data$ma_ub,taxon_names=sort(unique(fixed_coll_data$rock_no_sr)));

unique_rocks <- sort(unique(fixed_coll_data$rock_no_sr));
unique_rocks_plus_names <- data.frame(rock_no_sr=as.numeric(frak_it(fixed_coll_data$rock_no_sr)),rock_unit_senior=as.character(fixed_coll_data$rock_unit_senior));
unique_rocks_plus_names <- unique_rocks_plus_names[order(unique_rocks_plus_names$rock_no_sr),]
keepers <- match(unique_rocks,unique_rocks_plus_names$rock_no_sr);
unique_rocks_plus_names <- unique_rocks_plus_names[keepers,];
unrocks <- length(unique_rocks);
rock_ranges <- data.frame(array(0,dim=c(max(unique_rocks),6)));
colnames(rock_ranges) <- c("rock_no_sr","ma_lb","ma_ub","bin_lb","bin_ub","fuzz");
for (ur in 1:unrocks)	{
	rock_ranges$rock_no_sr[ur] <- rock_no <- unique_rocks[ur];
	rock_coll <- subset(fixed_coll_data,fixed_coll_data$rock_no_sr==unique_rocks[ur]);
	rock_ranges$ma_lb[ur] <- max(rock_coll$ma_lb);
	rock_ranges$ma_ub[ur] <- min(rock_coll$ma_ub);
	rock_ranges$bin_ub[ur] <- min(rock_coll$bin_ub);
	rock_ranges$bin_lb[ur] <- max(rock_coll$bin_lb);
	rock_ranges$bin_ub[ur] <- min(rock_coll$bin_ub);
	rock_ranges$fuzz[ur] <- max(rock_coll$bin_lb)-min(rock_coll$bin_ub);
	}
rock_ranges_limiting <- subset(rock_ranges,rock_ranges$fuzz==0);
unique_rocks_plus_names <- cbind(unique_rocks_plus_names,rock_ranges[1:unrocks,]);

coll_data$fuzz <- coll_data$bin_ub - coll_data$bin_lb;
obs_fuzz <- sort(unique(coll_data$fuzz))
obs_fuzz <- obs_fuzz[obs_fuzz>0];

# coll_data <- saved_coll_data <- coll_data
# find_data <- saved_find_data <- find_data
# taxon_stage_ranges <- saved_taxon_stage_ranges <- taxon_stage_ranges
# occr_by_bin <- saved_occr_by_bin <- occr_by_bin
inexactness <- cbind(inexactness,vector(length=nrow(inexactness)));
iteration <- ncol(inexactness);
inexactness[1,iteration] <- inexactness[1,(iteration-1)];
for (d in 1:length(obs_fuzz))	{
	dd <- obs_fuzz[d];
	c_d <- subset(coll_data,coll_data$fuzz==dd);
	c_d <- subset(c_d,c_d$rock_no_sr>0);
	ncd <- nrow(c_d);
	# separate out collections known from rocks otherwise known from only one bin
	limited_colls <- (1:ncd)[c_d$rock_no_sr %in% rock_ranges_limiting$rock_no_sr]
	# now get the rock numbers
	limiting_rocks <- c_d$rock_no_sr[limited_colls];
	# match the limiting rocks to the rock_ranges thesaurus
	lrn <- match(limiting_rocks,rock_ranges_limiting$rock_no_sr);
	# match the overall collections to the particular collections being refined here.
	lrc <- (1:ncolls)[coll_data$collection_no %in% c_d$collection_no[limited_colls]];
	coll_data$bin_lb[lrc] <- rock_ranges_limiting$bin_lb[lrn];
	coll_data$bin_ub[lrc] <- rock_ranges_limiting$bin_ub[lrn];
	coll_data$stage_lb[lrc] <- as.character(chronostrat_units$Stage[rock_ranges_limiting$bin_lb[lrn]]);
	coll_data$stage_ub[lrc] <- as.character(chronostrat_units$Stage[rock_ranges_limiting$bin_ub[lrn]]);
	adj_these_l <- lrc[coll_data$ma_lb[lrc] > rock_ranges_limiting$ma_lb[lrn]];
	coll_data$ma_lb[adj_these_l] <- rock_ranges_limiting$ma_lb[match(coll_data$rock_no_sr[adj_these_l],rock_ranges_limiting$rock_no_sr)];
	adj_these_u <- lrc[coll_data$ma_ub[lrc] < rock_ranges_limiting$ma_ub[lrn]];
	coll_data$ma_ub[adj_these_u] <- rock_ranges_limiting$ma_ub[match(coll_data$rock_no_sr[adj_these_u],rock_ranges_limiting$rock_no_sr)];
	adj_these <- sort(unique(c(adj_these_l,adj_these_u)));
	coll_data$fuzz[adj_these] <- coll_data$bin_ub[adj_these]-coll_data$bin_lb[adj_these];
	# adjust finds
	adj_coll <- coll_data$collection_no[lrc];
	inexactness[1,iteration] <- sum(coll_data$fuzz==0);
	inexactness[(dd+1),iteration] <- sum(coll_data$fuzz==dd);
#	ac <- 0;
	for (ac in 1:length(adj_coll))	{
#	while (ac < length(adj_coll))	{
#		ac <- ac+1;
		adj_occurrence_no <- (1:noccr)[find_data$collection_no %in% adj_coll[ac]];
		if (length(adj_occurrence_no)>0)	{
			find_data$bin_lb[adj_occurrence_no] <- coll_data$bin_lb[lrc[ac]];
			find_data$bin_ub[adj_occurrence_no] <- coll_data$bin_ub[lrc[ac]];
			find_data$ma_lb[adj_occurrence_no] <- coll_data$ma_lb[lrc[ac]];
			find_data$ma_ub[adj_occurrence_no] <- coll_data$ma_ub[lrc[ac]];
			find_data$stage_lb[adj_occurrence_no] <- coll_data$stage_lb[lrc[ac]];
			find_data$stage_ub[adj_occurrence_no] <- coll_data$stage_ub[lrc[ac]];
			find_data$fuzz[adj_occurrence_no] <- coll_data$fuzz[lrc[ac]];
			tx_no <- match(find_data$accepted_name[adj_occurrence_no],rownames(taxon_stage_ranges));
			bin_no <- match(coll_data$stage_lb[lrc[ac]],colnames(occr_by_bin));
			occr_by_bin[tx_no,bin_no] <- occr_by_bin[tx_no,bin_no]+1;
			new_oldest_bin <- unique(c(tx_no[taxon_stage_ranges$bin_lb[tx_no]>coll_data$bin_lb[lrc[ac]]],
								tx_no[taxon_stage_ranges$bin_lb[tx_no]==0]));
			new_youngest_bin <- unique(c(tx_no[taxon_stage_ranges$bin_ub[tx_no]<coll_data$bin_ub[lrc[ac]]],
								  tx_no[taxon_stage_ranges$bin_ub[tx_no]==0]));
			new_oldest_ma <- unique(c(tx_no[taxon_stage_ranges$ma_max[tx_no]<coll_data$ma_lb[lrc[ac]]],
							   	tx_no[taxon_stage_ranges$ma_max[tx_no]==0]));
			new_youngest_ma <- unique(c(tx_no[taxon_stage_ranges$ma_min[tx_no]>coll_data$ma_ub[lrc[ac]]],
								 tx_no[taxon_stage_ranges$ma_min[tx_no]==0]));
#			print(c(length(new_oldest_bin),length(new_youngest_bin),length(new_oldest_ma),length(new_youngest_ma)))
			if (length(new_oldest_bin)>0)	{
				taxon_stage_ranges$bin_lb[new_oldest_bin] <- coll_data$bin_lb[lrc[ac]];
				taxon_stage_ranges$stage_lb[new_oldest_bin] <- coll_data$stage_lb[lrc[ac]];
				}
			if (length(new_youngest_bin)>0)	{
				taxon_stage_ranges$bin_ub[new_youngest_bin] <- coll_data$bin_ub[lrc[ac]];
				taxon_stage_ranges$stage_ub[new_youngest_bin] <- coll_data$stage_ub[lrc[ac]];
				}
			if (length(new_oldest_ma)>0)	{
#				taxon_stage_ranges$max_ma[new_oldest_ma] <- frak_it(taxon_stage_ranges$max_max[new_oldest_ma]);
				taxon_stage_ranges$ma_max[new_oldest_ma] <- frak_it(coll_data$ma_lb[lrc[ac]]);
				}
			if (length(new_youngest_ma)>0)	{
				taxon_stage_ranges$ma_min[new_youngest_ma] <- coll_data$ma_ub[lrc[ac]];
				}
			}
		}
#	improvement <- improvement+length(lrc);
#	d <- d+1
#	dd <- obs_fuzz[d];
	print(dd)
	}

improvement <- inexactness[1,ncol(inexactness)] - inexactness[1,(ncol(inexactness)-1)];
while (improvement>0)	{
	### make note of which collections go from fuzz>0 to fuzz = 0
	###		separate out the occurrrences for those collections
	###		separate out taxa from those localities and update their range & sampling data.
	#### separate out unfixed localities: focus on these.
#	unfixed_coll <- subset(coll_data,coll_data$fuzz>0);
#	uncoll <- nrow(unfixed_coll);
#	unfixed_finds <- subset(find_data,find_data$fuzz>0);
#	unoccr <- nrow(unfixed_finds);

	print(paste("iteration",iteration,improvement));
	print(inexactness[,1:iteration]);

	if (iteration==1)	{
#		stage_coll <- array(0,dim=c(nstages,1))
		strt <- 1
		}	else	{
#		stage_coll <- cbind(stage_coll,stage_coll[,iteration-1])
		strt <- 2
		}
	inexactness <- cbind(inexactness,vector(length=(1+max_fuzz)));
#	for (dd in strt:length(obs_fuzz))	{
	dd <- 1;
	fixed_coll_info <- c();
	improved_coll_info <- c();
	while (dd < length(obs_fuzz))	{
		dd <- dd+1;		# skip obs_fuzz[1] = 0;
		d <- obs_fuzz[dd];
		# separate out both PaleoDB collection numbers & matrix row numbers
		relv_coll_no <- coll_data$collection_no[coll_data$fuzz==d];
		relv_coll <- (1:ncolls)[coll_data$fuzz==d];
		rc <- 0;
		while (rc < length(relv_coll))	{
			rc <- rc+1;
			relv_occr <- subset(find_data,find_data$collection_no==relv_coll_no[rc]);
			relv_occr <- relv_occr[(1:nrow(relv_occr))[match(unique(relv_occr$accepted_name),relv_occr$accepted_name)],];
			# fuzzy gives stage numbers of possible rock & collection age
			fuzzy <- c(frak_it(coll_data$bin_lb[relv_coll[rc]]),frak_it(coll_data$bin_ub[relv_coll[rc]]));
			fuzzy_range <- fuzzy[1]:fuzzy[2];
			assemblage <- relv_occr$taxon_no;
			index_ranges <- taxon_stage_ranges[assemblage,];
			unknowables <- subset(index_ranges,index_ranges$bin_lb==0);
			index_ranges <- subset(index_ranges,index_ranges$bin_lb>0);
			bin_imp <- bin_fx <- -1;
			if (nrow(index_ranges)>0)	{
				restriction <- c(max(frak_it(index_ranges$bin_lb)),min(frak_it(index_ranges$bin_ub)));
				restriction_range <- sort(restriction[1]:restriction[2]);
				# go through various cases of how well fossils restrict age of collection
				# redo this.  If fuzzy overlaps only partly with restriction
				#	reduce fuzzy that much.
				if (restriction[1]==restriction[2] && length(fuzzy_range[fuzzy_range %in% restriction_range])==1)	{
					### fixed
					bin_fx <- restriction[2];
#					fixed_coll_nos <- c(fixed_coll_nos,relv_coll_no[rc])
#					} else if (restriction[1] > restriction[2] && iteration > 2)	{
					} else if (restriction[1] > restriction[2])	{
					# latest first appearance is after earliest last appearance
					restriction <- sort(restriction);
					fuzzy_range <- fuzzy_range[fuzzy_range %in% restriction_range];
					if (length(fuzzy_range)==1)	{
						restriction[1] <- restriction[2] <- bin_fx <- fuzzy_range
						fuzzy <- c(min(fuzzy_range),max(fuzzy_range));
						} else if (length(fuzzy_range)==0)	{
						if (restriction[1]>fuzzy[2])	{
							bin_fx <- restriction[1] <- restriction[2] <- fuzzy[1] <- fuzzy[2];
							} else	{
							bin_fx <- restriction[1] <- restriction[2] <- fuzzy[2] <- fuzzy[1];
							}
						fuzzy_range <- fuzzy[1]:fuzzy[2];
						} else	{
						fuzzy[1] <- min(fuzzy_range);
						fuzzy[2] <- max(fuzzy_range);
						range <- cbind(as.numeric(index_ranges$bin_lb),as.numeric(index_ranges$bin_ub));
#					sapply(range,procrustean_binning,poss_range=sort(restriction),debug=TRUE);
						fits <- c();
						restriction <- sort(restriction);
						zzz <- rep(0,length(fuzzy_range));
						if (iteration > 3)	{
							for (rr in 1:nrow(range))	{
#								fits <- rbind(fits,procrustean_binning(range=range[rr,],poss_range=sort(restriction)));
								fits <- rbind(fits,procrustean_binning(range=range[rr,],poss_range=fuzzy));
								} # fits gives summed range extensions (in bins) we need to put the collection in this bin.
							xxx <- nrow(fits)+colSums(fits);
							if (min(xxx)<0)	xxx <- xxx-min(xxx);
							zzz <- xxx/max(1,sum(xxx));
							}
						if (max(zzz)>0)	{
							cutoff <- ((sum(xxx)/2)-(1/iteration))/sum(xxx);
							} else	{
							cutoff <- 0.0;
							}
						# shave off bad fits at extremes of possible ages
						while (iteration > 3 && length(zzz)>1 && min(zzz) < cutoff && (match(min(zzz),zzz)==1 || match(min(zzz),zzz)==length(zzz)))	{
							if (match(min(zzz),zzz)==1)	{
#								fuzzy[1] <- fuzzy[1]+1;
								fuzzy_range <- fuzzy_range[2:length(fuzzy_range)];
								zzz <- zzz[2:length(zzz)];
								} else {
#								fuzzy[2] <- fuzzy[2]-1;
								fuzzy_range <- fuzzy_range[1:(length(fuzzy_range)-1)];
								zzz <- zzz[1:((length(zzz)-1))];
								}
#							print(fuzzy);
							zzz <- zzz/sum(zzz);
							}
#						if ((length(zzz)==1 && iteration > 3) || length(fuzzy_range==1))	{
						if (length(fuzzy_range)==1)	{
							# if we have gone through 3+ iterations, then start assigning collections to best-fit bin
#							bin_fx <- fuzzy[1];
							bin_fx <- fuzzy[1] <- fuzzy[2] <- fuzzy_range[1];
							} else if (length(zzz)<(1+abs(fuzzy[2]-fuzzy[1])))	{
							bin_imp <- 1;
 							fuzzy <- restriction <- c(min(fuzzy_range),max(fuzzy_range));
							}
#						print(c(d,rc,relv_coll_no[rc]));
						}
					# hullo gives "weight" of each possible bin
					} else if (length((fuzzy[1]:fuzzy[2])[fuzzy[1]:fuzzy[2] %in% restriction[1]:restriction[2]])==0) {
					# if possible ages do not overlap with known ages, then assign the closest possible age
						if (fuzzy[2] < min(restriction))	{
					# if all species are younger, then assume youngest possible age
						bin_fx <- fuzzy[2];
#					fixed_coll_nos <- c(fixed_coll_nos,relv_coll_no[rc])
						} else if (fuzzy[1] > max(restriction))	{
					# if all species are older, then assume oldest possible age
						bin_fx <- fuzzy[1];
						}
#					fixed_coll_nos <- c(fixed_coll_nos,relv_coll_no[rc])
					} else if (restriction[1] < restriction[2]) {
					### see if taxa at least reduce "fuzziness" of collections at all
					if (fuzzy[1] < restriction[1])	{
					# earliest first appearance is after lower bound of uncertainty
						fuzzy[1] <- restriction[1];
						bin_imp <- 1;
						}
					if (fuzzy[2] > restriction[2])	{
					# latest last appearance is before upper bound of uncertainty
						fuzzy[2] <- restriction[2];
						bin_imp <- 1;
						}
					if (fuzzy[1]==fuzzy[2])	{
						bin_fx <- fuzzy[1];
						bin_imp <- 1;
						}
					} else if ((fuzzy[2]>=restriction[2] && fuzzy[1] <  restriction[1]) ||
							   (fuzzy[2]> restriction[2] && fuzzy[1] <= restriction[1])) {
					# case where taxa are known from only portion of possible rock range
					print(paste("this actually happened for dd =",dd,"& rc =",rc));
					if (restriction[1]==restriction[2])	{
						bin_fx <- restriction[1];
						bin_imp <- 1;
						} else	{
						fuzzy <- restriction;
						bin_imp <- 1;
						}
					}
				}
#			print(c(bin_fx,bin_imp));
#			taxon_stage_ranges[assemblage,]
			if (bin_fx != -1)	{
				# emend collection with bin_fx being the accepted bin
				fc <- match(relv_coll_no[rc],coll_data$collection_no);	# fixed collection
				fixed_coll_info <- rbind(fixed_coll_info,c(fc,bin_fx));
				coll_data$bin_lb[fc] <- coll_data$bin_ub[fc] <- bin_fx;
				coll_data$stage_lb[fc] <- coll_data$stage_ub[fc] <- as.character(chronostrat_units$Stage[bin_fx]);
				if (coll_data$ma_lb[fc] > frak_it(chronostrat_units$Ma_LB[bin_fx]))
					coll_data$ma_lb[fc] <- frak_it(chronostrat_units$Ma_LB[bin_fx]);
				if (coll_data$ma_ub[fc] < frak_it(chronostrat_units$Ma_UB[bin_fx]))
					coll_data$ma_ub[fc] <- frak_it(chronostrat_units$Ma_UB[bin_fx]);
				coll_data$fuzz[fc] <- 0;

				# emend strat ranges & find data
				scored_bin <- 1+bin_fx-min_stage;	# adjust for first bin being >1
				if (scored_bin <= ncol(occr_by_bin) && scored_bin>0)	{
					occr_by_bin[assemblage,scored_bin] <- occr_by_bin[assemblage,scored_bin]+1;
					ttl_taxon_occr[assemblage] <- ttl_taxon_occr[assemblage]+1;
					taxon_stage_ranges$bin_lb[assemblage[taxon_stage_ranges$bin_lb[assemblage] > bin_fx]] <- bin_fx;
					taxon_stage_ranges$bin_lb[assemblage[taxon_stage_ranges$bin_lb[assemblage] == 0]] <- bin_fx;
					taxon_stage_ranges$bin_ub[assemblage[taxon_stage_ranges$bin_ub[assemblage] < bin_fx]] <- bin_fx;
					newbies <- assemblage[taxon_stage_ranges$stage_lb[assemblage]==""];
					taxon_stage_ranges$stage_lb[assemblage]<- as.character(chronostrat_units$Stage[frak_it(taxon_stage_ranges$bin_lb[assemblage])]);
					taxon_stage_ranges$stage_ub[assemblage] <- as.character(chronostrat_units$Stage[frak_it(taxon_stage_ranges$bin_ub[assemblage])]);
					taxon_stage_ranges$ma_max[assemblage[taxon_stage_ranges$ma_max[assemblage] < coll_data$ma_lb[fc]]] <- frak_it(coll_data$ma_lb[fc]);
					taxon_stage_ranges$ma_min[assemblage[taxon_stage_ranges$ma_min[assemblage] > coll_data$ma_ub[fc]]] <- frak_it(coll_data$ma_ub[fc]);
					taxon_stage_ranges$ma_min[assemblage[taxon_stage_ranges$ma_min[assemblage] == 0]] <- frak_it(coll_data$ma_ub[fc]);
					coll_finds <- (1:noccr)[find_data$collection_no %in% relv_coll_no[rc]];
					find_data$bin_lb[coll_finds] <- frak_it(coll_data$bin_lb[fc]);
					find_data$bin_ub[coll_finds] <- frak_it(coll_data$bin_ub[fc]);
					find_data$ma_lb[coll_finds] <- frak_it(coll_data$ma_lb[fc]);
					find_data$ma_ub[coll_finds] <- frak_it(coll_data$ma_ub[fc]);
					find_data$stage_lb[coll_finds] <- coll_data$stage_lb[fc];
					find_data$stage_ub[coll_finds] <- coll_data$stage_ub[fc];
					find_data$fuzz[coll_finds] <- coll_data$fuzz[fc];
					}

#				taxon_stage_ranges[assemblage,]
				} else if (bin_imp==1)	{
#				print(c(dd,rc,relv_coll_no[rc]))
#				ic <- match(relv_coll_no[rc],coll_data$collection_no);	# fixed collection
				ic <- relv_coll[rc];
				improved_coll_info <- rbind(improved_coll_info,c(ic,fuzzy));
#				coll_data[ic,];
				# restrict the collection data
				if (coll_data$bin_lb[ic] < restriction[1])	{
					# collection is younger than currently allowed
					coll_data$bin_lb[ic] <- restriction[1];
					coll_data$stage_lb[ic] <- as.character(chronostrat_units$Stage[restriction[1]]);
#					if (coll_data$ma_lb[ic] > min(as.numeric(index_ranges$ma_max)));
#						coll_data$ma_lb[ic] <- min(as.numeric(index_ranges$ma_max));
					if (coll_data$ma_lb[ic] > chronostrat_units$Ma_LB[restriction[1]])
						coll_data$ma_lb[ic] <- chronostrat_units$Ma_LB[restriction[1]];
					}
				if (coll_data$bin_ub[ic] > restriction[2])	{
					# collection is older than currently allowed
					coll_data$bin_ub[ic] <- restriction[2];
					coll_data$stage_ub[ic] <- as.character(chronostrat_units$Stage[restriction[2]]);
#					if (coll_data$ma_ub[ic] > max(as.numeric(index_ranges$ma_min)))
#						coll_data$ma_ub[ic] <- max(index_ranges$ma_min);
					if (coll_data$ma_ub[ic] < chronostrat_units$Ma_UB[restriction[2]])
						coll_data$ma_ub[ic] <- chronostrat_units$Ma_UB[restriction[2]];
					}
				coll_data$fuzz[ic] <- coll_data$bin_ub[ic]-coll_data$bin_lb[ic];
				
				coll_finds <- (1:noccr)[find_data$collection_no %in% relv_coll_no[rc]];
				find_data$bin_lb[coll_finds] <- coll_data$bin_lb[ic];
				find_data$bin_ub[coll_finds] <- coll_data$bin_ub[ic];
				find_data$ma_lb[coll_finds] <- as.numeric(coll_data$ma_lb[ic]);
				find_data$ma_ub[coll_finds] <- as.numeric(coll_data$ma_ub[ic]);
				find_data$stage_lb[coll_finds] <- coll_data$stage_lb[ic];
				find_data$stage_ub[coll_finds] <- coll_data$stage_ub[ic];
				find_data$fuzz[coll_finds] <- coll_data$fuzz[ic];
				}
#			print(rc)
#			taxon_stage_ranges[assemblage,]
			}	# end case where we can reduce the viable bins
		} # end range of inexactness
	
#	fixed_coll <- subset(coll_data,coll_data$fuzz==0);
	obs_fuzz <- sort(unique(coll_data$fuzz[coll_data$fuzz>=0]));

	for (ddd in 1:(max_fuzz+1))	{
		inexactness[ddd,iteration+1] <- sum(coll_data$fuzz==(ddd-1));
		}
	improvement <- sum(abs(inexactness[,iteration]-inexactness[,iteration+1]));
	iteration <- iteration + 1;
#	test_taxon <- match("Leptaena (Ygdrasilomena) roomusoksi",rownames(taxon_stage_ranges))
#	print(taxon_stage_ranges[test_taxon,]);
	}

coll_data$fuzz <- coll_data$bin_ub-coll_data$bin_lb;
rock_no_sr <- coll_data$rock_no_sr[coll_to_find_key];
rock_no <- coll_data$rock_no[coll_to_find_key];
formation_no <- coll_data$formation_no[coll_to_find_key];
rock_unit_senior <- coll_data$rock_unit_senior[coll_to_find_key];
find_data <- cbind(find_data,rock_no_sr,rock_no,formation_no,rock_unit_senior);

output <- list(coll_data,find_data,unique_rocks_plus_names,taxon_stage_ranges,occr_by_bin);
names(output) <- c("Optimized_Collections","Optimized_Occurrences","Optimized_Rock_Units","Optimized_Taxon_Ranges","Optimized_Taxon_Finds_by_Bin");
return(output);
}

numerare_concha_ex_tempore <- function(taxon_no,stage_no,sample_no,taxon_names="",interval_names="")	{
# use unique combinations of taxa, collections & stages to eliminate 2+ versions of same taxon due to synonymization
# all_taxa: if false, then 
fossil_record <- data.frame(unique(cbind(taxon_no,stage_no,sample_no)));
if (length(taxon_names) == max(taxon_no))	{
#if (all_taxa)	{
	unique_taxa <- 1:max(max(taxon_no),length(taxon_names))
	} else	{
	unique_taxa <- sort(unique(taxon_no));
	}
ntaxa <- length(unique_taxa);
nstages <- length(unique(stage_no));
max_stage <- max(stage_no);
min_stage <- min(stage_no);
stage_finds <- array(0,dim=c(ntaxa,1+(max_stage-min_stage)));
if (length(taxon_names)>1)	rownames(stage_finds) <- taxon_names[unique_taxa];
if (length(interval_names)>1)	colnames(stage_finds) <- interval_names[min_stage:max_stage];
for (otu in 1:ntaxa)	{
	notu <- unique_taxa[otu];
	taxon_record <- unique(subset(fossil_record,fossil_record$taxon_no==notu));
	bins_pres <- sort(unique(taxon_record$stage_no));
	bp <- match(bins_pres,min_stage:max_stage);
	for (b in 1:length(bins_pres))
		stage_finds[otu,bp[b]] <- sum(taxon_record$stage_no==bins_pres[b]);
	}
cn <- colnames(stage_finds);
if (length(taxon_names)==nrow(stage_finds)) {
	stage_finds <- data.frame(stage_finds,row.names = as.character(taxon_names));
	} else	stage_finds <- data.frame(stage_finds);
colnames(stage_finds) <- cn;
return(stage_finds);
}

## count fossil finds per time unit
numerare_concha_ex_tempore_paleodb <- function(pbdb_finds,interval_names="",coll_unit="collections")	{
# use unique combinations of taxa, collections & stages to eliminate 2+ versions of same taxon due to synonymization
# all_taxa: if false, then 

notus <- sort(unique(pbdb_finds$taxon_no));
ntaxa <- length(notus);
setup_names <- unique(cbind(pbdb_finds$taxon_no,pbdb_finds$accepted_name));
setup_names <- setup_names[order(frak_it(setup_names[,1])),];
taxon_names <- as.character(setup_names[,2]);
min_stage <- min(pbdb_finds$bin_lb);
max_stage <- max(pbdb_finds$bin_ub);
#setup_stages <- unique(cbind(pbdb_finds$bin_lb,pbdb_finds$stage_lb))
#setup_stages <- setup_stages[order(frak_it(setup_stages[,1])),];

stage_finds <- array(0,dim=c(ntaxa,max_stage));
for (otu in 1:ntaxa)	{
	taxon_finds <- subset(pbdb_finds,pbdb_finds$taxon_no==notus[otu]);
	if (coll_unit=="collections")	{
		samples <- unique(cbind(taxon_finds$bin_lb,taxon_finds$collection_no));
		} else if (coll_unit=="rock units" || coll_unit=="rock_units")	{
		samples <- unique(cbind(taxon_finds$bin_lb,taxon_finds$rock_no_sr));
		} else if (coll_unit=="formation")	{
		samples <- unique(cbind(taxon_finds$bin_lb,taxon_finds$formation_no));
		}
	lbin <- min(samples[,1]);
	ubin <- max(samples[,1]);
	for (bb in lbin:ubin)	stage_finds[otu,bb] <- sum(samples[,1]==bb)
	}
rownames(stage_finds) <- taxon_names;
if (length(interval_names)>=ncol(stage_finds))
	colnames(stage_finds) <- interval_names[1:max_stage];
stage_finds <- stage_finds[,min_stage:max_stage];
return(stage_finds);
}

## count collection units (localities, formations, etc.)
numerare_contigere_ex_tempore <- function(stage_no,sample_no,interval_names="")	{
fossil_record <- data.frame(unique(cbind(sample_no,stage_no)));
fossil_record <- subset(fossil_record,fossil_record$sample_no>0);
unique_stages <- sort(unique(stage_no));
stage_samples <- array(0,dim=c(1,length(unique_stages)));
colnames(stage_samples) <- interval_names[unique_stages];
for (us in 1:length(unique_stages))	{
	stage_samples[1,us] <- sum(fossil_record$stage_no==unique_stages[us]);
	}
return(stage_samples);
}

old_and_slow <- function()	{
while (improvement>0)	{
	unfixed_coll <- subset(coll_data,coll_data$fuzz>0);
	unfixed_finds <- subset(find_data,find_data$fuzz>0);
	print(paste("iteration",iteration));
	if (iteration==1)	{
		stage_coll <- array(0,dim=c(nstages,1))
		strt <- 1
		}	else	{
		stage_coll <- cbind(stage_coll,stage_coll[,iteration-1])
		strt <- 2
		}
	inexactness <- cbind(inexactness,vector(length=length(obs_fuzz)));
#	for (dd in strt:length(obs_fuzz))	{
	dd <- 0;
	while (dd < length(obs_fuzz))	{
		dd <- dd+1;
		d <- obs_fuzz[dd];
		# separate out both PaleoDB collection numbers & matrix row numbers
		relv_coll_no <- coll_data$collection_no[fuzz==d];
		relv_coll <- (1:ncoll)[fuzz==d];
		rc <- 0;
		while (rc < length(relv_coll))	{
			rc <- rc+1;
			relv_occr <- subset(find_data,find_data$collection_no==relv_coll_no[rc]);
			hullo <- array(0,(d+1));
			## do two loops.  One to tie down age of collection (if needed).  Second to update species info
			if (d>0)	{
				ss <- 0;
				while (ss < nrow(relv_occr))	{
					ss <- ss+1;
					nspc <- relv_occr$taxon_no[ss];
					if (taxon_stage_ranges[nspc,2]>0 && taxon_stage_ranges[nspc,1]>0)	{
						fuzzy <- (coll_data$bin_lb[relv_coll[rc]]:coll_data$bin_ub[relv_coll[rc]]);
						if (taxon_stage_ranges[nspc,2]>=max(fuzzy) && taxon_stage_ranges[nspc,1]>min(fuzzy))	{
						# case where species appears within this fuzzy zone
							for (g in 1:(d+1))
								if (fuzzy[g] < taxon_stage_ranges[nspc,1])	
									hullo[g] <- hullo[g]+abs(fuzzy[g]-taxon_stage_ranges[nspc,1]);
							} else if (taxon_stage_ranges[nspc,1]<=min(fuzzy) && taxon_stage_ranges[nspc,1]<max(fuzzy))	{
						# case where species disappears within this fuzzy zone
							for (g in 1:(d+1))
								if (fuzzy[g] > taxon_stage_ranges[nspc,2])	
									hullo[g] <- hullo[g]+abs(taxon_stage_ranges[nspc,2]-fuzzy[g]);
							} else if (taxon_stage_ranges[nspc,1]>min(fuzzy) && taxon_stage_ranges[nspc,2]<max(fuzzy))	{
					# case where known range is within possible range of collection
						for (g in 1:(d+1))
							if (fuzzy[g] < taxon_stage_ranges[nspc,1])	{
								hullo[g] <- hullo[g]+abs(fuzzy[g]-taxon_stage_ranges[nspc,1])
								} else if (fuzzy[g] > taxon_stage_ranges[nspc,2])	{
								hullo[g] <- hullo[g]+abs(fuzzy[g]-taxon_stage_ranges[nspc,2])
								}
							}
						}
					# go through species and see how much distortion of ranges is necessary
					}
				# end case of inexact stratigraphic placement of collection
				}	else	{	
				fuzzy <- array(coll_data$bin_lb[relv_coll[rc]],dim=1)
				}
			bin <- fuzzy[hullo %in% min(hullo)];
			if (d==0 || length(bin)==1)	{
				ss <- 0;
				while (ss < nrow(relv_occr))	{
					ss <- ss+1;
					nspc <- relv_occr$taxon_no[ss];
					occr_by_bin[nspc,bin] <- occr_by_bin[nspc,bin]+1;
					ttl_taxon_occr[nspc] <- ttl_taxon_occr[nspc]+1;
					if (ttl_taxon_occr[nspc]==1)	{
						taxon_stage_ranges[nspc,1] <- taxon_stage_ranges[nspc,2] <- bin
						}	else if (bin < taxon_stage_ranges[nspc,1])	{
						taxon_stage_ranges[nspc,1] <- bin
						}	else if (bin > taxon_stage_ranges[nspc,2])	{
						taxon_stage_ranges[nspc,2] <- bin
						}
					}
				stage_coll[bin,iteration] <- stage_coll[bin,iteration]+1
#				if (!is.na(coll_data$rock_unit_senior[relv_coll[rc]]))	{
#					fn <- match(coll_data$formation_no[relv_coll[rc]],form_nos);
#					forms_by_bin[fn,bin] <- 1;
#					rn <- match(coll_data$rock_no_sr[relv_coll[rc]],rock_unit_nos);
#					rocks_by_bin[rn,bin] <- 1;
#					}
				coll_data$bin_lb[relv_coll[rc]] <- coll_data$bin_ub[relv_coll[rc]] <- bin
				coll_data$stage_lb[relv_coll[rc]] <- coll_data$stage_ub[relv_coll[rc]] <- as.character(chronostrat_units$Stage[bin]);
				if (coll_data$ma_lb[relv_coll[rc]] > chronostrat_units$Ma_LB[bin])
					coll_data$ma_lb[relv_coll[rc]] <- chronostrat_units$Ma_LB[bin];
				if (coll_data$ma_ub[relv_coll[rc]] < chronostrat_units$Ma_UB[bin])
					coll_data$ma_ub[relv_coll[rc]] <- chronostrat_units$Ma_UB[bin];
				coll_data$fuzz[relv_coll[rc]] <- fuzz[relv_coll[rc]] <- 0;
				# end case of single best fit
				}	else if (length(bin)<=d && sum(bin[2:length(bin)]-bin[1:(length(bin)-1)])==(length(bin)-1))	{
				coll_data$bin_lb[relv_coll[rc]] <- min(bin);
				coll_data$bin_ub[relv_coll[rc]] <- max(bin);
				#### adjust other strat info
				coll_data$stage_lb[relv_coll[rc]] <- as.character(chronostrat_units$Stage[min(bin)]); 
				coll_data$stage_ub[relv_coll[rc]] <- as.character(chronostrat_units$Stage[max(bin)]);
				if (coll_data$ma_lb[relv_coll[rc]] > max(chronostrat_units$Ma_LB[bin]))
					coll_data$ma_lb[relv_coll[rc]] <- max(chronostrat_units$Ma_LB[bin]);
				if (coll_data$ma_ub[relv_coll[rc]] < min(chronostrat_units$Ma_UB[bin]))
					coll_data$ma_ub[relv_coll[rc]] <- max(chronostrat_units$Ma_UB[bin]);
				coll_data$fuzz[relv_coll[rc]] <- fuzz[relv_coll[rc]] <- max(bin)-min(bin);
				}
#			print(rc);
			}	# end case where we can reduce the viable bins
		} # end range of inexactness
	for (ddd in 1:length(obs_fuzz))
		inexactness[ddd,iteration+1] <- sum(fuzz==obs_fuzz[ddd]);

#	obs_fuzz <- sort(unique(fuzz));	# number of different ranges of uncertainty in collection age
#	obs_fuzz <- obs_fuzz[obs_fuzz>=0]

	improvement <- sum(abs(inexactness[,iteration]-inexactness[,iteration+1]));
	iteration <- iteration + 1;
	}
}

sepkoskify_occurrence_data_old <- function(taxon_no,stage_no,sample_no,sample_age_lb,sample_age_ub,taxon_names="",interval_names="")  {
fossil_record <- data.frame(taxon_no,stage_no,sample_no,sample_age_lb,sample_age_ub);
unique_taxa <- sort(unique(taxon_no));
ntaxa <- length(unique_taxa);
compendium <- c();
for (otu in 1:ntaxa)	{
	notu <- unique_taxa[otu];
	taxon_record <- unique(subset(fossil_record,fossil_record$taxon_no==notu));
	dim(taxon_record)
	ma_max <- max(taxon_record$sample_age_lb);
	ma_min <- min(taxon_record$sample_age_ub);
	bin_lb <- min(taxon_record$stage_no);
	bin_ub <- max(taxon_record$stage_no);
	if (length(interval_names)>0)	{
		stage_lb <- interval_names[bin_lb];
		stage_ub <- interval_names[bin_ub];
		compendium <- rbind(compendium,c(stage_lb,stage_ub,bin_lb,bin_ub,ma_max,ma_min));
		} else	{
		compendium <- rbind(compendium,c(bin_lb,bin_ub,ma_max,ma_min));
		}
	}
if (length(interval_names)>0)	{
	colnames(compendium) <- c("stage_lb","stage_ub","bin_lb","bin_ub","ma_max","ma_min");
	} else	{
	colnames(compendium) <- c("bin_lb","bin_ub","ma_max","ma_min");
	}
if (length(taxon_names)==length(unique_taxa))	rownames(compendium) <- taxon_names;
return(data.frame(compendium));
}

# find gaps within taxon ranges
revelio_range_gaps_in_finds_per_bin <- function(taxon_finds_per_bin)	{
nbins <- length(taxon_finds_per_bin);
bin_finds <- (1:nbins)[taxon_finds_per_bin>0];
output <- c(0,0);
if (length(bin_finds)>1)	{
	ind_gaps <- c();
	for (f in 2:length(bin_finds))
		ind_gaps <- c(ind_gaps,abs(bin_finds[f]-bin_finds[f-1])-1);
	output <- (c(max(ind_gaps),bin_finds[match(max(ind_gaps),ind_gaps)]));
	}
return(output);
}

# split taxa with overlong gaps into two taxa
diffindo_gaps_within_ranges <- function(finds_per_bin,taxon_ranges="",bins="",taxa="",max_gap=3)	{
ntaxa <- nrow(finds_per_bin);
nbins <- ncol(finds_per_bin);
if (length(taxa) < ntaxa)	taxa <- rownames(finds_per_bin);
if (length(bins) < nbins)	bins <- colnames(finds_per_bin); 
jack_fix <- FALSE;
min_bin <- 0
if (is.matrix(taxon_ranges) || is.data.frame(taxon_ranges))	{
	jack_fix <- TRUE;
	min_bin <- min(taxon_ranges)-1;
	}
#taxon_finds_per_bin <- (finds_per_bin);
#gps <- base::t(sapply(taxon_finds_per_bin,revelio_range_gaps_in_finds_per_bin));
#gps[,2]
gappies <- 1;

while (length(gappies) > 0 )	{
	gps <- c()
	for (tx in 1:ntaxa)	{
		taxon_finds_per_bin <- finds_per_bin[tx,];
	 	gps <- rbind(gps,revelio_range_gaps_in_finds_per_bin(taxon_finds_per_bin))
		}
	rownames(gps) <- rownames(finds_per_bin);
	gappies <- (1:ntaxa)[gps[,1]>=max_gap];
	if (length(gappies)>0)	{
#		finds_per_bin_alt <- finds_per_bin;
		gappies <- sort(gappies,decreasing = TRUE);
		for (gs in 1:length(gappies))	{
			gg <- gappies[gs];
			finds_per_bin_alt1 <- finds_per_bin_alt2 <- finds_per_bin[gg,]
			if (jack_fix)	{
				new_ranges_1 <- new_ranges_2 <- taxon_ranges[gg,]
				new_ranges_1[2] <- gps[gg,2]+min_bin;
				new_ranges_2[1] <- gps[gg,2]+gps[gg,1]+1+min_bin;
				new_ranges <- rbind(new_ranges_1,new_ranges_2);
				rownames(new_ranges) <- c(rownames(taxon_ranges)[gg],paste(rownames(taxon_ranges)[gg],"1",sep=""));
				taxon_ranges <- xxxx <- rbind(taxon_ranges[1:(gg-1),],new_ranges,taxon_ranges[(gg+1):ntaxa,]);
				}
			finds_per_bin_alt2[1:(gps[gg,2])] <- 0;
			finds_per_bin_alt1[(gps[gg,2]+1):nbins] <- 0;
			finds_per_bin_alt <- rbind(finds_per_bin_alt1,finds_per_bin_alt2)
			finds_per_bin <- xxx <- rbind(finds_per_bin[1:(gg-1),],finds_per_bin_alt,finds_per_bin[(gg+1):ntaxa,])
			ntaxa <- nrow(finds_per_bin);
			}
		}
	}
if (jack_fix)	{
	output <- list(finds_per_bin,taxon_ranges);
	names(output) <- c("finds_per_bin","taxon_ranges");
	} else	{
	output <- list(finds_per_bin);
	names(output) <- c("finds_per_bin");
	}

return(output);
}

# get basic numbers for boundary crossing
summarize_data_for_boundary_crossing_analysis <- function(finds_per_bin,taxon_ranges,interval_names="")	{
ntaxa <- nrow(taxon_ranges)
## make sure that the first column of finds_per_bin corresponds to bin 1 in taxon_ranges
a_fpb <- min((1:ncol(finds_per_bin))[colSums(finds_per_bin)>0]);
stage_1 <- min(taxon_ranges);
if (a_fpb < stage_1)	{
	dummy <- array("",dim=c(ntaxa,stage_1-a_fpb));
	colnames(dummy) <- rep("",stage_1-a_fpb);
	finds_per_bin <- cbind(dummy,finds_per_bin);
	}
if (length(interval_names) < ncol(finds_per_bin))
	interval_names <- colnames(finds_per_bin);

taxon_ranges <- data.frame(fa=as.numeric(taxon_ranges[,1]),la=as.numeric(taxon_ranges[,2]));
rownames(taxon_ranges) <- rownames(finds_per_bin);
stage_nos <- min(taxon_ranges):max(taxon_ranges);
ttl_bins <- ncol(finds_per_bin);

z_fpb <- ncol(finds_per_bin);
Nbt <- NFt <- NbL <- Nb <- Nt <- c();
n_bt <- n_bb <- n_ot <- c();
for (b in (1+a_fpb):(z_fpb-1))	{
	range_throughs <- ((1:ntaxa)[taxon_ranges$fa<b])[(1:ntaxa)[taxon_ranges$fa<b] %in% (1:ntaxa)[taxon_ranges$la>b]];
	Nbt <- c(Nbt,nrow(taxon_ranges[range_throughs,]));
	first_appears <- ((1:ntaxa)[taxon_ranges$fa==b])[(1:ntaxa)[taxon_ranges$fa==b] %in% (1:ntaxa)[taxon_ranges$la>b]];
	NFt <- c(NFt,nrow(taxon_ranges[first_appears,]));
	last_appears <- ((1:ntaxa)[taxon_ranges$fa<b])[(1:ntaxa)[taxon_ranges$fa<b] %in% (1:ntaxa)[taxon_ranges$la==b]];
	NbL <- c(NbL,nrow(taxon_ranges[last_appears,]));
	n_bt <- c(n_bt,sum(finds_per_bin[range_throughs,b]>0));					# range throughs sampled in bin
	n_bb <- c(n_bb,sum(finds_per_bin[range_throughs,b-1]>0));					# enterers sampled before bottom
	n_ot <- c(n_ot,sum(finds_per_bin[range_throughs,b+1]>0));	# exiters sampled after top
	}
Nb <- NbL+Nbt;
Nt <- NFt+Nbt;

output <- data.frame(n_bt=as.numeric(n_bt),Nbt=as.numeric(Nbt),Nb=as.numeric(Nb),NbL=as.numeric(NbL),n_ot=as.numeric(n_ot),Nt=as.numeric(Nt),NFt=as.numeric(NFt),n_bb=as.numeric(n_bb));
rownames(output) <- interval_names[(1+a_fpb):(z_fpb-1)];
return(output);
}

# get basic numbers for classic synoptic analyses
summarize_data_for_all_taxon_analysis <- function(finds_per_bin,taxon_ranges,interval_names="")	{
ntaxa <- nrow(taxon_ranges)
## make sure that the first column of finds_per_bin corresponds to bin 1 in taxon_ranges
a_fpb <- min((1:ncol(finds_per_bin))[colSums(finds_per_bin)>0]);
stage_1 <- min(taxon_ranges);
if (a_fpb < stage_1)	{
	dummy <- array("",dim=c(ntaxa,stage_1-a_fpb));
	colnames(dummy) <- rep("",stage_1-a_fpb);
	finds_per_bin <- cbind(dummy,finds_per_bin);
	}
if (length(interval_names) < ncol(finds_per_bin))
	interval_names <- colnames(finds_per_bin);

taxon_ranges <- data.frame(fa=as.numeric(taxon_ranges[,1]),la=as.numeric(taxon_ranges[,2]));
rownames(taxon_ranges) <- rownames(finds_per_bin);
stage_nos <- min(taxon_ranges):max(taxon_ranges);
ttl_bins <- ncol(finds_per_bin);

z_fpb <- ncol(finds_per_bin);
Si <- Ss <- Sp <- ni <- ns <- np <- uno <- c();
#n_bt <- n_bb <- n_ot <- c();
for (b in (1+a_fpb):(z_fpb-1))	{
	standing <- ((1:ntaxa)[taxon_ranges$fa<=b])[(1:ntaxa)[taxon_ranges$fa<=b] %in% (1:ntaxa)[taxon_ranges$la>=b]];
	surviving <- standing[taxon_ranges$la[standing]>b];
	preceding <- standing[taxon_ranges$fa[standing]<b];
	singletons <- ((1:ntaxa)[taxon_ranges$fa==b])[(1:ntaxa)[taxon_ranges$fa==b] %in% (1:ntaxa)[taxon_ranges$la==b]];
	sampled_during <- finds_per_bin[standing,b];
	sampled_after <- finds_per_bin[standing,b+1];
	sampled_before <- finds_per_bin[standing,b-1];
	Si <- c(Si,length(standing));
	Ss <- c(Ss,length(surviving));
	Sp <- c(Sp,length(preceding));
	ni <- c(ni,sum(sampled_during>0));
	ns <- c(ns,sum(sampled_after>0));
	np <- c(np,sum(sampled_before>0));
	uno <- c(uno,length(singletons));
	}

output <- data.frame(Si=as.numeric(Si),
	ni=as.numeric(ni),
	singletons=as.numeric(uno),
	Ss=as.numeric(Ss),
	ns=as.numeric(ns),
	Sp=as.numeric(Sp),
	np=as.numeric(np));
rownames(output) <- interval_names[(1+a_fpb):(z_fpb-1)];
return(output);
}