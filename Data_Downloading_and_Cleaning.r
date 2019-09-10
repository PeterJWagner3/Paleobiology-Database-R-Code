#### Setup Program ####
common_source_folder <- "~/Documents/R_Projects/Common_R_Source_Files/";		# directory to folder where you keep common source
data_for_R_folder <- "~/Documents/R_Projects/Data_for_R/";						# directory to folder where you keep R-Data files
local_directory <- ("~/Documents/R_Projects/Diversification_for_RevBayes/");	# directory from which you are operating this program
source('~/Documents/R_Projects/Common_R_Source_Files/Chronos.r')  #
source('~/Documents/R_Projects/Common_R_Source_Files/Data_Downloading_v4.R')  #
source('~/Documents/R_Projects/Common_R_Source_Files/Estimating_Abundance_Distributions.r');
source('~/Documents/R_Projects/Common_R_Source_Files/General_Plot_Templates.R')  #
source('~/Documents/R_Projects/Common_R_Source_Files/Historical_Diversity_Metrics.r')  #
source('~/Documents/R_Projects/Common_R_Source_Files/Occurrence_Data_Routines.r')  #
source('~/Documents/R_Projects/Common_R_Source_Files/Sampling_and_Occupancy_Distributions.r')  #
source('~/Documents/R_Projects/Common_R_Source_Files/Wagner_Kluges.r')  #
source('~/Documents/R_Projects/Common_R_Source_Files/Wagner_Stats_and_Probability_101.r');  #

load(paste(data_for_R_folder,"Rock_Unit_Database.RData",sep=""));  		# data for rock-units including biozanation & superposition
load(paste(data_for_R_folder,"Gradstein_2012_Augmented.RData",sep="")); # refined Gradstein 2012 timescale & biozonations
load(paste(data_for_R_folder,"PaleoDB_Edits.RData",sep=""));			# edits to collections I cannot edit. # refined Gradstein 2012 timescale & biozonations

#### Load External Databases ####
rock_database <- rock_unit_data$rock_unit_database;				# Basic information about Rock Units
rock_to_zone_database <- rock_unit_data$rock_to_zone_database;	# Rock to Zone data	
time_scale <- gradstein_2012_emended$time_scale;				# Gradstein's 2012 scale plus a lot of emendations	
zone_database <- gradstein_2012_emended$zones;					# Biozonatin time scale

### Wagner's personal PaleoDB edits....
fossilworks_collections <- paleodb_fixes$fossilworks_collections;
paleodb_rock_reidentifications <- paleodb_fixes$paleodb_rock_reidentifications;
paleodb_collection_edits <- paleodb_fixes$paleodb_collection_edits;
if (!is.na(match("X",colnames(paleodb_collection_edits))))	{
	paleodb_collection_edits[,match("X",colnames(paleodb_collection_edits)):ncol(paleodb_collection_edits)] <- NULL;
	}

#### Tell R how to get your data! ####
analysis_name <- "Metazoans"					# This will be used to name files
species_only <- T;								# if T, the Bellerophon sp. occurrences are culled.
basic_environments <- c("marine","unknown");	# use "terr" for terrestrial "marine" for marine & "unknown" for unspecificed
onset <- "Cambrian";							# oldest collections & occurrences
end <- "Silurian"								# youngest collections & occurrences
taxon <- c("Trilobita","Porifera","Rugosa","Tabulata","Mollusca","Echinodermata","Brachiopoda","Conodonta","Graptolithina");
save_files <- F;								# This will save the data for each taxonomic group as csv files if TRUE; leave it F for now. 
lump_intervals <- F;							# if you want to lump intervals, then enter T; you will need to provide a key

time_scale_stratigraphic_scale <- "Stage Slice"; # This must match a value in time_scale$scale
temporal_precision <- 0.05;
library(rlist);		#install.packages("rlist", dependencies=TRUE);

#### Download Data ####
# The data will be collected in two ways.  
all_finds <- all_sites <- list();	# this will make a list for collections & occurrences for each taxonomic group
tx <- 0;
while (tx <  length(taxon))	{
	tx <- tx + 1;
	print(paste("Doing",taxon[tx]));
	taxon_finds <- accio_occurrence_data(taxa=taxon[tx],onset,end,basic_environments,species_only,clean_entered_taxa=TRUE,
					directory=local_directory,save_files,output_type=".csv");
	taxon_sites <- accio_collection_data(taxa=taxon[tx],onset,end,basic_environments,standardize_members=F,
					directory=local_directory,save_files,species_only=F,output_type=".csv");
	taxon_sites <- taxon_sites[taxon_sites$collection_no %in% taxon_finds$collection_no,];

	taxon_finds <- add_subgenus_names_to_paleodb_finds(paleodb_finds = taxon_finds);
	taxon_finds <- subset(taxon_finds,!is.na(taxon_finds$subgenus));
	taxon_finds <- subset(taxon_finds,!is.na(taxon_finds$genus));
	if (taxon[tx]=="Rugosa" || taxon[tx]=="Tabulata")
		taxon_finds$class <- taxon[tx];
	if (tx==1)	{
		lumped_finds <- taxon_finds;
		lumped_sites <- taxon_sites;
		} else	{
		lumped_finds <- rbind(lumped_finds,taxon_finds);
		lumped_sites <- rbind(lumped_sites,taxon_sites[!taxon_sites$collection_no %in% lumped_sites$collection_no,]);
		}
	
	all_finds <- rlist::list.append(all_finds,taxon_finds);
	all_sites <- rlist::list.append(all_sites,taxon_sites);
	}
names(all_finds) <- names(all_sites) <- taxon;

### output data for individual groups
if (length(taxon)>1)	{
	for (af in 1:length(all_finds))	{
		filename <- paste(local_directory,names(all_finds)[af],"_Finds.csv",sep="");
		write.csv(all_finds[[af]],file=filename,row.names = F);
		}
	}

#### Prepare data for Fine Tuning with External Databases ####
lumped_finds_init <- lumped_finds;	# for debugging...
lumped_finds <- lumped_finds_init[order(lumped_finds_init$collection_no,lumped_finds_init$collection_no),];

lumped_sites_init <- lumped_sites;
lumped_sites <- lumped_sites_init[order(lumped_sites_init$collection_no),];

lumped_sites <- reparo_paleodb_paleogeography_with_fossilworks_data(paleodb_collections=lumped_sites,fossil_works_geography=fossilworks_collections);
direct_dates <- data.frame(direct_ma=as.numeric(fossilworks_collections$direct_ma[match(lumped_sites$collection_no,fossilworks_collections$collection_no)]),
						   direct_ma_error=as.numeric(fossilworks_collections$direct_ma_error[match(lumped_sites$collection_no,fossilworks_collections$collection_no)]),
						   direct_ma_method=as.character(fossilworks_collections$direct_ma_method[match(lumped_sites$collection_no,fossilworks_collections$collection_no)]),
						   stringsAsFactors=hell_no);
direct_dates <- evanesco_na_from_matrix(direct_dates,"");
lumped_sites <- cbind(lumped_sites,direct_dates);

# this provides edits for paleodb collections that cannot currently be edited.
lumped_sites <- reparo_unedittable_paleodb_rock_identification(paleodb_collections=lumped_sites,paleodb_rock_reidentifications=paleodb_rock_reidentifications);
lumped_sites <- reparo_unedittable_paleodb_collections(paleodb_collections=lumped_sites,paleodb_collection_edits=paleodb_collection_edits);

# deal with time scale
hierarchical_chronostrat <- accio_hierarchical_timescale(chronostrat_units=unique(c(unique(as.character(lumped_sites$early_interval)),unique(as.character(lumped_sites$late_interval)))),time_scale=gradstein_2012_emended$time_scale,regional_scale=time_scale_stratigraphic_scale);
hierarchical_chronostrat$ma_lb <- temporal_precision*round(hierarchical_chronostrat$ma_lb/temporal_precision,0);
hierarchical_chronostrat$ma_ub <- temporal_precision*round(hierarchical_chronostrat$ma_ub/temporal_precision,0);
finest_chronostrat <- hierarchical_chronostrat[hierarchical_chronostrat$bin_first==hierarchical_chronostrat$bin_last,];
finest_chronostrat <- finest_chronostrat[match(finest_chronostrat$bin_first,finest_chronostrat$bin_first),];

print(paste("Redoing PaleoDB intervals to ",time_scale_stratigraphic_scale,"...",sep=""));
reset_sites <- reset_paleodb_intervals_to_desired_time_scale(collections=lumped_sites,finest_chronostrat = finest_chronostrat,time_scale);
reset_sites$min_ma <- round(temporal_precision*round(reset_sites$min_ma/temporal_precision,0),floor(-log10(temporal_precision)));
reset_sites$max_ma <- round(temporal_precision*round(reset_sites$max_ma/temporal_precision,0),floor(-log10(temporal_precision)));
ncolls <- nrow(reset_sites);

##### Use External Database to refine collections based on Rock Unit ####
print("Refining PaleoDB data with rock-unit and biozonation databases...");
#redone_collections <- match_paleodb_collections_to_external_stratigraphic_database(collections=reset_sites,wagner_rocks=rock_database);
zone_database <- gradstein_2012_emended$zones;
relv_slices <- unique(c(unique(reset_sites$early_interval),unique(reset_sites$late_interval)));
young_enough <- (1:nrow(zone_database))[zone_database$interval_lb %in% relv_slices];
old_enough <- (1:nrow(zone_database))[zone_database$interval_ub %in% relv_slices];
just_right <- young_enough[young_enough %in% old_enough];

refined_data <- refine_collection_dates_with_external_database(study=analysis_name,collections=reset_sites,rock_database,zone_database=zone_database[just_right,],rock_to_zone_database,time_scale=finest_chronostrat,directory=local_directory);
refined_sites <- refined_data$Recalibrated_Collections;

# now, rebin the data with the revised time scale 
age <- temporal_precision*round(refined_sites$ma_lb/temporal_precision,0);
refined_sites$interval_lb <- sapply(age,rebin_collection_with_time_scale,onset_or_end = "onset",fine_time_scale = finest_chronostrat);
age <- temporal_precision*round(refined_sites$ma_ub/temporal_precision,0);
refined_sites$interval_ub <- sapply(age,rebin_collection_with_time_scale,onset_or_end = "end",fine_time_scale = finest_chronostrat);

lumped_finds$ma_lb <- refined_sites$ma_lb[match(lumped_finds$collection_no,lumped_sites$collection_no)];
lumped_finds$ma_ub <- refined_sites$ma_ub[match(lumped_finds$collection_no,lumped_sites$collection_no)];

ncoll <- nrow(refined_sites);
bin_lb <- match(refined_sites$interval_lb,finest_chronostrat$interval);
bin_ub <- match(refined_sites$interval_ub,finest_chronostrat$interval);
# use quantitative biostratigraphy 101 to refine dates here.
print("Using basic biostratigraphy to minimize gaps for uncertainly aged collections...");
#unique(zone_database$interval_lb)
optim_sites_orig <- optim_sites <- optimo_paleodb_collection_and_occurrence_stratigraphy(paleodb_finds=lumped_finds,paleodb_collections=refined_sites,hierarchical_chronostrat=finest_chronostrat,zone_database=zone_database[just_right,],update_search=T);
file_starter <- paste("Optimized",onset,end,analysis_name,sep="_");
write.csv(optim_sites_orig,paste(local_directory,file_starter,"_Finest_Scale.csv",sep=""),row.names = F);

if (lump_intervals)	{
	interval_thesaurus_file <- file.choose();
	interval_thesaurus <- read.csv(interval_thesaurus_file,header=T,stringsAsFactors=hell_no);
	interval_thesaurus$bin_old <- as.numeric(1:nrow(interval_thesaurus));
	interval_thesaurus$bin <- match(interval_thesaurus$new_interval,unique(interval_thesaurus$new_interval));
	finest_chronostrat <- accio_coarsened_time_scale(hierarchical_chronostrat,interval_thesaurus);
	finest_chronostrat <- subset(finest_chronostrat,!is.na(finest_chronostrat$interval))
	optim_sites$interval_lb <- interval_thesaurus$new_interval[match(optim_sites$interval_lb,interval_thesaurus$orig_interval)];
	optim_sites$interval_ub <- interval_thesaurus$new_interval[match(optim_sites$interval_ub,interval_thesaurus$orig_interval)];
#	write.csv(optim_sites,paste(local_directory,"Optimized_C_S_Metazoans_lumped_slices.csv",sep=""),row.names = F);
	write.csv(optim_sites,paste(local_directory,file_starter,"_Lumped_Intervals.csv",sep=""),row.names = F);
	}
#youngest_zone_age <- min(abs(time_scale$ma_ub[match(interval_thesaurus$orig_interval,time_scale$interval)]));

optim_sites <- name_unnamed_rock_units(paleodb_collections=optim_sites,finest_chronostrat);
#write.csv(optim_sites,paste(local_directory,"Optimized_C_S_Metazoans_Added_rocks.csv",sep=""),row.names = F);
write.csv(optim_sites,paste(local_directory,file_starter,"_Added_Rocks.csv",sep=""),row.names = F);

ncolls <- nrow(optim_sites);
finest_chronostrat$ma_lb <- temporal_precision*round(finest_chronostrat$ma_lb/temporal_precision,0)
finest_chronostrat$ma_ub <- temporal_precision*round(finest_chronostrat$ma_ub/temporal_precision,0)
paleodb_collections <- completely_rebin_collections_with_uniform_time_scale(collections=optim_sites,uniform_time_scale = finest_chronostrat);
problematic_collections <- (1:nrow(paleodb_collections))[paleodb_collections$bin_ub<paleodb_collections$bin_lb];
paleodb_collections <- paleodb_collections[!(1:nrow(paleodb_collections)) %in% problematic_collections,];
write.csv(paleodb_collections,paste(local_directory,file_starter,"_Collections_Analyzed.csv",sep=""),row.names = F);
#write.csv(paleodb_collections,paste(local_directory,"C_S_Metazoan_Collections_Analyzed.csv",sep=""),row.names = F);
lumped_finds$phylum[lumped_finds$class %in% c("Lingulata","Micrina","Paterinata","Setatella")] <- "Linguliformea";
write.csv(lumped_finds,paste(local_directory,file_starter,"_Finds_Analyzed.csv",sep=""),row.names = F);

#### Get Basic Data Summaries ####
#paste(paleodb_collections$collection_no[problematic_collections],collapse = ",");
#cbind(paleodb_collections$interval_lb[problematic_collections],paleodb_collections$interval_ub[problematic_collections]);
#paleodb_collections$interval_ub[problematic_collections] <- paleodb_collections$interval_lb[problematic_collections];
#paleodb_collections$bin_ub[problematic_collections] <- paleodb_collections$bin_lb[problematic_collections];
#paleodb_collections$ma_ub[problematic_collections] <- 456.5;
taxa_and_controls <- read.csv(paste(local_directory,"Taxon_Controls.csv",sep=""),header=T,stringsAsFactors = hell_no);
higher_taxa <- unique(taxa_and_controls$Taxon);
for (cl in 1:length(higher_taxa))	{
#	clade_finds <- all_finds[[cl]];
	clade_finds <- rbind(subset(lumped_finds,lumped_finds$class==higher_taxa[cl]),subset(lumped_finds,lumped_finds$phylum==higher_taxa[cl]));
	clade_finds <- clade_finds[clade_finds$collection_no %in% paleodb_collections$collection_no,];
	control_taxa <- taxa_and_controls$Control[taxa_and_controls$Taxon==higher_taxa[cl]];
	control_taxa <- control_taxa[control_taxa!=higher_taxa[cl]];
	if (length(control_taxa)>0)	{
		control_finds <- c();
		for (cn in 1:length(control_taxa))
			control_finds <- rbind(control_finds,rbind(subset(lumped_finds,lumped_finds$class==control_taxa[cn]),subset(lumped_finds,lumped_finds$phylum==control_taxa[cn])));
		relv_colls <- sort(unique(c(match(clade_finds$collection_no,paleodb_collections$collection_no),match(control_finds$collection_no,paleodb_collections$collection_no))));
		} else	{
		relv_colls <- sort(unique(match(clade_finds$collection_no,paleodb_collections$collection_no)));
		}
	
#	clade_collections <- paleodb_collections[paleodb_collections$collection_no %in% clade_finds$collection_no,];
	clade_collections <- paleodb_collections[relv_colls,];
	if (cl==1)	{
		slice_rocks_clade <- tally_rock_units_occupied_by_subinterval(taxon_collections=clade_collections,finest_chronostrat);
		slice_sites_clade <- tally_collections_occupied_by_subinterval(taxon_collections=clade_collections,finest_chronostrat);
		} else	{
		slice_rocks_clade <- rbind(slice_rocks_clade,tally_rock_units_occupied_by_subinterval(taxon_collections=clade_collections,finest_chronostrat));
		slice_sites_clade <- rbind(slice_sites_clade,tally_collections_occupied_by_subinterval(taxon_collections=clade_collections,finest_chronostrat));
		}
	}
rownames(slice_sites_clade) <- rownames(slice_rocks_clade) <- higher_taxa;
write.csv(slice_sites_clade,paste(local_directory,"Sites_per_Clade_per_Slice.csv",sep=""),row.names=T);
write.csv(slice_rocks_clade,paste(local_directory,"Rocks_per_Clade_per_Slice.csv",sep=""),row.names=T);
