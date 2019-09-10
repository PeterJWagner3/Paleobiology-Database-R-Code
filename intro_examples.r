#install.packages("devtools")
### Test case.  Let's say we want to find which Ordovician graptolite
###		species have taxonomic authority data and which do not.  This will 
###		look at all graptolite records entered after the 19th of November 1998
###		(i.e., the entire history of the PaleoDB) and return three files.
###			1. A file of species with authority data plus all of the genus-
###				species combinations under which it has occurrences
###			2. A file of species without authority data, including the epithet
###				in a separate column (allowing you to sort all "smithi" together)
###				and the number of times that species occurs in the PaleoDB.
###			3. A file giving the sources for the species that do not have authority
###				data.  
taxa <- "Graptolithina"
onset <- "Ordovician"
end <- "Ordovician"
oldest_entry <- "1998-11-19"
# use today's date
latest_entry <- strsplit(as.character(Sys.time())," ")[[1]][1]
# If "TRUE" then the routine automatically outputs relevant files
save_files <- TRUE
# If you are saving the files, then you can specify how to append the files.
#	I like to append with ".xls" so that Excel will open it directly.
#	If you use ".cvs" then the output is comma-delimited rather than tab delimited
output_type <- ".xls"
	
data_breakdowns <- paleodb_entered_species_for_curation(taxa,onset,end,oldest_entry,latest_entry,save_files,output_type)

### Test case #2.  Let's say we want information about marine rock units from
###		the Ordovician.  This routine downloads all collections with animals,
###		plants or protists from the PaleoDB, and then separates out those with
###		marine environments.  It then returns for each formation-member combination
###			1. A file giving all of the intervals to which that rock unit is
###				assigned;
###			2. A file giving all of the zones to which that rock unit is assigned.
###			3. A file of rock units that are considered both formations and 
###				members by different workers, with the sources providing the
###				different rank opinions.  
end <- onset <- "Ordovician"
system <- "Marine"
delete_rock_type <- TRUE	# removes rock types from names for comparisons
save_files <- TRUE			# saves three files for stages, zones & inconsistent formations & members
output_type <- ".xls"		# 
strata_summaries <- paleodb_rock_units_for_curation(onset,end,taxa=c("Animalia","Protista","Plantae"),system,delete_rock_type,save_files,output_type)

end <- onset <- "Cenozoic"
system <- "Terrestrial"
taxa <- "Squamata"
delete_rock_type <- TRUE	# removes rock types from names for comparisons
save_files <- TRUE			# saves three files for stages, zones & inconsistent formations & members
output_type <- ".xls"		# 
strata_summaries <- paleodb_rock_units_for_curation(onset,end,taxa,system,delete_rock_type,save_files,output_type)
