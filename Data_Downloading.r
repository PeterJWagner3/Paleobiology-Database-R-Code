install.packages("paleobioDB", dependencies=TRUE)
library(paleobioDB)

# this simply gets all of the occurrence, taxonomic and locality data for a taxon.
get_a_group<-function(taxon,file_format)
{
httpO<-paste("http://paleobiodb.org/data1.1/occs/list.txt?base_name=",taxon,"&show=time,loc,crmod&limit=all",sep="")
occurrences <-read.table(httpO, sep=',', header=T)
#finds<-dim(occurrences)[1]
# get occurrence data: if this can be modified to make it species-only, then do that.
occurrences$taxon_name<-paleodb_species_occurrence_name_cleaning(occurrences$taxon_name)

paleodb_species_occurrence_name_cleaning(occurrences$taxon_name)
occ_species<-rbind(subset(occurrences, occurrences$taxon_rank=="species"))
colnames(occ_species)<-colnames(occurrences)
  
httpC<-paste("http://paleobiodb.org/data1.1/colls/list.txt?base_name=",taxon,"&show=time,stratext,geo,lithext,paleoloc,crmod&limit=all",sep="")
collections <-read.table(httpC, sep=',', header=T)
#locales<-dim(collections)[1]
  
httpT<-paste("http://paleobiodb.org/data1.1/taxa/list.txt?name=",taxon,"&rel=all_children&show=attr&limit=99999&rank=genus,subgenus,species",sep="")
taxonomy<-read.table(httpT, sep=',', header=T)
interval_names<-unique(collections$early_interval)
stages<-length(interval_names)
strat_names <- pbdb_intervals(vocab="pbdb",limit="all") #read all intervals at once
for (s in 1:stages) {
	i<-match(interval_names[s],strat_names$interval_name)
	if (s==1) {
		strat<-strat_names[i,]
    } else {
    strat<-rbind(strat,strat_names[i,])  
    }
  }
  #mid<-(early_late[,1]+early_late[,2])/2
filename1<-paste(taxon,"Records",sep="_")
filename1<-paste(filename1,file_format,sep="")
  
filename2<-paste(taxon,"Localities",sep="_")
filename2<-paste(filename2,file_format,sep="")
  
filename3<-paste(taxon,"Chronostratigraphic_Intervals",sep="_")
filename3<-paste(filename3,file_format,sep="")
  
sortdata<-c("early_age","late_age")
strat<-strat[do.call("order",strat[sortdata]),]
  
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
paleodb_data_download<-function(taxon,oldest,youngest,get_records,get_localities,get_abundances,get_taxonomy,file_format)
{
  
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
	time_range<-paste(oldest,"-",youngest,sep="")
  } else {
  time_range<-oldest
  }
  
if (get_records==TRUE)	{
	if (get_abundances=="TRUE") {
      httpO<-paste("http://paleobiodb.org/data1.1/occs/list.txt?base_name=",taxon,"&max_ma=",start_ma,"&min_ma=",end_ma,"&show=abund,time,loc,crmod&limit=all",sep="")
    } else {
      httpO<-paste("http://paleobiodb.org/data1.1/occs/list.txt?base_name=",taxon,"&max_ma=",start_ma,"&min_ma=",end_ma,"&show=time,loc,crmod&limit=all",sep="")  
    }
    occurrences <-read.table(httpO, sep=',', header=T)
    #finds<-dim(occurrences)[1]
    # get occurrence data: if this can be modified to make it species-only, then do that.
    occurrences$taxon_name<-paleodb_species_occurrence_name_cleaning(occurrences$taxon_name)
    filename<-paste(taxon,time_range,"Records",sep="_")
    filename<-paste(filename,file_format,sep="")
    if (file_format!=".csv")	{
      write.table(occurrences,filename,sep="\t",eol="\n",row.names=FALSE)
    } else if (file_format==".csv")	{
      write.table(occurrences,filename,sep=",",eol="\n",row.names=FALSE)
    }
  }	
  
  if (get_localities==TRUE)	{
    httpC<-paste("http://paleobiodb.org/data1.1/colls/list.txt?base_name=",taxon,"&max_ma=",start_ma,"&min_ma=",end_ma,"&show=time,stratext,geo,lithext,paleoloc,crmod&limit=all",sep="")
    collections <-read.table(httpC, sep=',', header=T)
    collections$formation<-clean_rock_unit_names(collections$formation,delete_rock_type=TRUE)
    collections$member<-clean_rock_unit_names(collections$member,delete_rock_type=FALSE)
    collections$late_age<-clear_na_from_vector(collections$late_age, 0)
    collections$paleomodel<-clear_na_from_vector(collections$paleomodel, "")
    collections$paleomodel2<-clear_na_from_vector(collections$paleomodel2, "")
    collections$paleomodel3<-clear_na_from_vector(collections$paleomodel3, "")
    collections$paleomodel4<-clear_na_from_vector(collections$paleomodel4, "")
    collections$paleolng2<-clear_na_from_vector(collections$paleolng2, "")
    collections$paleolng3<-clear_na_from_vector(collections$paleolng3, "")
    collections$paleolng4<-clear_na_from_vector(collections$paleolng4, "")
    collections$paleolat2<-clear_na_from_vector(collections$paleolat2, "")
    collections$paleolat3<-clear_na_from_vector(collections$paleolat3, "")
    collections$paleolat4<-clear_na_from_vector(collections$paleolat4, "")
    collections$geoplate2<-clear_na_from_vector(collections$geoplate2, "")
    collections$regionalsection<-clear_na_from_vector(collections$regionalsection, "")
    collections$regionalbed<-clear_na_from_vector(collections$regionalbed, "")
    collections$regionalorder<-clear_na_from_vector(collections$regionalorder, "")
    collections$collection_subset<-clear_na_from_vector(collections$collection_subset, "")
    filename<-paste(taxon,time_range,"Collections",sep="_")
    filename<-paste(filename,file_format,sep="")
    if (file_format!=".csv")	{
      write.table(collections,filename,sep="\t",eol="\n",row.names=FALSE)
    } else if (file_format==".csv")	{
      write.table(collections,filename,sep=",",eol="\n",row.names=FALSE)
    }
  }
  
  if (get_taxonomy==TRUE)	{
    httpT<-paste("http://paleobiodb.org/data1.1/taxa/list.txt?name=",taxon,"&rel=all_children&show=attr&limit=99999&rank=genus,subgenus,species",sep="")
    taxonomy<-read.table(httpT, sep=',', header=T)
    filename<-paste(taxon,time_range,"Taxonomy",sep="_")
    filename<-paste(filename,file_format,sep="")
    if (file_format!=".csv")	{
      write.table(taxonomy,filename,sep="\t",eol="\n",row.names=FALSE)
    } else if (file_format==".csv")	{
      write.table(taxonomy,filename,sep=",",eol="\n",row.names=FALSE)
    }
  }
}

# this gets occurrence, taxonomic and/or locality data for a taxon.  Enter 'TRUE' to get records, localities and/or taxonomy.  Note that it gets only occurrences with specimen or individual occurrences
paleodb_abundance_data_download<-function(taxon,oldest,youngest,get_localities,file_format)
{
	
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
	
	if (oldest!=youngest)	{
		time_range<-paste(oldest,"-",youngest,sep="")
	} else {
		time_range<-oldest
	}
	
	httpO<-paste("http://paleobiodb.org/data1.1/occs/list.txt?base_name=",taxon,"&max_ma=",start_ma,"&min_ma=",end_ma,"&show=abund,time,loc,crmod&limit=all",sep="")
	occurrences <-read.table(httpO, sep=',', header=T)
	abundances<-subset(occurrences,occurrences$abund_unit=="specimens")
	abundances<-rbind(abundances,subset(occurrences,occurrences$abund_unit=="individuals"))
	
	abundances$taxon_name<-paleodb_species_occurrence_name_cleaning(abundances$taxon_name)
	filename<-paste(taxon,time_range,"Abundance_Records",sep="_")
	filename<-paste(filename,file_format,sep="")
	if (file_format!=".csv")	{
		write.table(occurrences,filename,sep="\t",eol="\n",row.names=FALSE)
	} else if (file_format==".csv")	{
		write.table(occurrences,filename,sep=",",eol="\n",row.names=FALSE)
	}
	
if (get_localities==TRUE)  {
	localities<-unique(abundances$collection_no)
	coll<-length(localities)
	httpC<-paste("http://paleobiodb.org/data1.1/colls/list.txt?id=",localities[1],"&show=loc,time",sep="")
	collections<-read.table(httpC, sep=',', header=T)
	for (c in 2:coll) {
		httpC<-paste("http://paleobiodb.org/data1.1/colls/list.txt?id=",localities[c],"&show=loc,time",sep="")
	collections<-rbind(collections,read.table(httpC, sep=',', header=T))
	}
		
	collections$formation<-clean_rock_unit_names(collections$formation,delete_rock_type=TRUE)
	collections$member<-clean_rock_unit_names(collections$member,delete_rock_type=FALSE)
	collections$late_age<-clear_na_from_vector(collections$late_age, 0)
	collections$paleomodel<-clear_na_from_vector(collections$paleomodel, "")
	collections$paleomodel2<-clear_na_from_vector(collections$paleomodel2, "")
	collections$paleomodel3<-clear_na_from_vector(collections$paleomodel3, "")
	collections$paleomodel4<-clear_na_from_vector(collections$paleomodel4, "")
	collections$paleolng2<-clear_na_from_vector(collections$paleolng2, "")
	collections$paleolng3<-clear_na_from_vector(collections$paleolng3, "")
	collections$paleolng4<-clear_na_from_vector(collections$paleolng4, "")
	collections$paleolat2<-clear_na_from_vector(collections$paleolat2, "")
	collections$paleolat3<-clear_na_from_vector(collections$paleolat3, "")
	collections$paleolat4<-clear_na_from_vector(collections$paleolat4, "")
	collections$geoplate2<-clear_na_from_vector(collections$geoplate2, "")
	collections$regionalsection<-clear_na_from_vector(collections$regionalsection, "")
	collections$regionalbed<-clear_na_from_vector(collections$regionalbed, "")
	collections$regionalorder<-clear_na_from_vector(collections$regionalorder, "")
	collections$collection_subset<-clear_na_from_vector(collections$collection_subset, "")
		filename<-paste(taxon,time_range,"Collections",sep="_")
		filename<-paste(filename,file_format,sep="")
		if (file_format!=".csv")	{
			write.table(collections,filename,sep="\t",eol="\n",row.names=FALSE)
		} else if (file_format==".csv")	{
			write.table(collections,filename,sep=",",eol="\n",row.names=FALSE)
		}
	}
}

# this gets occurrence, taxonomic and/or locality data for a taxon.  Enter 'TRUE' to get records, localities and/or taxonomy. This retrieves occurrences only from one country.
paleodb_data_download_by_country<-function(taxon,oldest,youngest,get_records,get_localities,get_abundances,get_taxonomy,country,file_format)
{
	
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
	time_range<-paste(oldest,"-",youngest,sep="")
	} else {
	time_range<-oldest
	}
	
if (get_records==TRUE)	{
	if (get_abundances=="TRUE") {
		httpO<-paste("http://paleobiodb.org/data1.1/occs/list.txt?base_name=",taxon,"&max_ma=",start_ma,"&min_ma=",end_ma,"&show=abund,time,loc,crmod&limit=all",sep="")
		} else {
		httpO<-paste("http://paleobiodb.org/data1.1/occs/list.txt?base_name=",taxon,"&max_ma=",start_ma,"&min_ma=",end_ma,"&show=time,loc,crmod&limit=all",sep="")
		}
	oc<-read.table(httpO, sep=',', header=T)
	occurrences<-subset(oc,oc$cc==country)
		#finds<-dim(occurrences)[1]
		# get occurrence data: if this can be modified to make it species-only, then do that.
	occurrences$taxon_name<-paleodb_species_occurrence_name_cleaning(occurrences$taxon_name)
	filename<-paste(taxon,time_range,"Records_from",country,sep="_")
	filename<-paste(filename,file_format,sep="")
	if (file_format!=".csv")	{
		write.table(occurrences,filename,sep="\t",eol="\n",row.names=FALSE)
		} else if (file_format==".csv")	{
		write.table(occurrences,filename,sep=",",eol="\n",row.names=FALSE)
		}
	}	
	
if (get_localities==TRUE)	{
	httpC<-paste("http://paleobiodb.org/data1.1/colls/list.txt?base_name=",taxon,"&max_ma=",start_ma,"&min_ma=",end_ma,"&show=time,stratext,geo,lithext,loc,paleoloc,crmod&limit=all",sep="")
	coll<-read.table(httpC, sep=',', header=T)
	collections<-subset(coll, coll$cc==country)
	collections$formation<-clean_rock_unit_names(collections$formation,delete_rock_type=TRUE)
	collections$member<-clean_rock_unit_names(collections$member,delete_rock_type=FALSE)
	collections$late_interval<-clear_na_from_vector(collections$late_interval,"")
	collections$paleomodel<-clear_na_from_vector(collections$paleomodel, "")
	collections$paleomodel2<-clear_na_from_vector(collections$paleomodel2, "")
	collections$paleomodel3<-clear_na_from_vector(collections$paleomodel3, "")
	collections$paleomodel4<-clear_na_from_vector(collections$paleomodel4, "")
	collections$paleolng2<-clear_na_from_vector(collections$paleolng2, "")
	collections$paleolng3<-clear_na_from_vector(collections$paleolng3, "")
	collections$paleolng4<-clear_na_from_vector(collections$paleolng4, "")
	collections$paleolat2<-clear_na_from_vector(collections$paleolat2, "")
	collections$paleolat3<-clear_na_from_vector(collections$paleolat3, "")
	collections$paleolat4<-clear_na_from_vector(collections$paleolat4, "")
	collections$geoplate2<-clear_na_from_vector(collections$geoplate2, "")
	collections$regionalsection<-clear_na_from_vector(collections$regionalsection, "")
	collections$regionalbed<-clear_na_from_vector(collections$regionalbed, "")
	collections$regionalorder<-clear_na_from_vector(collections$regionalorder, "")
	collections$collection_subset<-clear_na_from_vector(collections$collection_subset, "")
	filename<-paste(taxon,time_range,"Collections_from",country,sep="_")
	filename<-paste(filename,file_format,sep="")
	if (file_format!=".csv")	{
		write.table(collections,filename,sep="\t",eol="\n",row.names=FALSE)
		} else if (file_format==".csv")	{
		write.table(collections,filename,sep=",",eol="\n",row.names=FALSE)
		}
	}
	
if (get_taxonomy==TRUE)	{
	httpT<-paste("http://paleobiodb.org/data1.1/taxa/list.txt?name=",taxon,"&rel=all_children&show=attr&limit=99999&rank=genus,subgenus,species",sep="")
	taxonomy<-read.table(httpT, sep=',', header=T)
	filename<-paste(taxon,time_range,"Taxonomy",sep="_")
	filename<-paste(filename,file_format,sep="")
	if (file_format!=".csv")	{
		write.table(taxonomy,filename,sep="\t",eol="\n",row.names=FALSE)
		} else if (file_format==".csv")	{
		write.table(taxonomy,filename,sep=",",eol="\n",row.names=FALSE)
		}
	}
}


# this gets occurrence, taxonomic and/or locality data for a taxon.  Enter 'TRUE' to get records, localities and/or taxonomy. This retrieves occurrences only from one country.
paleodb_data_download_by_country<-function(taxon,oldest,youngest,get_records,get_localities,get_abundances,get_taxonomy,environment,file_format)
{

if (environment=="terrestrial") {
  environs<-c("terrestrial indet.", "fluvial indet.", "alluvial fan", "channel lag", "coarse channel fill", "fine channel fill", "channel", "wet floodplain", "dry floodplain", "&quot;floodplain&quot;", "crevasse splay", "levee", "mire/swamp", "fluvial-lacustrine indet.", "delta plain", "fluvial-deltaic indet.", "lacustrine - large", "lacustrine - small", "pond", "crater lake", "lacustrine delta plain", "lacustrine interdistributary bay", "lacustrine delta front", "lacustrine prodelta", "lacustrine deltaic indet.", "lacustrine indet.", "dune", "interdune", "loess", "eolian indet.", "cave", "fissure fill", "sinkhole", "karst indet.", "tar", "mire/swamp", "spring", "glacial")
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
  time_range<-paste(oldest,"-",youngest,sep="")
  } else {
  time_range<-oldest
  }

if (get_records==TRUE)	{
  if (get_abundances=="TRUE") {
    httpO<-paste("http://paleobiodb.org/data1.1/occs/list.txt?base_name=",taxon,"&max_ma=",start_ma,"&min_ma=",end_ma,"&show=abund,time,loc,crmod&limit=all",sep="")
    } else {
    httpO<-paste("http://paleobiodb.org/data1.1/occs/list.txt?base_name=",taxon,"&max_ma=",start_ma,"&min_ma=",end_ma,"&show=time,loc,geo,crmod&limit=all",sep="")
    }
  oc<-read.table(httpO, sep=',', header=T)
  oc$environment<-clear_na_from_vector(oc$environment, "")
  occurrences<-subset(oc,oc$environment==environs)
  occurrences$reid_no<-clear_na_from_vector(occurrences$reid_no,"")
  occurrences$superceded<-clear_na_from_vector(occurrences$superceded,"")
  #finds<-dim(occurrences)[1]
  # get occurrence data: if this can be modified to make it species-only, then do that.
  occurrences$taxon_name<-paleodb_species_occurrence_name_cleaning(occurrences$taxon_name)
  filename<-paste(taxon,time_range,"Records_from",environment,sep="_")
  filename<-paste(filename,file_format,sep="")
  if (file_format!=".csv")	{
    write.table(occurrences,filename,sep="\t",eol="\n",row.names=FALSE)
    } else if (file_format==".csv")	{
    write.table(occurrences,filename,sep=",",eol="\n",row.names=FALSE)
    }
  }


if (get_localities==TRUE)	{
  httpC<-paste("http://paleobiodb.org/data1.1/colls/list.txt?base_name=",taxon,"&max_ma=",start_ma,"&min_ma=",end_ma,"&show=time,stratext,geo,lithext,loc,paleoloc,crmod&limit=all",sep="")
  coll<-read.table(httpC, sep=',', header=T)
  collections<-subset(coll, coll$cc==country)
  collections$formation<-clean_rock_unit_names(collections$formation,delete_rock_type=TRUE)
  collections$member<-clean_rock_unit_names(collections$member,delete_rock_type=FALSE)
  collections$late_interval<-clear_na_from_vector(collections$late_interval,"")
  collections$paleomodel<-clear_na_from_vector(collections$paleomodel, "")
  collections$paleomodel2<-clear_na_from_vector(collections$paleomodel2, "")
  collections$paleomodel3<-clear_na_from_vector(collections$paleomodel3, "")
  collections$paleomodel4<-clear_na_from_vector(collections$paleomodel4, "")
  collections$paleolng2<-clear_na_from_vector(collections$paleolng2, "")
  collections$paleolng3<-clear_na_from_vector(collections$paleolng3, "")
  collections$paleolng4<-clear_na_from_vector(collections$paleolng4, "")
  collections$paleolat2<-clear_na_from_vector(collections$paleolat2, "")
  collections$paleolat3<-clear_na_from_vector(collections$paleolat3, "")
  collections$paleolat4<-clear_na_from_vector(collections$paleolat4, "")
  collections$geoplate2<-clear_na_from_vector(collections$geoplate2, "")
  collections$regionalsection<-clear_na_from_vector(collections$regionalsection, "")
  collections$regionalbed<-clear_na_from_vector(collections$regionalbed, "")
  collections$regionalorder<-clear_na_from_vector(collections$regionalorder, "")
  collections$collection_subset<-clear_na_from_vector(collections$collection_subset, "")
  filename<-paste(taxon,time_range,"Collections_from",country,sep="_")
  filename<-paste(filename,file_format,sep="")
  if (file_format!=".csv")	{
    write.table(collections,filename,sep="\t",eol="\n",row.names=FALSE)
    } else if (file_format==".csv")	{
    write.table(collections,filename,sep=",",eol="\n",row.names=FALSE)
    }
	}

if (get_taxonomy==TRUE)	{
  httpT<-paste("http://paleobiodb.org/data1.1/taxa/list.txt?name=",taxon,"&rel=all_children&show=attr&limit=99999&rank=genus,subgenus,species",sep="")
  taxonomy<-read.table(httpT, sep=',', header=T)
  filename<-paste(taxon,time_range,"Taxonomy",sep="_")
  filename<-paste(filename,file_format,sep="")
  if (file_format!=".csv")	{
    write.table(taxonomy,filename,sep="\t",eol="\n",row.names=FALSE)
	  } else if (file_format==".csv")	{
    write.table(taxonomy,filename,sep=",",eol="\n",row.names=FALSE)
	  }
	}
#occs.marine <- occs.raw %>% filter(cx_int_no %in% stage$interval_no) %>%filter(environment %in% marine_env)
}


taxon<-"Gastropoda"
oldest<-"Ordovician"
youngest<-"Ordovician"
file_format<-".xls"
taxon_level<-"species"
# this gets occurrence, taxonomic and/or locality data for a taxon.  Enter 'TRUE' to get records, localities and/or taxonomy.
paleodb_taxonomic_data_download<-function(taxon,oldest,youngest,taxon_level,file_format)
{
  
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
  
if (oldest!=youngest)  {
  time_range<-paste(oldest,"-",youngest,sep="")
  } else {
  time_range<-oldest
  }
  
httpO<-paste("http://paleobiodb.org/data1.1/occs/list.txt?base_name=",taxon,"&max_ma=",start_ma,"&min_ma=",end_ma,"&show=time,loc,crmod&limit=all",sep="")  
occurrences <-read.table(httpO, sep=',', header=T)
    #finds<-dim(occurrences)[1]
    # get occurrence data: if this can be modified to make it species-only, then do that.
occurrences$taxon_name<-paleodb_species_occurrence_name_cleaning(occurrences$taxon_name)
filename<-paste(taxon,time_range,"Records",sep="_")
filename<-paste(filename,file_format,sep="")
if (file_format!=".csv")	{
  write.table(occurrences,filename,sep="\t",eol="\n",row.names=FALSE)
  } else if (file_format==".csv")	{
  write.table(occurrences,filename,sep=",",eol="\n",row.names=FALSE)
  }

httpT<-paste("http://paleobiodb.org/data1.1/taxa/list.txt?name=",taxon,"&rel=all_children&show=attr&limit=99999&rank=",taxon_level,sep="")
taxonomy<-read.table(httpT, sep=',', header=T)
taxonomy$associated_records<-clear_na_from_vector(taxonomy$associated_records, "")
taxonomy$common_name<-clear_na_from_vector(taxonomy$common_name, "")
taxonomy$attribution<-clear_na_from_vector(taxonomy$attribution, "")
taxonomy$pubyr<-clear_na_from_vector(taxonomy$pubyr, "")
spc<-dim(taxonomy)[1]

### START HERE!!!
for (s in 1:spc)  {
  if (s==1) { opinions<-pbdb_ref_taxa(id=taxonomy$orig_no[s],vocab="pbdb")
    } else {
    ops<-pbdb_ref_taxa(id=taxonomy$orig_no[s],vocab="pbdb")  
    opinions<-rbind(opinions,ops)
    }
  }
#http://testpaleodb.geology.wisc.edu/data1.2/taxa/opinions.txt?base_id=58103&select=all
filename<-paste(taxon,time_range,"Taxonomy",sep="_")
filename<-paste(filename,file_format,sep="")
if (file_format!=".csv")	{
  write.table(taxonomy,filename,sep="\t",eol="\n",row.names=FALSE)
  } else if (file_format==".csv")	{
  write.table(taxonomy,filename,sep=",",eol="\n",row.names=FALSE)
  }
}

get_references_for_collections<-function(collections,project_name)
{
N<-length(collections)
#http://paleobiodb.org/data1.2/occs/list.txt?coll_id=50068&show=ref
#http://paleobiodb.org/data1.2/colls/list.txt?id=50068&show=ref
for (c in 1:N)  {
    coll<-collections[c]
    httpO<-paste("http://paleobiodb.org/data1.2/occs/list.txt?coll_id=",coll,"&show=ref",sep="")
    finds<-read.table(httpO, sep=',', header=T)
    new_refs<-unique(finds$reference_no)
    if (c==1)   {
        references<-new_refs
        }   else {
        references<-c(references,new_refs)
        }
    }
#ttlrefs<-max(references)
#impact<-vector(length=ttlrefs)
#for (t in 1:ttlrefs)    impact[references[t]]<-impact[references[t]]+1

novel_refs<-sort(unique(subset(references,references<100000)))
refs<-length(novel_refs)

for (r in 1:refs)   {
    cit<-novel_refs[r]
    httpR<-paste("http://paleobiodb.org/data1.2/refs/single.ris?id=",cit,"&textresult",sep="")
    article<-read.table(httpR, sep="\n", header=T)
    ll<-dim(article)[1]
    art<-article[3:ll,1]
    art<-c(as.character(art),"")
    if (r==1)   {
        alexandria<-art
        }   else {
        alexandria<-c(as.character(alexandria),as.character(art))
        }
    }
filename<-paste(project_name,"References.ris",sep="_")
write.table(alexandria,filename,sep="",quote=FALSE, eol="\n",row.names=FALSE,col.names=FALSE)
return(alexandria)
}
