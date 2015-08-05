#install.packages("paleobioDB", dependencies=TRUE)
library(paleobioDB)
# data curation scripts

# taxonomy curation scripts
# Arguments:
# 	taxon: proper taxonomic name
# 	oldest: oldest geological interval from which you want new records
# 	youngest: youngest geological interval from which you want new records
# 	start_date: the oldest entered records that you want to include
# 	end_date: the most recent records that you want to include
# 	file_format: the end tag on the output files: '.xls' for Excel, '.txt', '.csv" for comma-delimited
paleodb_occurrence_curation<-function(taxon,oldest,youngest,start_date,end_date,file_format)
{

# get temporal information for search
strat_names <- pbdb_intervals(limit="all") #read all intervals at once
#sn<-read.table("http://paleobiodb.org/data1.1/intervals/list.txt?scale=1", sep=',', header=T)

if(oldest %in% strat_names$nam) {
	start_ma <- strat_names$eag[match(oldest, strat_names$nam)]
} else 
{ print("Error! ",oldest," is not a term in our chronostratigraphic lexicon!")}

if(youngest %in% strat_names$nam) {
	end_ma <- strat_names$lag[match(youngest, strat_names$nam)]
} else 
{ print("Error! ",youngest," is not a term in our chronostratigraphic lexicon!")}

http<-paste("http://paleobiodb.org/data1.1/occs/list.txt?base_name=",taxon,"&max_ma=",start_ma,"&min_ma=",end_ma,"&show=time,loc,crmod&limit=all&modified_since=",start_date,"&modified_before=",end_date,sep="")
occurrences <-read.table(http, sep=',', header=T)
finds<-dim(occurrences)[1]
# get occurrence data: if this can be modified to make it species-only, then do that.
occurrences$taxon_name<-paleodb_species_occurrence_name_cleaning(occurrences$taxon_name)

httpT<-paste("http://paleobiodb.org/data1.1/occs/taxa.txt?base_name=",taxon,"&max_ma=",start_ma,"&min_ma=",end_ma,"&show=attr&limit=all&modified_since=",start_date,"&modified_before=",end_date,sep="")
occurring_taxa<-read.table(httpT, sep=',', header=T)
#write.table(occurring_taxa,"Test.xls",sep="\t",eol="\n",row.names=TRUE)
	
dd<-dim(occurring_taxa)
occurring_taxa_no<-dd[1]
# 0: not a species record of any sort
# 1: species record with taxonomic data
# 2: species record without taxonomic data
entry_level<-vector(length=occurring_taxa_no)
for (i in 1:dd[1])  {
	if (as.character(occurring_taxa$rank[i])=="species" || as.character(occurring_taxa$rank[i])=="subspecies")	{
		entry_level[i]<-1
		}	else if (as.character(occurring_taxa$rank[i])=="genus" || as.character(occurring_taxa$rank[i])!="subgenus")	{
		words<-length(strsplit(as.character(occurring_taxa$taxon_name[i])," ")[[1]])
	last_name<-strsplit(as.character(occurring_taxa$taxon_name[i])," ")[[1]][words]
	if (last_name=="sp.")	{
			entry_level[i]<-0
			} else if (last_name=="indet.")	{
			entry_level[i]<-0
			} else if (last_name=="spp.")	{
			entry_level[i]<-0
			} else if (last_name=="nov.")	{
			entry_level[i]<-0
			} else if (last_name=="columnals")	{
			entry_level[i]<-0
			} else if (last_name=="holdfasts")	{
			entry_level[i]<-0
			} else if (words<2)	{
			entry_level[i]<-0	
			}	else {
			entry_level[i]<-2
			}	# end possibilities for a genus or subgenus rank entry
		}	# end case of a genus or subgenus rank of taxon
	}	# end search of occurring taxa

# set up matrices of species with and without taxonomic data
entered_species<-rbind(subset(occurring_taxa, entry_level==1))
unentered_species<-rbind(subset(occurring_taxa, entry_level==2))

entered_species$attribution<-clear_na_from_vector(entered_species$attribution,"")
# find the author names & latest opinions for entered species
entered<-dim(entered_species)[1]
latest_opinion<-vector(length=entered)
for (sp in 1:entered)	{
	opinions<-pbdb_ref_taxa(id=entered_species$taxon_no[sp],vocab="pbdb")
  # For some reason, the routine below introduces strange errors for those species that do not get an attribution.
#  words<-length(strsplit(as.character(entered_species$taxon_name[sp])," ")[[1]])
#  species_whole<-species_name<-strsplit(as.character(entered_species$taxon_name[sp])," ")[[1]][1]
#  for (n in 2:words)  {
#    a<-species_name<-strsplit(as.character(entered_species$taxon_name[sp])," ")[[1]][n]
#    species_whole<-paste(species_whole,a,sep="%20")
#    }
#  httpTR<-paste("http://paleobiodb.org/data1.1/taxa/refs.txt?base_name=",species_whole,"&textresult",sep="")
#  opinions<-read.table(httpTR, sep=',', header=T)  #### HERE!
	cites<-dim(opinions)[1]
	if (cites>0) {
		if (entered_species$attribution[sp]=="")	{
      b<-paste(as.character(opinions$author1last[1]),as.character(opinions$pubyr[1]),sep=" ")
      c<-paste(as.character(entered_species$attribution[sp]),b,sep="")
      entered_species$attribution[sp]<-as.character(b)
		  }
		latest_opinion[sp]<-paste(opinions$author1last[cites],opinions$pubyr[cites],sep=" ")
		}	else {
			if (entered_species$attribution[sp]=="")	{
				entered_species$attribution[sp]<-"Information Lost"
			}
			latest_opinion[sp]<-"Information Lost"	
		}
	}

# find the reference responsible for the most finds for entered & unentered species
entered_species<-cbind(entered_species,latest_opinion)
dominant_ref_entered_sp<-vector(length=entered)
for (sp in 1:entered)	{
	if (entered_species$associated_records[sp]==1)	{
		dominant_ref_entered_sp[sp]<-occurrences$reference_no[match(entered_species$taxon_name[sp], occurrences$taxon_name)]
		if (is.na(dominant_ref_entered_sp[sp]))
			dominant_ref_entered_sp[sp]<-occurrences$reference_no[match(entered_species$taxon_name[sp], occurrences$matched_name)]
	} else {
		temp_refs<-subset(occurrences$reference_no, occurrences$taxon_name==entered_species$taxon_name[sp])
		aaa<-table(temp_refs)
		if (dim(aaa)>1)	{
			bref<-temp_refs[1]
			mref<-aaa[[1]]
			m<-1
			for (j in 2:dim(aaa))	{
				if (aaa[[j]]>mref)	{
					bref<-temp_refs[m]
					mref<-aaa[[j]]
					}
				m<-m+aaa[[j]]
				}
			dominant_ref_entered_sp[sp]<-bref
			} else {
			dominant_ref_entered_sp[sp]<-temp_refs[1]
			}
		}	# end case where there are multiple occurrences
	}

dominant_ref_entered_sp<-clear_na_from_vector(dominant_ref_entered_sp, 0)
uniq_refs_ent<-unique(subset(dominant_ref_entered_sp, dominant_ref_entered_sp>1))
urefs<-length(uniq_refs_ent)
#refs_ent<- pbdb_references(id=uniq_refs_ent, vocab="pbdb", show="formatted",limit="all")
for (u in 1:urefs)  {
  httpR<-paste("http://paleobiodb.org/data1.1/refs/single.txt?id=",uniq_refs_ent[u],"&show=both",sep="")
  if (u==1) {
    refs_ent<-read.table(httpR, sep=',', header=T)
    } else {
    nr<-read.table(httpR, sep=',', header=T)
    refs_ent<-rbind(refs_ent,nr)
    }
  }
refs_ent$doi<-clear_na_from_vector(refs_ent$doi,"unentered")
refs_ent<-clear_na_from_matrix(refs_ent,"")
dominant_occurrence_ref_entered<-vector(length=entered)

for (sp in 1:entered)	{
	if (dominant_ref_entered_sp[sp]>0)  {
    dominant_occurrence_ref_entered[sp]<-as.character(refs_ent$formatted[match(dominant_ref_entered_sp[sp], uniq_refs_ent)])
	  } else {
    dominant_occurrence_ref_entered[sp]<-"Information Missing"  
	  }
  }
#dominant_ref_data_ent <- as.character(pbdb_references(id=dput(unique(dominant_ref_entered_sp)), show="formatted", limit="all"))
entered_species<-cbind(entered_species,dominant_occurrence_ref_entered)

# NOW, do this for the un-entered species!!!
unentered<-dim(unentered_species)[1]
#unentered_species<-cbind(unentered_species,latest_opinion)
dominant_ref_unentered_sp<-vector(length=unentered)
for (sp in 1:unentered)  {
  if (unentered_species$associated_records[sp]==1)  {
    dominant_ref_unentered_sp[sp]<-occurrences$reference_no[match(unentered_species$taxon_name[sp], occurrences$taxon_name)]
    if (is.na(dominant_ref_unentered_sp[sp]))
      dominant_ref_unentered_sp[sp]<-occurrences$reference_no[match(unentered_species$taxon_name[sp], occurrences$matched_name)]
    } else {
    temp_refs<-subset(occurrences$reference_no, occurrences$taxon_name==unentered_species$taxon_name[sp])
    aaa<-table(temp_refs)
    if (dim(aaa)>1)	{
      bref<-temp_refs[1]
      mref<-aaa[[1]]
      m<-1
      for (j in 2:dim(aaa))	{
        if (aaa[[j]]>mref)	{
          bref<-temp_refs[m]
          mref<-aaa[[j]]
          }
        m<-m+aaa[[j]]
        }
      dominant_ref_unentered_sp[sp]<-bref
      } else {
      dominant_ref_unentered_sp[sp]<-temp_refs[1]
      }
    }	# end case where there are multiple occurrences
  }
dominant_ref_unentered_sp<-clear_na_from_vector(dominant_ref_unentered_sp, 0)
uniq_refs_unent<-unique(subset(dominant_ref_unentered_sp, dominant_ref_unentered_sp>1))
urefs<-length(uniq_refs_unent)
#refs_unent<-pbdb_references(id=uniq_refs_unent, show="formatted",limit="all")
for (u in 1:urefs)  {
  httpR<-paste("http://paleobiodb.org/data1.1/refs/single.txt?id=",uniq_refs_ent[u],"&show=both",sep="")
  if (u==1) {
    refs_unent<-read.table(httpR, sep=',', header=T)
    } else {
    nr<-read.table(httpR, sep=',', header=T)
    refs_unent<-rbind(refs_unent,nr)
    }
  }
refs_unent$doi<-clear_na_from_vector(refs_unent$doi,"unentered")
refs_unent<-clear_na_from_matrix(refs_unent,"")

dominant_occurrence_ref_unentered<-vector(length=unentered)
for (sp in 1:unentered)  {
  if (dominant_ref_unentered_sp[sp]>0) {
    dominant_occurrence_ref_unentered[sp]<-as.character(refs_unent$formatted[match(dominant_ref_unentered_sp[sp], uniq_refs_unent)])
    } else {
    dominant_occurrence_ref_unentered[sp]<-"Information Lost"  
    }
  }
#dominant_ref_data_unent <- as.character(pbdb_references(id=dput(unique(dominant_ref_unentered_sp)), show="formatted", limit="all"))
unentered_species<-cbind(unentered_species,dominant_occurrence_ref_unentered)

# set up file name.  If just one time unit used, then use that name instead of oldest and youngest time unit names
if (oldest!=youngest)	{
	time_range<-paste(oldest,"-",youngest,sep="")
	} else {
		time_range<-oldest
	}
filename1<-paste("Entered",taxon,time_range,start_date,"to",end_date,sep="_")
filename2<-paste("Unentered",taxon,time_range,start_date,"to",end_date,sep="_")
filename1<-paste(filename1,file_format,sep="")
filename2<-paste(filename2,file_format,sep="")
# tab delimited or comma delimited
if (file_format!=".csv")	{
	write.table(entered_species,filename1,sep="\t",eol="\n",row.names=FALSE)
	write.table(unentered_species,filename2,sep="\t",eol="\n",row.names=FALSE)
	}
if (file_format==".csv")	{
	write.table(output1,filename1,sep=",",eol="\n",row.names=FALSE)
	write.table(output2,filename2,sep=",",eol="\n",row.names=FALSE)
	}
}

paleodb_find_gaps<-function(taxon,oldest,youngest,signif_gap,file_format)
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
	
	# get genus names
	httpG<-paste("http://paleobiodb.org/data1.1/taxa/list.txt?name=",taxon,"&rel=all_children&show=attr&rank=genus,subgenus&limit=9999",sep="")
	genera<-read.table(httpG, sep=',', header=T)
	httpF<-paste("http://paleobiodb.org/data1.1/occs/list.txt?base_name=",taxon,"&max_ma=",start_ma,"&min_ma=",end_ma,"&show=time,loc,crmod&limit=all",sep="")
	occurrences <-read.table(httpF, sep=',', header=T)
	
	# remove n. sp., aff., etc.
  occurrences$taxon_name<-paleodb_species_occurrence_name_cleaning(occurrences$taxon_name)
	occurrences$matched_name<-paleodb_species_occurrence_name_cleaning(occurrences$matched_name)
	occurrences$late_age<-clear_na_from_vector(occurrences$late_age,0)
	
	# sometimes taxon_name comes up NA; replace it with matched_name
  occurrences<-clear_matrix_na_with_another_cell_value(occurrences,match("taxon_name",colnames(occurrences)),match("matched_name",colnames(occurrences)))
	
	mid_age<-(occurrences$early_age+occurrences$late_age)/2
	occurrences<-cbind(occurrences,mid_age)
	finds<-length(occurrences$matched_name)
	found_genus<-vector(length=finds)
	# separate genus name
  for (i in 1:finds)  {
		found_genus[i]<-strsplit(as.character(occurrences$matched_name[i])," ")[[1]][1]
  	}
	occurrences<-cbind(occurrences,found_genus)
	gap_ma<-vector(length=finds)
	occurrences<-cbind(occurrences,gap_ma)
	collection_no_older<-occurrences$collection_no
	matched_name_older<-occurrences$matched_name
	early_age_older<-occurrences$early_age
	late_age_older<-occurrences$late_age
	early_interval_older<-occurrences$early_interval
	occurrences<-cbind(occurrences,collection_no_older,matched_name_older,early_age_older,late_age_older,early_interval_older)
	sortdata<-c("found_genus","mid_age")
	# flip-flop age so that we can get oldest members first
	occurrences$mid_age<--1*occurrences$mid_age
	occurrences<-occurrences[do.call("order",occurrences[sortdata]),]
	occurrences$mid_age<--1*occurrences$mid_age
	
	gaps<-1
	k<-0
	for (i in 1:(finds-1))  {
		j<-i+1
		k<-k+1
		if ((occurrences$taxon_rank[i]=="genus" || occurrences$taxon_rank[i]=="subgenus") || occurrences$taxon_rank[i]=="species")	{
			if (as.character(occurrences$found_genus[i])==as.character(occurrences$found_genus[j]) && (occurrences$late_age[i]-occurrences$early_age[j])>signif_gap)  {
				gap_length<-occurrences$late_age[j]-occurrences$early_age[i]
				if (gaps==1)  {
					occurrences$gap_ma[i]<-as.numeric(occurrences$late_age[j])-as.numeric(occurrences$early_age[i])
					occurrences$collection_no_older[i]<-occurrences$collection_no[j]
					# add "sp." to genus or subgenus only ids.
					if (strsplit(as.character(occurrences$matched_name_older[i])," ")[[1]][1]==1)	{
						occurrences$matched_name_older[i]<-paste(occurrences$matched_name[j],"sp.",sep=" ")
					} else if (strsplit(as.character(occurrences$matched_name_older[i])," ")[[1]][1]==2 && occurrences$matched_rank[j]=="subgenus")	{
						occurrences$matched_name_older[i]<-paste(occurrences$matched_name[j],"sp.",sep=" ")	
					}	else {
						occurrences$matched_name_older[i]<-occurrences$matched_name[j]
					}
					occurrences$early_age_older[i]<-occurrences$early_age[j]
					occurrences$late_age_older[i]<-occurrences$late_age[j]
					occurrences$early_interval_older[i]<-occurrences$early_interval[j]   	
					genus_gaps<-occurrences[i,]
					gaps<-gaps+1
				} else {
					occurrences$gap_ma[i]<-as.numeric(occurrences$late_age[j])-as.numeric(occurrences$early_age[i])
					occurrences$gap_ma[i]<-as.numeric(occurrences$late_age[j])-as.numeric(occurrences$early_age[i])
					occurrences$collection_no_older[i]<-occurrences$collection_no[j]
					occurrences$matched_name_older[i]<-occurrences$matched_name[j]
					occurrences$early_age_older[i]<-occurrences$early_age[j]
					occurrences$late_age_older[i]<-occurrences$late_age[j]
					occurrences$early_interval_older[i]<-occurrences$early_interval[j]   	
					genus_gaps<-rbind(genus_gaps,occurrences[i,])
					gaps<-gaps+1
				}
				last_gap<-j
				#		Debugging routine
				#    flurb<-paste(k,i,gaps,genus_gaps$gap_ma[gaps-1],genus_gaps$found_genus[gaps-1],occurrences$found_genus[i],sep=" ")
				#    print(flurb)
			}
		}  
	}
	output<-as.character(genus_gaps$found_genus)
	chch<-"genus"
	output<-cbind(output,as.character(genus_gaps$matched_name_older))
	chch<-cbind(chch,"taxon_name_younger")
	output<-cbind(output,genus_gaps$collection_no_older)
	chch<-cbind(chch,"collection_no_younger")
	output<-cbind(output,as.character(genus_gaps$early_interval_older))
	chch<-cbind(chch,"interval_younger")
	output<-cbind(output,as.character(genus_gaps$early_age_older))
	chch<-cbind(chch,"early_age_younger")
	output<-cbind(output,as.character(genus_gaps$late_age_older))
	chch<-cbind(chch,"late_age_younger")
	output<-cbind(output,as.character(genus_gaps$matched_name))
	chch<-cbind(chch,"taxon_name_older")
	output<-cbind(output,genus_gaps$collection_no)
	chch<-cbind(chch,"collection_no_older")
	output<-cbind(output,as.character(genus_gaps$early_interval))
	chch<-cbind(chch,"interval_older")
	output<-cbind(output,as.character(genus_gaps$early_age))
	chch<-cbind(chch,"early_age_older")
	output<-cbind(output,as.character(genus_gaps$late_age))
	colnames(output)<-cbind(chch,"late_age_older")
	#output<-cbind(output,genus_gaps$collection_no_older,as.character(genus_gaps$early_interval_older),as.character(genus_gaps$early_age_older),as.character(genus_gaps$late_age_older),as.character(genus_gaps$matched_name),genus_gaps$collection_no,as.character(genus_gaps$early_interval),as.character(genus_gaps$early_age),as.character(genus_gaps$late_age))
	
	#colnames(output)<-cbind("genus","taxon_name_older","collection_no_older","interval_older","early_age_older","late_age_older","taxon_name_younger","collection_no_younger","interval_younger","early_age_younger","late_age_younger")
	
	if (oldest!=youngest)	{
		time_range<-paste(oldest,"-",youngest,sep="")
	} else {
		time_range<-oldest	
	}
	filename<-paste("Suspicious_Gaps",taxon,time_range,sep="_")
	filename<-paste(filename,file_format,sep="")
	if (file_format!=".csv")	
		write.table(output,filename,sep="\t",eol="\n",row.names=FALSE)
	if (file_format==".csv")	
		write.table(output1,filename1,sep=",",eol="\n",row.names=FALSE)
}

# Collect information about rock units to make stratigraphy consistent
paleodb_vett_formations<-function(taxon,oldest,youngest,today,file_format)
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
	
	request<-paste("base_name=",as.character(taxon[1]),sep="")
	gr<-length(taxon)
	if (gr>1)	{
		for (t in 2:gr)	{
			request<-paste(request,"&base_name=",sep="")
			request<-paste(request,as.character(taxon[t]),sep="")
			}
		}
	
#	http<-paste("http://paleobiodb.org/data1.1/colls/list.txt?base_name=",taxon,"&max_ma=",start_ma,"&min_ma=",end_ma,"&show=time,stratext,crmod&limit=all",sep="")
	http<-paste("http://paleobiodb.org/data1.1/colls/list.txt?",request,"&max_ma=",start_ma,"&min_ma=",end_ma,"&show=time,stratext,crmod&limit=all",sep="")
	collections <-read.table(http, sep=',', header=T)
	collections<-clear_na_from_matrix(collections, "")
	lists<-dim(collections)[1]
	mid_age<-(collections$early_age+collections$late_age)/2
	collections<-cbind(collections,mid_age)
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
	
	rock_unit<-paste(collections$formation,collections$member,sep=" - ")
	rock_unit <- gsub(", NA", "",rock_unit)
	rock_unit <- gsub(" \\?" ,"",rock_unit)
	collections<-cbind(collections,rock_unit)
	
	sortdata<-c("rock_unit","mid_age","early_interval","zone")
	collections<-collections[do.call("order",collections[sortdata]),]
	
	reduced_collections<-cbind(as.character(collections$formation),as.character(collections$member),as.character(collections$early_interval),as.character(collections$late_interval),as.character(collections$zone),as.numeric(collections$early_age),as.numeric(collections$late_age),as.numeric(collections$mid_age))
	colnames(reduced_collections)<-cbind("formation","member","early_interval","late_interval","zone","early_age","late_age","mid_age")
	b<-1
	while (collections$formation[b]=="" && collections$member[b]=="")
		b<-b+1
	
	orphans<-collections[1:b,]
	
	reduced_collections[,7]<-clear_na_from_vector(reduced_collections[,7],0)
	reduced_collections[,8]<-clear_na_from_vector(reduced_collections[,8],0)
	formation_info<-reduced_collections[b,]
	k<-1
	for (c in (b+1):lists)	{
		d<-c-1
		if (sum(as.numeric(reduced_collections[c,]!=reduced_collections[d,]))>0 && reduced_collections[c,1]!=" - ")	{
			formation_info<-rbind(formation_info,reduced_collections[c,])
			k<-k+1		
			}	# end case of rock units with different information
		}	# end search for unique lists
	
	if (oldest!=youngest)	{
		time_range<-paste(oldest,"-",youngest,sep="")
		} else {
		time_range<-oldest
		}
	
	filename1<-paste("Rock_Unit_Information",taxon,time_range,sep="_")
	filename1<-paste(filename1,file_format,sep="")
	if (file_format!=".csv")	
		write.table(formation_info,filename1,sep="\t",eol="\n",row.names=FALSE)
	if (file_format==".csv")	
		write.table(formation_info,filename1,sep=",",eol="\n",row.names=FALSE)
	
	filename2<-paste("Formationless_Collections",taxon,time_range,sep="_")
	filename2<-paste(filename2,file_format,sep="")
	if (file_format!=".csv")	
		write.table(orphans,filename2,sep="\t",eol="\n",row.names=FALSE)
	if (file_format==".csv")	
		write.table(orphans,filename2,sep=",",eol="\n",row.names=FALSE)
}

#taxon_list_name<-"Hone_Species.txt"
#taxon_set<-"Hone"
#file_format<-".xls"
get_paleodb_occurrence_for_taxon_list<-function(taxon_set, taxon_list_name, file_format)
{ 
taxon_list<-read.table(file=taxon_list_name, sep="\t", header=FALSE, stringsAsFactors=TRUE)
species<-dim(taxon_list)[1]
samp<-0
taxon_records<-vector(length=species)
for (sp in 3:species) {
	spc<-as.character(taxon_list[sp,1])
	spc<-gsub(" ","%20",spc)
  http<-paste("http://paleobiodb.org/data1.1/occs/list.txt?base_name=",spc,"&show=time,loc&limit=all",sep="")
  taxon_finds<-read.table(http, sep=',', header=T)
	taxon_records[sp]<-dim(taxon_finds)[1]
  if (taxon_records[sp]>0)  {    
    if (samp==0)  {
      occurrences<-read.table(http, sep=',', header=T)
      } else {
      occurrences<-rbind(occurrences,taxon_finds)  
      }
    samp<-samp+taxon_records[sp]
		}
  }

localities<-unique(occurrences$collection_no)
localities<-sort(localities,FALSE)
locales<-length(localities)
for (l in 1:locales)  {
  httpC<-paste("http://www.paleobiodb.org/data1.1/colls/list.txt?id=",localities[l],"&show=loc,time,strat,stratext",sep="")
  if (l==1) {
    collection_info<-read.table(httpC, sep=',', header=T)
    } else {
      nc<-read.table(httpC, sep=',', header=T)
      collection_info<-rbind(collection_info,nc)
    }
  }

recordless_taxa<-subset(taxon_list, taxon_records==0)

filename1<-paste("Recordless",taxon_set,sep="_")
filename1<-paste(filename1,file_format,sep="")

filename2<-paste(taxon_set,"Localities",sep="_")
filename2<-paste(filename2,file_format,sep="")

filename3<-paste(taxon_set,"Records",sep="_")
filename3<-paste(filename3,file_format,sep="")

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


find_species_name_within_higher_taxon<-function(clade,species,sub)
{
httpM<-paste("http://www.paleobiodb.org/data1.1/taxa/list.txt?name=",clade,"&rel=all_children&show=attr,app",sep="")
members<-read.table(httpM, sep=',', header=T)
member_species<-subset(members, members$rank=="species")
if (sub==TRUE)
	member_species<-rbind(member_species,subset(members, members$rank=="subspecies"))

sp_ct<-dim(member_species)[1]

matches<-0
if (sp_ct>0){
  for (s in 1:sp_ct)	{
  	words<-length(strsplit(as.character(member_species$taxon_name[s])," ")[[1]])
  	species_name<-strsplit(as.character(member_species$taxon_name[s])," ")[[1]][words]
  	if (species_name==species)
  		matches<-matches+1
  	}

  if (matches>0)	{	
  	poss_species<-vector(length=matches)
  	m<-0
  	for (s in 1:sp_ct)	{
  		words<-length(strsplit(as.character(member_species$taxon_name[s])," ")[[1]])
  		species_name<-strsplit(as.character(member_species$taxon_name[s])," ")[[1]][words]
  		if (species_name==species)	{
  			m<-m+1
  			poss_species[m]<-as.character(member_species$taxon_name[s])
  			}
  		}
  	} else {
  	  poss_species<-vector(length=1)
  	  poss_species[1]<-paste("no matches")
  	  
  	}
  } else {
    poss_species<-vector(length=1)
    poss_species[1]<-paste(clade," not in database")
  }

return(poss_species)
}


find_taxa_assigned_to_invalid_higher_taxon<-function(taxon,file_format)
{
httpO<-paste("http://paleobiodb.org/data1.2/taxa/opinions.txt?base_name=",taxon,"&order=hierarchy&limit=all",sep="")
opinions<-read.table(httpO, sep=',', header=T)
#opinions[,13]<-clear_na_from_vector(opinions[,13],"")
#opinions[,13]<-as.numeric(opinions[,13])
taxa<-dim(opinions)[1]
statuses<-vector(length=max(opinions$orig_no))
parent_status<-vector(length=taxa)
parent<-vector(length=max(opinions$orig_no))
parent_no<-vector(length=max(opinions$orig_no))
authority<-vector(length=max(opinions$orig_no))
probs<-0
for (t in 1:taxa) {
  tn<-as.numeric(opinions$child_spelling_no[t])
  statuses[tn]<-as.character(opinions$status[t])
  tn<-as.numeric(opinions$orig_no[t])
  statuses[tn]<-as.character(opinions$status[t])
  parent[tn]<-as.character(opinions$parent_name[t])
  parent_no[tn]<-as.numeric(opinions$parent_no[t])
  authority[tn]<-paste(as.character(opinions$author[t]),as.character(opinions$pubyr[t]),sep=" ")
  if (t>1 && opinions$status[t]=="belongs to") {
    ht<-as.numeric(opinions$parent_no[t])
    if (statuses[ht]=="belongs to") {
      parent_status[t]<-"No worries"
      } else if (statuses[ht]=="replaced by")  {
      pt<-as.numeric(opinions$parent_no[t])
      parent_status[t]<-paste(as.character(opinions$parent_name[t]),"has been replaced by",parent[pt],"by",authority[ht],sep=" ")
      probs<-probs+1
      } else if (statuses[ht]=="subjective synonym of")  {
      pt<-as.numeric(opinions$parent_no[t])
      parent_status[t]<-paste(as.character(opinions$parent_name[t]),"was subjectively synonymized with",parent[pt],"by",authority[ht],sep=" ")
      probs<-probs+1
      } else if (statuses[ht]=="objective synonym of")  {
      pt<-as.numeric(opinions$parent_no[t])
      parent_status[t]<-paste(as.character(opinions$parent_name[t]),"was objectively synonymized with",parent[pt],"by",authority[ht],sep=" ")
      probs<-probs+1
      } else if (statuses[ht]=="nomen dubium" || (statuses[ht]=="nomen nudum" || statuses[ht]=="nomen oblitum")) {
        parent_status[t]<-paste(as.character(opinions$parent_name[t]),"is now a",statuses[ht],"according to",authority[ht],sep=" ")
        probs<-probs+1
      } else if (statuses[ht]=="invalid subgroup of") {
      parent_status[t]<-paste(as.character(opinions$parent_name[t]),"is now an invalid subgroup of",parent[pt],"according to",authority[ht],sep=" ")
      probs<-probs+1
      } else {
      parent_status[t]<-paste(as.character(opinions$parent_name[t])," is ",statuses[ht],sep="")
      probs<-probs+1
      }
    }
  }
output<-matrix(0,probs,3)
prc<-0
for (t in 2:taxa) {
  if (opinions$status[t]=="belongs to" && parent_status[t]!="No worries"){
    prc<-prc+1
    output[prc,1]<-as.numeric(opinions$orig_no[t])
    output[prc,2]<-as.character(opinions$taxon_name[t])
    output[prc,3]<-as.character(parent_status[t])
#    print(c(as.numeric(opinions$orig_no[t]),as.character(opinions$taxon_name[t]),as.character(parent_status[t])))
    }
  }
colnames(output)<-c("orig_no","taxon_name","issue")
problem_children<-paste(taxon,"Problem_Children",sep="_")
problem_children<-paste(problem_children,file_format,sep="")

if (file_format!=".csv")  {
  write.table(output,problem_children,sep="\t",eol="\n",row.names=FALSE, col.names=TRUE)
  }
if (file_format==".csv")	{
  write.table(output,problem_children,sep=",",eol="\n",row.names=FALSE,col.names=TRUE)
  }

}


find_taxa_assigned_to_invalid_higher_taxa<-function(taxon,file_format=".xls")
{

httpO<-paste("http://paleobiodb.org/data1.2/taxa/opinions.txt?base_name=",taxon,"&order=hierarchy&limit=all",sep="")
opinions<-read.csv(httpO)
x<-opinions

val <- as.numeric(sapply(opinions$parent_no, function(x) which(x==opinions$orig_no)))

name_valid <- data.frame(id=opinions$orig_no, taxon_name=opinions$taxon_name, issue=paste(opinions$taxon_name[val],opinions$status[val],opinions$parent_name[val],"according to",opinions$author[val],opinions$pubyr[val]))

output <- name_valid[which(opinions$status[val] != "belongs to"),]

opinion_status <- opinions$status[match(output$taxon_name, opinions$taxon_name)]

output <- output[which(opinion_status=="belongs to"),]
  
problem_children<-paste(taxon,"Problem_Children",sep="_")
problem_children<-paste(problem_children,file_format,sep="")
  
if (file_format!=".csv")  {
	write.table(output,problem_children,sep="\t",eol="\n",row.names=FALSE, col.names=TRUE)
  }
if (file_format==".csv")        {
	write.csv(output,problem_children, row.names=F)
  }
  
}