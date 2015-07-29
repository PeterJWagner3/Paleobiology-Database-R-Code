clear_na_from_matrix<-function(data, replacement)  {
  size<-dim(data)
  for (i in 1:size[1])	{
    for (j in 1:size[2]) if (is.na(data[i,j]))	data[i,j]<-replacement
  }
  return(data)
}

clear_matrix_na_with_another_cell_value<-function(data,j, k)	{
  size<-dim(data)
  for (i in 1:size[1])	{
    if (is.na(data[i,j]))	data[i,j]<-data[i,k]
  }
  return(data)
}

clear_na_from_vector<-function(data, replacement)	{
size<-length(data)
for (i in 1:size[1])	if (is.na(data[i]))	data[i]<-replacement
return(data)
}

paleodb_species_occurrence_name_cleaning<-function(taxon_name)
{
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
  taxon_name <- gsub("  " ," ",taxon_name)
  return(taxon_name)
}


clean_rock_unit_names<-function(named_rock_unit,delete_rock_type)
{
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
  #	named_rock_unit <- gsub("’","\\'",,named_rock_unit)
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
  #	named_rock_unit <- gsub("â€™","\\'",,named_rock_unit)
  
  if (delete_rock_type==TRUE || delete_rock_type==T)	{
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
}