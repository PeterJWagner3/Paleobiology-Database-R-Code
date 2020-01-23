# brute force kluges to make R what I effin' well want it to do..... 

matrix_to_list_by_columns <- function(data_mat)	{
output <- tapply(data_mat,rep(1:ncol(data_mat),each=nrow(data_mat)),function(i)i);
if (length(colnames(data_mat))>0)	
	names(output) <- colnames(data_mat);
return(output);
}

matrix_to_list_by_rows <- function(data_mat)	{
output <- tapply(data_mat,rep(1:nrow(data_mat),each=ncol(data_mat)),function(i)i);
if (length(colnames(data_mat))>0)	
	names(output) <- rownames(data_mat);
return(output);
}

scourgify_web_text <- function(web_text)	{
web_text <- gsub("–","-",web_text);
web_text <- gsub("“","\"",web_text);
web_text <- gsub("”","\"",web_text);
web_text <- gsub("‘","\'",web_text);
web_text <- gsub("’","\'",web_text);
web_text <- gsub("‚Äì","–",web_text);
web_text <- gsub("‚Äî","-",web_text);
web_text <- gsub("‚Äú","“",web_text);
web_text <- gsub("‚Äù","”",web_text);
web_text <- gsub("‚Äò","‘",web_text);
web_text <- gsub("‚Äô","’",web_text);
web_text <- gsub("é","e",web_text);
web_text <- gsub("√©","e",web_text);
#web_text <- gsub("ö","&ouml;",web_text);
#web_text <- gsub("ü","&uuml;",web_text);
web_text <- gsub("ü","u",web_text);
web_text <- gsub("√º","u",web_text);
web_text <- gsub("ý","y",web_text);
web_text <- gsub("√Ω","y",web_text);
return(web_text);
}

frak_it <- function(x)	{
return(as.numeric(as.character(x)))
}

frak_it_character <- function(x)	{
if (is.factor(x)) {
	return(as.character(x));
	}	else {
	return(x);
	}
}

transmogrify_to_title_case <- function(name) {
name_part <- strsplit(name, " ")[[1]]
return (paste(toupper(substring(name_part, 1,1)), substring(name_part, 2),sep="", collapse=" "));
}

transmogrify_diacritics <- function(funky_text)	{
j <- strsplit(as.character(funky_text),split="",fixed=TRUE)[[1]];
eek <- c("à","á","â","ã","ā","ă","ȧ","ä","ả","å","ǎ","ȁ","ȃ","ą","ạ","ḁ","ẚ","ầ","ấ","ẫ","ẩ","ằ","ắ","ẵ","ẳ","ǡ","ǟ","ǻ","ậ","ặ");
eekeek <- (1:length(j))[j %in% eek];
j[eekeek] <- "a";
eekeek <- (1:length(j))[j %in% toupper(eek)];
j[eekeek] <- "A";
eek <- c("æ","ǽ","ǣ");
eekeek <- (1:length(j))[j %in% eek];
j[eekeek] <- "ae";
eekeek <- (1:length(j))[j %in% toupper(eek)];
j[eekeek] <- "Ae";
eek <- c("ć","ĉ","ċ","č","ƈ","ç","ḉ","ȼ");
eekeek <- (1:length(j))[j %in% eek];
j[eekeek] <- "c";
eekeek <- (1:length(j))[j %in% toupper(eek)];
j[eekeek] <- "C";
eek <- c("ḋ","ɗ","ḍ","ḏ","ḑ","ḓ","ď");
eekeek <- (1:length(j))[j %in% eek];
j[eekeek] <- "d";
eekeek <- (1:length(j))[j %in% toupper(eek)];
j[eekeek] <- "D";
eek <- c("è","é","ê","ẽ","ē","ĕ","ė","ë","ẻ","ě","ȅ","ȇ","ẹ","ȩ","ę","ḙ","ḛ","ề","ế","ễ","ể","ḕ","ḗ","ệ","ḝ");
eekeek <- (1:length(j))[j %in% eek];
j[eekeek] <- "e";
eekeek <- (1:length(j))[j %in% toupper(eek)];
j[eekeek] <- "E";
eek <- c("ǵ","ĝ","ḡ","ğ","ġ","ǧ","ɠ","ģ","ǥ");
eekeek <- (1:length(j))[j %in% eek];
j[eekeek] <- "g";
eekeek <- (1:length(j))[j %in% toupper(eek)];
j[eekeek] <- "G";
eek <- c("ĥ","ḣ","ȟ","ḥ","ḩ","ḫ","ẖ","ħ","ⱨ");
eekeek <- (1:length(j))[j %in% eek];
j[eekeek] <- "h";
eekeek <- (1:length(j))[j %in% toupper(eek)];
j[eekeek] <- "H";
eek <- c("ì","í","î","ĩ","ī","ĭ","ı","ï","ỉ","ǐ","ị","į","ȉ","ȋ","ḭ","ɨ","ḯ");
eekeek <- (1:length(j))[j %in% eek];
j[eekeek] <- "i";
eekeek <- (1:length(j))[j %in% toupper(eek)];
j[eekeek] <- "I";
eek <- c("ĵ","ǰ","ȷ","ɉ");
eekeek <- (1:length(j))[j %in% eek];
j[eekeek] <- "j";
eekeek <- (1:length(j))[j %in% toupper(eek)];
j[eekeek] <- "J";
eek <- c("ḱ","ǩ","ḵ","ƙ","ḳ","ĸ","ⱪ");
eekeek <- (1:length(j))[j %in% eek];
j[eekeek] <- "k";
eekeek <- (1:length(j))[j %in% toupper(eek)];
j[eekeek] <- "K";
eek <- c("ĺ","ľ","ŀ","ł","ƚ");
eekeek <- (1:length(j))[j %in% eek];
j[eekeek] <- "l";
eekeek <- (1:length(j))[j %in% toupper(eek)];
j[eekeek] <- "L";
eek <- c("ḿ","ṁ","ṃ","ɱ");
eekeek <- (1:length(j))[j %in% eek];
j[eekeek] <- "m";
eekeek <- (1:length(j))[j %in% toupper(eek)];
j[eekeek] <- "M";
eek <- c("ǹ","ń","ñ","ṅ","ň");
eekeek <- (1:length(j))[j %in% eek];
j[eekeek] <- "n";
eekeek <- (1:length(j))[j %in% toupper(eek)];
j[eekeek] <- "N";
eek <- c("ò","ó","ô","õ","ō","ŏ","ȯ","ö","ỏ","ő","ǒ","ȍ","ȏ","ơ","ǫ","ọ","ɵ","ø");
eekeek <- (1:length(j))[j %in% eek];
j[eekeek] <- "o";
eekeek <- (1:length(j))[j %in% toupper(eek)];
j[eekeek] <- "O";
eek <- c("ṕ","ṗ");
eekeek <- (1:length(j))[j %in% eek];
j[eekeek] <- "p";
eekeek <- (1:length(j))[j %in% toupper(eek)];
j[eekeek] <- "P";
eek <- c("ŕ","ṙ","ř","ȑ","ȓ","ṛ","ŗ","ṟ","ṝ","ɍ");
eekeek <- (1:length(j))[j %in% eek];
j[eekeek] <- "r";
eekeek <- (1:length(j))[j %in% toupper(eek)];
j[eekeek] <- "R";
eek <- c("ś","ŝ","ṡ","š","ṣ","ș","ş","ȿ","ṥ","ṧ");
eekeek <- (1:length(j))[j %in% eek];
j[eekeek] <- "s";
eekeek <- (1:length(j))[j %in% toupper(eek)];
j[eekeek] <- "S";
eek <- c("ß");
eekeek <- (1:length(j))[j %in% eek];
j[eekeek] <- "ss";
eek <- c("ṫ","ẗ","ť","ƫ","ṭ","ț","ţ","ṱ","ṯ","ŧ");
eekeek <- (1:length(j))[j %in% eek];
j[eekeek] <- "t";
eekeek <- (1:length(j))[j %in% toupper(eek)];
j[eekeek] <- "T";
eek <- c("ù","ú","û","ũ","ū","ŭ","ü","ủ","ű","ǔ","ȕ","ȗ","ụ","ṳ","ų","ṷ","ṵ");
eekeek <- (1:length(j))[j %in% eek];
j[eekeek] <- "u";
eekeek <- (1:length(j))[j %in% toupper(eek)];
j[eekeek] <- "U";
eek <- c("ṽ","ṿ","ⱱ");
eekeek <- (1:length(j))[j %in% eek];
j[eekeek] <- "v";
eekeek <- (1:length(j))[j %in% toupper(eek)];
j[eekeek] <- "V";
eek <- c("ẁ","ẃ","ŵ","ẇ","ẅ","ẘ","ẉ");
eekeek <- (1:length(j))[j %in% eek];
j[eekeek] <- "w";
eekeek <- (1:length(j))[j %in% toupper(eek)];
j[eekeek] <- "W";
eek <- c("ẋ","ẍ");
eekeek <- (1:length(j))[j %in% eek];
j[eekeek] <- "x";
eekeek <- (1:length(j))[j %in% toupper(eek)];
j[eekeek] <- "X";
eek <- c("ỳ","ý","ŷ","ȳ","ẏ","ÿ","ỷ","ẙ","ɏ");
eekeek <- (1:length(j))[j %in% eek];
j[eekeek] <- "y";
eekeek <- (1:length(j))[j %in% toupper(eek)];
j[eekeek] <- "Y";
eek <- c("ź","ẑ","ż","ž","ȥ","ẓ","ẕ","ƶ");
eekeek <- (1:length(j))[j %in% eek];
j[eekeek] <- "z";
eekeek <- (1:length(j))[j %in% toupper(eek)];
j[eekeek] <- "Z";
j <- j[j!="̈"];
j <- j[j!="̀"];
j <- j[j!="̌"];
return(paste(j,collapse=""));
}

evanesco_na_from_matrix <- function(data, replacement="")  {
for (i in 1:ncol(data))	{
	if(sum(is.na(data[,i]))>0)	{
		duds <- (1:nrow(data))[is.na(data[,i])]
		data[duds,i] <- replacement
		}
	}
return(data)
}

evanesco_na_from_vector <- function(data, replacement="")	{
if(sum(is.na(data))>0)	{
	duds <- (1:length(data))[is.na(data)]
	data[duds] <- replacement
	}
return(data)
}

match_value_to_matrix_cell <- function(go_seek,hide)	{
row_matches <- c()
for (cc in 1:ncol(hide))	{
	rc <- match(go_seek,hide[,cc])
	if (!is.na(rc))
		row_matches <- c(row_matches,rc)
	}
if (is.null(row_matches))
	row_matches <- "";	
row_matches <- unique(row_matches);
return(row_matches);
}

insert_cell_into_vector_x <- function(x,new_value,cell_no)	{
# routine to put take vector x and add a new value at some point (cell_no)
# x = c(0,1,2,3,4,5), new_value=10, cell_no=3 gives x=c(0,1,10,2,3,4,5)
if (cell_no==1)	{
	return(c(new_value,x));
	} else if (cell_no==(length(x)+1))	{
	return(c(x,new_value));
	} else	{
	return(c(x[1:(cell_no-1)],new_value,x[cell_no:length(x)]));
	}
}

insert_row_into_matrix_x <- function(x,new_row,row_no)	{
# routine to put take matrix x and add a new row at that bumps other rows down one
# x = c(0,1,2,3,4,5), new_value=10, cell_no=3 gives x=c(0,1,10,2,3,4,5)
if (row_no==1)	{
	xx <- rbind(new_row,x);
	colnames(xx) <- colnames(x);
	return(xx);
	} else if (row_no==(nrow(x)+1))	{
	return(rbind(x,new_row));
	} else	{
	return(rbind(x[1:(row_no-1),],new_row,x[(row_no:nrow(x)),]));
#	return(c(x[1:(cell_no-1)],new_value,x[cell_no:length(x)]));
	}
}

#find_me <- paleodb_clean_member_basic_no_rock[need_info_member[3]]
matrix_row_match <- function(find_me,examined_matrix)	{
xxx <- which(examined_matrix==find_me,arr.ind=TRUE);
#examined_matrix[unique(xxx[,1]),]
if (nrow(xxx)==0)	{
	return(-1);
	} else if (length(unique(xxx[,1]))==1)	{
	return(xxx[1,1])
	} else {
	vvv <- unique(xxx[,1]);
	mx <- 0;
	mj <- 0;
	for (j in 1:length(vvv))	{
		if (sum(vvv[j]==xxx[,1])>mx)	{
			mj <- j;
			mx <- sum(vvv[j]==xxx[,1])
			} else if (sum(vvv[j]==xxx[,1])==mx)	{
			mj <- c(mj,j)
			}
		}
	if (length(mj)==1)	{
		return(vvv[mj]);
		} else {
		return(-2);
		}
	}
}

#character_array <- strsplit(rock_unit_name,split = "");
simplify2vector <- function(character_array)	{
return(simplify2array(character_array)[,1])
}

specialis_revelio_file_type <- function(filename)	{
file_info <- stringr::str_split(filename,pattern="")[[1]];
f_i <- length(file_info);
f_i_s <- 1+(1:f_i)[file_info %in% "."];
file_type <- c();
for (i in f_i_s:f_i)	file_type <- paste(file_type,file_info[i],sep="");
return(file_type);
}

# transform dataframe into a list based on some criterion.
#	list_criterion gives the column name with the variables that will be used to divide the data.frame into a list
listifor <- function(dataframe_to_transform,list_criterion)	{
criterion_col <- match(list_criterion,colnames(dataframe_to_transform));
list_categories <- sort(unique(dataframe_to_transform[,criterion_col]));
list_categories <- list_categories[list_categories!=""];

listed_dataframe <- list();
for (lc in 1:length(list_categories))	{
	this_subset <- subset(dataframe_to_transform,dataframe_to_transform[,criterion_col]==list_categories[lc]);
	listed_dataframe <- rlist::list.append(listed_dataframe,this_subset);
	}
names(listed_dataframe) <- list_categories;
return(listed_dataframe);
}

choose_with_logs <- function(N,n)	{
if (n==0 || n==N)	{
	return(1);
	} else	{
	return(exp(sum(log(1:N))-(sum(log(1:n))+sum(log(1:(N-n))))));
	}
}
