default_minor_time_axis_size <- 4.285714285;

make_font_Arial <- function()	{
arial <- quartzFonts(arial = c("Arial","Arial-Bold","Arial-Italic","Avenir-BoldItalic"))
}

makeTransparent <- function(someColor, alpha=100)	{
newColor <- col2rgb(someColor)
apply(newColor, 2, function(curcoldata)X = {
	rgb(red=curcoldata[1], green=curcoldata[2],blue=curcoldata[3],alpha=alpha, maxColorValue=255)}
	)
}

shape_size_rosetta_stone <- function()	{
return(rosetta_stone_shape_size <- data.frame(size = c(1.5,1.0,2/3),circles=c(8.11,5.4,3.6),squares=c(7.18,4.79,3.19),triangles=c(10.905,7.27,4.846667),diamonds=c(10.15,6.765,4.51)));
}

accio_cex_size_for_shape <- function(cex_size,shape_pch)	{
rosetta_stone <- shape_size_rosetta_stone();
circle_pchs <- c(1,8,10,16,19,20,21);
square_pchs <- c(0,3,4,7,12,13,14,15,22);
triangle_pchs <- c(2,6,11,17,24,25);
diamond_pchs <- c(5,9,18,23);
area_to_match <- (rosetta_stone$squares[2]*cex_size)^2
if (!is.na(match(shape_pch,circle_pchs)))	{
	rosetta_circle <- match("circles",colnames(rosetta_stone));
	new_cex <- 2*sqrt(area_to_match/pi)/rosetta_stone[2,rosetta_circle];
	} else if (!is.na(match(shape_pch,triangle_pchs)))	{
	rosetta_triangle <- match("triangles",colnames(rosetta_stone));
#	h <- sqrt(((rosetta_stone[2,rosetta_triangle]/2)^2)+(rosetta_stone[2,rosetta_triangle]^2))
	b <- rosetta_stone[2,rosetta_triangle];
	tri_area <- sqrt(3)*(b^2)/4
	new_cex <- sqrt(area_to_match/tri_area)
#	tri_area2 <- sqrt(3)*((new_cex*b)^2)/4
#	sqrt(3)*(((sqrt(area_to_match/tri_area)*rosetta_stone[2,rosetta_triangle]))^2)/4
	} else if (!is.na(match(shape_pch,diamond_pchs)))	{
	rosetta_diamond <- match("diamonds",colnames(rosetta_stone));
	diagon <- rosetta_stone[2,rosetta_diamond];
	dia_area <- (diagon^2)/2
	new_cex <- sqrt(area_to_match/dia_area)
#	side <- sqrt((rosetta_stone[2,rosetta_diamond]^2)/2)
#	dia_area <- (side^2)
	} else	{
	new_cex <- cex_size;
	}
return(new_cex);
#rossetta_stone_shape_size[1,]/rossetta_stone_shape_size[2,]}
}

draw_time_scale <- function(time_scale,onset,end,strat_colors="",strat_names="",line_color="black",ord,height,strat_label_size,axis=1)	{
if (onset[1]<0 && time_scale[1]>0)
	time_scale <- -1*time_scale;
first_period <- sum(abs(time_scale)>=abs(onset));
last_period <- sum(abs(time_scale)>=abs(end));
if (length(strat_colors)==0)	strat_colors <- rep("white",last_period-first_period);

for (b in first_period:last_period)
	rect(max(onset,time_scale[b]),ord,min(end,time_scale[b+1]),height,col=as.character(strat_colors[b]),border=as.character(strat_colors[b]))

for (s in first_period:last_period)  segments(time_scale[s],ord,time_scale[s],height,lwd=0.5,col=line_color)
# put in decisive final line
if (!is.na(match(end,time_scale)))	segments(end,ord,end,height,col=line_color);

segments(onset,ord,end,ord,col=line_color,lwd=2)
segments(onset,height,end,col=line_color,height,lwd=2)

if (length(strat_names)>0)	{
	dummy <- mid <- c();
	for (i in first_period:last_period)	{
		mid <- c(mid,(max(onset,time_scale[i])+min(end,time_scale[(i+1)]))/2)
		dummy <- c(dummy,(height+ord)/2);
		}
	#add something
	text(mid[first_period:last_period],dummy[first_period:last_period],label=strat_names[first_period:last_period],cex=strat_label_size);
	}
}

# use this for simulated expectations (Supplementary figures)
Phanerozoic_Timescale_Plot <- function(onset=-541, end=0, time_scale, mxy, mny, ordinate, xsize=6, ysize=4.285714285, hues=FALSE, alt_back=FALSE, main_time_tick=100)	{
if (time_scale[1]>0)	time_scale <-  -1*time_scale;
if (onset>0)			onset <- -1*onset;
if (end>0)				end <- -1*end;
lst <- length(time_scale)
par(pin=c(xsize,ysize));
ages <- vector(length=1+abs(ceiling(onset/100))-abs(ceiling(end/100)))
st <- abs(ceiling(onset/100))+1
for (i in 1:length(ages))	{
	ages[i] <- 100*(st-i)
	}
plot(NA,type='n',axes=FALSE,main="",xlab="",ylab=ordinate,xlim=c(onset,end),ylim=c(mny,mxy));
#axis(1,at=seq(onset,end,by=1),tcl=-0.0,labels=FALSE,lwd=4/3)
axis(1,at=seq(100*ceiling(onset/100),end,by=100),tcl=-0.3,labels=ages,lwd=0,lwd.ticks=4/3)
axis(1,at=seq(100*ceiling(onset/100),end,by=50)[!seq(100*ceiling(onset/100),end,by=50)%in%seq(100*ceiling(onset/100),end,by=100)],tcl=-0.2,labels=FALSE,lwd=0,lwd.ticks=4/3)
axis(1,at=seq(10*round(onset/10),end,by=10)[!seq(10*round(onset/10),end,by=10)%in%seq(100*ceiling(onset/100),end,by=50)],tcl=-0.1,labels=FALSE,lwd=0,lwd.ticks=4/3)

exy <- (mxy-mny)/25
if (alt_back!=FALSE)	{
	for (b in 1:(lst-1))	{
		if (b%%2==0)	{
			rect(time_scale[b],(mny+exy),time_scale[b+1],mxy,col=alt_back,border=alt_back)
			}
		}
	}

#segments(onset,mny+exy,onset,mny-exy)
if (hues)	{
	per_colors <- Phanerozoic_Period_Colors();
	for (b in 1:(lst-1))	{
		rect(time_scale[b],(mny+exy),time_scale[b+1],(mny-exy),col=per_colors[b],border=per_colors[b])
		}
	}
for (s in 1:lst)  segments(max(onset,time_scale[s]),mny+exy,max(onset,time_scale[s]),mny-exy)
segments(onset,mny+exy,end,mny+exy,lwd=2)
segments(onset,mny-exy,end,mny-exy,lwd=2)

mid <- vector(length=(lst-1))
for (i in 1:(lst-1))    mid[i] <- (max(onset,time_scale[i])+time_scale[(i+1)])/2
text(mid,y=c(mny,mny,mny,mny,mny,mny,mny,mny,mny,mny,mny),label=c("Cm","O","S","D","C","P","Tr","J","K","Pg","Ng"),cex=0.9)
}

Phanerozoic_Timescale_Plot_Partial <- function(onset=-541, end=0, time_scale, mxy, mny, ordinate, stage_names, xsize=6, ysize=4.285714285, hues=FALSE, alt_back=FALSE, main_time_tick=100)	{
#quartzFonts(arial = c("Arial", "Arial Black", "Arial Oblique","Arial Black Oblique"))
if (time_scale[1]>0)	time_scale <-  -1*time_scale;
if (onset>0)			onset <- -1*onset;
if (end>0)				end <- -1*end;
time_scale <- time_scale[time_scale>=onset]
time_scale <- time_scale[time_scale<=end]
lst <- length(time_scale)

par(pin=c(xsize,ysize));
#ages <- vector(length=1+abs(ceiling(onset/main_time_tick))-abs(ceiling(end/main_time_tick)))
st <- abs(ceiling(onset/main_time_tick))+1	# get increment distance
#ages <- c()
#for (i in 1:length(ages))	{
#	ages[i] <- main_time_tick*(st-i)
#	ages <- c(ages,main_time_tick*(st-i))
#	}
ages1 <- seq(main_time_tick*floor(abs(onset)/main_time_tick),main_time_tick*ceiling(abs(end)/main_time_tick),by=-1*abs(main_time_tick))
med_time_tick <- main_time_tick/2
ages2 <- seq(med_time_tick*floor(abs(onset)/med_time_tick),med_time_tick*ceiling(abs(end)/med_time_tick),by=-1*abs(med_time_tick))[!seq(med_time_tick*floor(abs(onset)/med_time_tick),med_time_tick*ceiling(abs(end)/med_time_tick),by=-1*abs(med_time_tick)) %in% ages1]
minor_time_tick <- main_time_tick/10
ages3 <- seq(minor_time_tick*floor(abs(onset)/minor_time_tick),minor_time_tick*ceiling(abs(end)/minor_time_tick),by=-1*abs(minor_time_tick))
ages3 <- ages3[!ages3 %in% c(ages1,ages2)]
exy <- (mxy-mny)/25

plot(NA,type='n',axes=FALSE,main="",xlab="",ylab=ordinate,xlim=c(onset,end),ylim=c(mny,mxy))
#axis(1,at=seq(onset,end,by=abs(onset-end)),tcl=-0.0,labels=FALSE,lwd=4/3)
#axis(1,at=seq(main_time_tick*ceiling(onset/main_time_tick),end,by=main_time_tick),tcl=-0.3,labels=ages,lwd=0,lwd.ticks=4/3)
#axis(1,at=seq(main_time_tick*ceiling(onset/main_time_tick),end,by=(main_time_tick/2))[!seq(main_time_tick*ceiling(onset/main_time_tick),end,by=(main_time_tick/2))%in%seq(main_time_tick*ceiling(onset/main_time_tick),end,by=main_time_tick)],tcl=-0.2,labels=FALSE,lwd=0,lwd.ticks=4/3)
#axis(1,at=seq(10*round(onset/10),end,by=10)[!seq(10*round(onset/10),end,by=10)%in%seq(main_time_tick*ceiling(onset/main_time_tick),end,by=(main_time_tick/2))],tcl=-0.1,labels=FALSE,lwd=0,lwd.ticks=4/3)

if (hues)	{
#	per_colors <- Phanerozoic_Period_Colors()
	for (b in 1:(lst-1))	{
		per_col <- hermione_interval_color(abs(time_scale[b]),abs(time_scale[b+1]))
		rect(time_scale[b],(mny+exy),time_scale[b+1],(mny-exy),col=per_col,border=per_col)
#		rect(time_scale[b],(mny+exy),time_scale[b+1],(mny-exy),col=per_colors[b],border=per_colors[b])
		}
	}
#segments(onset,mny+exy,onset,mny-exy)
for (s in 1:lst)  segments(time_scale[s],mny+exy,time_scale[s],mny-exy)
segments(onset,mny+exy,end,mny+exy,lwd=2)
segments(onset,mny-exy,end,mny-exy,lwd=2)
axis(1,at=-1*ages1,tcl=-0.3,labels=ages1,lwd=0,lwd.ticks=4/3)
axis(1,at=-1*ages2,tcl=-0.2,labels=FALSE,lwd=0,lwd.ticks=4/3)
axis(1,at=-1*ages3,tcl=-0.1,labels=FALSE,lwd=0,lwd.ticks=4/3)

if (alt_back!=FALSE)	{
	for (b in 1:(lst-1))	{
		if (b%%2==0)	{
			rect(time_scale[b],(mny+exy),time_scale[b+1],mxy,col=alt_back,border=alt_back)
			}
		}
	}

#mid <- c()
for (i in 1:(lst-1)) {
#	mid <- c(mid,(time_scale[i]+time_scale[(i+1)])/2)
	mid <- (time_scale[i]+time_scale[(i+1)])/2
	text(mid,mny,label=stage_names[i],cex=0.9)
	}
}

Phanerozoic_Timescale_Plot_Flexible <- function(onset=-541,end=0,time_scale_to_plot,mxy,mny,use_strat_labels=TRUE,strat_names,strat_colors,plot_title="",ordinate="",abscissa="Ma",yearbreaks,xsize=6, ysize=4.285714285,hues=TRUE,colored="base",alt_back=FALSE,strat_label_size=1)	{
# UPDATED 2016-12-31
#	Smaller updates 2018-08-31
# onset: when to start x-axis (use -100 for 100 Million years ago)
# end: when to end x-axis 0
# time_scale: vector giving start of each bin, with one more than bins
# mxy: maximum value for y-axis
# mny: minimum value for y-axis
# use_strat_labels: true means that you put labels on stratigraphic unts
# strat_names: names to use for those labels
# strat_colors
# ordinate: label for Y-axis
# yearbreaks: where to put minor & major breaks; c(10,50,100) will put minor
#	medium and major breaks at 10, 50 & 100 Ma marks
# xsize: dimension of X-axis for plot for par(pin=c(xsize,ysize)) command
# ysize: dimension of Y-axis for plot for par(pin=c(xsize,ysize)) command
# hues: we want colors
# colored: where the colors should be: base or background
# alt_back: alternating white and light gray backgrouns
# use_strat_labels added 2016/09/23
if (time_scale_to_plot[1]>0)	time_scale_to_plot <-  -1*time_scale_to_plot;
if (onset>0)			onset <- -1*onset;
if (end>0)				end <- -1*end;
lst <- length(time_scale_to_plot)
first_period <- sum(time_scale_to_plot<=onset)
last_period <- sum(time_scale_to_plot<=end)
#if (is.na(match(end,time_scale)))	{
#	last_period <- match(end,time_scale)
#	}	else {
#	last_period <- match(end,time_scale)-1
#	}

draws <- length(yearbreaks)
stx <- vector(length=draws)
for (i in 1:draws)	stx[i] <- yearbreaks[i]*ceiling(onset/yearbreaks[i])

exy <- (mxy-mny)/25;
#print(c(xsize,ysize));
par(pin=c(xsize,ysize));
plot(NA,type='n',axes=FALSE,main=plot_title,xlab=abscissa,ylab=ordinate,xlim=c(onset,end),ylim=c(mny-exy,mxy))
#axis(1,at=seq(onset,end,by=abs(onset-end)),tcl=0.0,labels=FALSE,lwd=2)
if (draws==2)	{
	pts <- c(-0.15,-0.30)
	}	else if (draws==3)	{
	pts <- c(-0.1,-0.2,-0.3)
	}
for (i in 1:draws)	{
	on <- yearbreaks[i]*ceiling(onset/yearbreaks[i])
	en <- yearbreaks[i]*floor(end/yearbreaks[i])
	tcs <- seq(on,en,by=yearbreaks[i])
	if (i < draws)	{
		on2 <- yearbreaks[i+1]*ceiling(onset/yearbreaks[i+1])
		en2 <- yearbreaks[i+1]*floor(end/yearbreaks[i+1])
		tcs <- tcs[!tcs %in% seq(on2,en2,by=yearbreaks[i+1])]
		axis(1,at=tcs,tcl=pts[i],labels=FALSE,lwd=0,lwd.ticks=4/3)
		}	else	{
		axis(1,at=tcs,tcl=pts[i],labels=abs(tcs),lwd=0.0,lwd.ticks=4/3)
		}
	}

if (alt_back && colored!="backdrop")
	for (b in first_period:last_period)
		if (gtools::even(b))	rect(max(onset,time_scale_to_plot[b]),mny,min(end,time_scale_to_plot[b+1]),mxy,col="gray90",border="gray90")

if (hues)	{
	if (colored=="base")	{
		for (b in first_period:last_period)
			rect(max(onset,time_scale_to_plot[b]),mny,min(end,time_scale_to_plot[b+1]),(mny-2*exy),col=as.character(strat_colors[b]),border=as.character(strat_colors[b]))
		} else {
		for (b in first_period:last_period)
			rect(max(onset,time_scale_to_plot[b]),mny,min(end,time_scale_to_plot[b+1]),mxy,col=strat_colors[b],border=strat_colors[b])
		rect(onset,mny+exy,end,mxy,col=makeTransparent("white",50),border=makeTransparent("white",50))
		}
	}

# put in lines separating stages
for (s in first_period:last_period)  segments(time_scale_to_plot[s],mny,time_scale_to_plot[s],mny-2*exy,lwd=0.5,col="gray25")
# put in decisive final line
if (!is.na(match(end,time_scale_to_plot)))	segments(end,mny,end,mny-2*exy);
segments(onset,mny,end,mny,lwd=2)
segments(onset,mny-2*exy,end,mny-2*exy,lwd=2)

if (use_strat_labels)	{
	thin_bin <- min(abs(time_scale_to_plot[2:length(time_scale_to_plot)]-time_scale_to_plot[1:(length(time_scale_to_plot)-1)])/abs(onset-end));
	mid <- vector(length=(lst-1));
	for (i in 1:(lst-1))    mid[i] <- (max(onset,time_scale_to_plot[i])+min(end,time_scale_to_plot[(i+1)]))/2;
	dummy <- vector(length=(lst-1));
	for (i in 1:(lst-1))    dummy[i] <- mny-exy;
	#add something
	text(mid[first_period:last_period],dummy[first_period:last_period],label=strat_names[first_period:last_period],cex=min(2,thin_bin*xsize*strat_label_size));
	}
}
														
Phanerozoic_Timescale_Plot_Flexible_Ordinate <- function(onset=-541,end=0,time_scale_to_plot,mxx,mnx,use_strat_labels=TRUE,strat_names,strat_colors,plot_title="",abscissa="",ordinate="Ma",yearbreaks,ysize=6,xsize=4.285714285,hues=TRUE,colored="base",alt_back=FALSE,strat_label_size=1,stage_box_width=0.025)	{
# UPDATED 2016-12-31
#	Smaller updates 2018-08-31
#	Smaller updates 2019-09-16
#	Smaller updates 2019-10-18
# onset: when to start x-axis (use -100 for 100 Million years ago)
# end: when to end x-axis 0
# time_scale: vector giving start of each bin, with one more than bins
# mxy: maximum value for y-axis
# mny: minimum value for y-axis
# use_strat_labels: true means that you put labels on stratigraphic unts
# strat_names: names to use for those labels
# strat_colors
# ordinate: label for Y-axis
# yearbreaks: where to put minor & major breaks; c(10,50,100) will put minor
#	medium and major breaks at 10, 50 & 100 Ma marks
# xsize: dimension of X-axis for plot for par(pin=c(xsize,ysize)) command
# ysize: dimension of Y-axis for plot for par(pin=c(xsize,ysize)) command
# hues: we want colors
# colored: where the colors should be: base or background
# alt_back: alternating white and light gray backgrouns
# use_strat_labels added 2016/09/23
if (time_scale_to_plot[1]>0)	time_scale_to_plot <-  -1*time_scale_to_plot;
if (onset>0)			onset <- -1*onset;
if (end>0)				end <- -1*end;
lst <- length(time_scale_to_plot)
first_period <- sum(time_scale_to_plot<=onset)
last_period <- sum(time_scale_to_plot<=end)

draws <- length(yearbreaks)
stx <- vector(length=draws)
for (i in 1:draws)	stx[i] <- yearbreaks[i]*ceiling(onset/yearbreaks[i])

exx <- stage_box_width*(mxx-mnx);
#print(c(xsize,ysize));
par(pin=c(xsize,ysize));
#plot(NA,type='n',axes=FALSE,main=plot_title,ylab=ordinate,xlab=abscissa,ylim=c(onset,end),xlim=c(mnx-exx,mxx));
plot(NA,type='n',axes=FALSE,main=plot_title,ylab=ordinate,xlab=abscissa,ylim=c(onset,end),xlim=c(mnx-exx,mxx),xaxs="i",yaxs="i");
#axis(2,at=seq(onset,end,by=abs(onset-end)),tcl=0.0,labels=FALSE,lwd=2);
# set size of ticks
if (draws==2)	{
	pts <- c(-0.15,-0.30)
	}	else if (draws==3)	{
	pts <- c(-0.1,-0.2,-0.3)
	}
for (i in 1:draws)	{
	on <- yearbreaks[i]*ceiling(onset/yearbreaks[i]);
	en <- yearbreaks[i]*floor(end/yearbreaks[i]);
	tcs <- seq(on,en,by=yearbreaks[i]);
	if (i < draws)	{
		on2 <- yearbreaks[i+1]*ceiling(onset/yearbreaks[i+1]);
		en2 <- yearbreaks[i+1]*floor(end/yearbreaks[i+1]);
		tcs <- tcs[!tcs %in% seq(on2,en2,by=yearbreaks[i+1])]
		axis(2,at=tcs,tcl=pts[i],labels=FALSE,lwd=0,lwd.ticks=4/3);
		}	else	{
		axis(2,at=tcs,tcl=pts[i],labels=abs(tcs),lwd=0.0,lwd.ticks=4/3,las=2)
		}
	}

if (alt_back && colored!="backdrop")
	for (b in first_period:last_period)
		if (gtools::even(b))	rect(max(onset,time_scale_to_plot[b]),mnx,min(end,time_scale_to_plot[b+1]),mxx,col="gray90",border="gray90")

if (hues)	{
	if (colored=="base")	{
		for (b in first_period:last_period)
			rect(0,max(onset,time_scale_to_plot[b]),(mnx-exx),min(end,time_scale_to_plot[b+1]),col=as.character(strat_colors[b]),border=as.character(strat_colors[b]));
#			rect(mnx,max(onset,time_scale_to_plot[b]),0,min(end,time_scale_to_plot[b+1]),col=as.character(strat_colors[b]),border=as.character(strat_colors[b]));
		} else {
		for (b in first_period:last_period)
			rect(mnx,max(onset,time_scale_to_plot[b]),mxx,min(end,time_scale_to_plot[b+1]),col=strat_colors[b],border=strat_colors[b])
		rect(mnx+exx,onset,mxx,end,col=makeTransparent("white",50),border=makeTransparent("white",50))
		}
	}

# put in lines separating stages
for (s in first_period:last_period)  segments(0,time_scale_to_plot[s],mnx-exx,time_scale_to_plot[s],lwd=0.5,col="gray25");
# put in decisive final line
if (!is.na(match(end,time_scale_to_plot)))	segments(0,end,mnx-exx,end);
segments(0,onset,0,end,lwd=2);
axis(2,at=seq(onset,end,by=abs(onset-end)),tcl=0.0,labels=FALSE,lwd=2);
#points(0,-450)
#points(mnx,-450)
#points(mnx-exx,-450)
if (use_strat_labels)	{
	thin_bin <- min(abs(time_scale_to_plot[2:length(time_scale_to_plot)]-time_scale_to_plot[1:(length(time_scale_to_plot)-1)])/abs(onset-end));
	mid <- vector(length=(lst-1));
	for (i in 1:(lst-1))    mid[i] <- (max(onset,time_scale_to_plot[i])+min(end,time_scale_to_plot[(i+1)]))/2;
	dummy <- vector(length=(lst-1));
	for (i in 1:(lst-1))    dummy[i] <- (mnx-exx)/2;
	#add something
	text(dummy[first_period:last_period],mid[first_period:last_period],label=strat_names[first_period:last_period],cex=min(2,thin_bin*xsize*strat_label_size));
	}
}

#### routines to get IGS colors for geological ages
Phanerozoic_Period_Colors <- function(Carboniferous=TRUE)	{
# if Carboniferous is TRUE, then use that; otherwise, use Miss. & Penn.
if (Carboniferous)	{
	period_col <- vector(length=11)
	names(period_col) <- c("Cambrian","Ordovician","Silurian","Devonian","Carboniferous","Permian","Triassic","Jurassic","Cretaceous","Paleogene","Neogene")
	period_col[5] <- "#67A599"	# Carboniferous
	period_col[6] <- "#F04028"	# Permian
	period_col[7] <- "#812B92"	# Triassic
	period_col[8] <- "#34B2C9"	# Jurassic
	period_col[9] <- "#7FC64E"	# Cretaceous
	period_col[10] <- "#F4B470"	# Neogene
	period_col[11] <- "#FFE619"	# Paleogene
	}	else	{
	names(period_col) <- c("Cambrian","Ordovician","Silurian","Devonian","Mississippian","Pennsylvanian","Permian","Triassic","Jurassic","Cretaceous","Paleogene","Neogene")
	period_col <- vector(length=12)
	period_col[5] <- "#678F66"	# Mississippian
	period_col[6] <- "#99C2B5"	# Pennsylvanian
	period_col[7] <- "#F04028"	# Permian
	period_col[8] <- "#812B92"	# Triassic
	period_col[9] <- "#34B2C9"	# Jurassic
	period_col[10] <- "#7FC64E"	# Cretaceous
	period_col[11] <- "#F4B470"	# Neogene
	period_col[12] <- "#FFE619"	# Paleogene
	}
period_col[1] <- "#7FA056"	# Cambrian
period_col[2] <- "#009270"	# Ordovician
period_col[3] <- "#B3E1B6"	# Silurian
period_col[4] <- "#CB8C37"	# Devonian
return (period_col)
}

Phanerozoic_Interval_Colors <- function()	{
interval_colors <- c("Cenozoic","#F2F91D")
interval_colors <- rbind(interval_colors,c("Mesozoic","#67C5CA"))
interval_colors <- rbind(interval_colors,c("Paleozoic","#99C08D"))
return(interval_colors)
}

accio_time_scale_color <- function()	{
interval <- c("Eoarchean","Archean","Precambrian","Precambrian-Cambrian","Precambrian-Paleozoic","Precambrian-Phanerozoic","Paleoarchean","Mesoarchean","Kenoran","Neoarchean","Siderian","Aphebian","Paleoproterozoic","Proterozoic","Proterozoic-Cambrian","Proterozoic-Paleozoic","Rhyacian","Orosirian","Hudsonian","Statherian","Calymmian","Paleohelikian","Mesoproterozoic","Elsonian","Ectasian","Neohelikian","Stenian","Tonian","Neoproterozoic","Neoproterozoic-Cambrian","Neoproterozoic-Paleozoic","Cryogenian","Hadrynian");
interval <- c(interval,"Ediacaran","Fortunian","Terreneuvian","Begadean","Early Cambrian","Waucoban","Cambrian","Cambrian-Ordovician","Paleozoic","Phanerozoic","Stage 2","Wyattia","Fritzaspis","Stage 3","Epoch 2","Series 2","Fallotaspis","Montezuman","Nevadella","Olenellus","Dyeran","Stage 4","Middle Cambrian","Albertan","Eokochaspis nodosa","Delamaran","Amecephalus arrojosensis","Plagiura-Poliella","Wuliuan","Stage 5","Epoch 3","Series 3","Albertella","Ptychagnostus praecurrens","Topazan","Ptychagnostus gibbus","Drumian","Bolaspidella","Marjuman","Dresbachian","Late Cambrian","Guzhangian","Cedaria","Crepicephalus","Aphelaspis","Paibian","Steptoean","Furongian","Franconian","Dunderbergia","Jiangshanian","Elvinia","Taenicephalus","Sunwaptan","Pseudoyuepingia asaphoides","Trempealeauan","Ellipsocephaloides-Idahoia","Saukiella junia/Saukiella pyrene","Nelegerian","Stage 10","Saukiella serotina","Eurekia apopsis","Missisquoia","Skullrockian","Symphysurina brevispicata","Symphysurina bulbosa","Iapetognathus fluctivagus","Tremadocian","Early Ordovician","Ordovician","Ordovician-Silurian","Cordylodus angulatus","Rossodus manitouensis","Macerodus dianae","Stairsian","Acodus delatus/Oneotodus costatus","Tulean","Floian","Arenigian","Oepikodus communis","Blackhillsian","Reutterodus andinus","Microzarkodina flabellum/Tripodus laevus","Dapingian","Middle Ordovician","Whiterock","Histiodella altifrons","Histiodella sinuosa","Darriwilian","Llanvirnian","Histiodella holodentata","Chazyan","Phragmodus polonicus","Cahabagnathus friendsvillensis","Llandeilo","Cahabagnathus sweeti","Blackriverian","Sandbian","Caradocian","Late Ordovician","Plectodina aculeata","Rocklandian-Kirkfield","Rocklandian-Shermanian","Rocklandian","Mohawkian","Erismodus quadridactylus","Kirkfield","Belodina compressa","Shermanian","Phragmodus undatus","Edenian-Maysvillian","Plectodina tenuis","Katian","Belodina confluens","Edenian","Oulodus velicuspis","Maysvillian","Richmondian-Hirnantian","Oulodus robustus","Aphelognathus grandis","Richmondian","Ashgill","Aphelognathus divergens","Aphelognathus shatzeri","Hirnantian","Gamachian","Distomodus kentuckyensis","Rhuddanian","Alexandrian","Llandovery","Early Silurian","Silurian","Silurian-Devonian","Aspelunda expansa","Aeronian","Pterospathodus tenuis","Distomodus staurognathoides","Telychian","Niagaran","Pterospathodus eopennatus","Pterospathodus amorphognathoides angulatus","Pterospathodus amorphognathoides lennarti","Pterospathodus amorphognathoides lithuanicus","Pterospathodus amorphognathoides amorphognathoides","Pterospathodus pennatus procerus","Sheinwoodian","Wenlock","Kockelella ranuliformis","Ozarkodina sagitta rhenana","Kockelella walliseri","Kockelella ortus ortus","Ozarkodina sagitta sagitta","Lockportian","Homerian","Ozarkodina bohemica longa","Kockelella ortus obsidata","Kockelella crassa","Gorstian","Ludlow","Late Silurian","Kockelella variabilis variabilis","Ancoradella ploeckensis","Ludfordian","Cayugan","Polygnathoides siluricus","Ozarkodina snajdri","Ozarkodina crispa","Ozarkodina eosteinhornensis","Pridoli","Oulodus elegans detortus","Gedinnian","Lochkovian","Helderbergian","Early Devonian","Ulsterian","Devonian","Devonian-Mississippian","Caudicriodus hesperius","Caudicriodus postwoschmidti","Lanea omoalpha","Lanea eleanorae","Leanea transitans","Ancryodelloides trigonicus","Masaraella pandora morpho. Beta","Pedavis gilberti","Gondwania irregularis","Siegenian","Pragian","Gondwania kindlei","Eocostapolygnathus pireneae","Eocostapolygnathus kitabicus","Emsian","Sawkillian","Deerparkian","Eocostapolygnathus excavatus");
interval <- c(interval,"Eocostapolygnathus gronbergi","Eocostapolygnathus nothoperbonus","Polygnathus inversus","Linguipolygnathus serotinus","Polygnathus costatus patulus","Polygnathus costatus partitus","Eifelian","Southwoodian","Erian","Cazenovian","Middle Devonian","Cazenovia","Polygnathus costatus costatus","Tortodus knockelianus knockelianus","Polygnathus ensensis","Polygnathus hemiansatus","Givetian","Tioughniogan","Polygnathus varcus","Schmidtognathus hermanni","Senecan","Klapperina disparilis","Mesotaxis guanwushanensis","Fingerlakesian","Frasnian","Late Devonian","Palmatolepis transitans","Palmatolepis punctata","Palmatolepis hassi","Chemungian","Palmatolepis rhenana","Palmatolepis linguiformis","Palmatolepis triangularis","Cassadagan","Famennian","Palmatolepis crepida","Chatauquan","Palmatolepis rhomboidea","Palmatolepis marginifera","Conewangan","Palmatolepis rugosa trachytera","Palmatolepis perlobata postera","Palmatolepis gracilis expansa","Siphonodella praesulcata","Siphonodella sulcata","Kinderhookian","Tournaisian","Mississippian","Carboniferous","Siphonodella duplicata","Siphonodella sandbergi-Siphonodella belkai","Siphonodella quadruplicata-Patrognathus andersoni","Gnathodus typicus-Siphonodella isosticha","Osagean","Dollimae bouckaerti","Gnathodus semiglaber-Polygnathus communis","Gnathodus pseudosemiglaber-Scallioganthus anchoralis","Visean","Gnathodus texanus","Meramecian","Gnathodus praebillineatus","Gnathodus bilineatus","Lochriea mononodosa","Chesterian","Lochriea nodosa","Lochriea ziegleri","Serpukhovian","Namurian","Lochriea cruciformis","Gnathodus bollandensis","Gnathodus postbilineatus","Declinognathodus noduliferus","Morrowan","Bashkirian","Pennsylvanian","Pennsylvanian-Permian","Idiognathoides sinuatus","Neognathodus askynensis","Idiognathodus sinuosus","Atokan","Declinognathodus marginodosus","Neognathodus atokaensis","Declinognathodus donetzianus","Moscovian","Westphalian","Neognathodus uralicus","Streptognathodus dissectus","Desmoinian","Neoghanthodus medexultimus-Streptognathodus-concinnus","Neoghanthodus round-Streptognathodus cancellosus","Streptognathodus subexcelsus","Kasimovian","Stephanian","Missourian","Idiognathodus sagittalis","Streptognathodus cancellosus","Idiognathodus toretzianus","Streptognathodus firmus","Streptognathodus simulator","Gzhelian","Virgilian","Streptognathodus vitali","Streptognathodus virgilicus","Streptognathodus simplex-Streptognathodus bellus","Streptognathodus wabaunsensis","Streptognathodus isolatus","Asselian","Wolfcampian","Early Permian","Cisuralian","Permian","Streptognathodus sigmoidalis-Streptognathodus cristellaris","Streptognathodus constrictus-Mesogondolella belladontae","Streptognathodus fusus","Streptognathodus postfusus","Sweetognathus merrilli-Mesogondolella uralensis","Sakmarian","Sweetognathus binodosus","Sweetognathus anceps-Mesogondolella bisselli","Sweetognathus whitei","Artinskian","Sweetognathus clarki","Neostreptognathodus pequopensis","Leonardian","Neostreptognathodus pnevi","Kungurian","Neostreptognathodus prayi","Sweetognathus guizhouensis","Mesogondolella lamberti-Neostreptognathodus sulcoplicatus","Roadian","Ufimian","Jinogondolella nankingensis","Kazanian","Guadalupian","Wordian","Jinogondolella aserrata","Jinogondolella postserrata","Capitanian","Jinogondolella shannoni","Jinogondolella altudaensis","Jinogondolella prexuanhanensis","Jinogondolella granti","Clarkina postitteri hongshuiensis","Ochoan","Clarkina postitteri postitteri","Dzhulfian","Wuchiapingian","Tatarian","Late Permian","Lopingian","Clarkina dukouensis","Clarkina asymmetrica","Clarkina leveni","Clarkina guangyuanensis","Clarkina transcaucasica","Clarkina orientalis","Clarkina longicuspidata","Clarkina wangi","Changhsingian","Clarkina subarinata","Clarkina changxingensis","Clarkina deflecta-Clarkina yini","Clarkina zhejiangensis-Clarkina meishanensis","lower Otoceras boreale","Hindeodus parvus","upper Otoceras boreale","Griesbachian","Induan","Scythian","Early Triassic","Triassic","Triassic-Jurassic","Mesozoic","Ophiceras commune");
interval <- c(interval,"Isarcicella isarcica","Proptychites rosenkrantzi strigatus","Neogondollela krystyni","Sweetospathodus kummeli","Proptychites candidus","Dienerian","Neospathodus dieneri Morph 3","Vavilovites sverdrupi","Neospathodus waageni","Hedenstroemia hedenstroemi","Smithian","Olenekian","Euflemingites romunderi","Borinella buurensis-Scythogondolella milleri","Anawasatchites tardus","Bajarunia eumphala","Spathian","Neospathodus pingdingshanensis","Icriospathodus collinsoni","Olenikites pilaticus","Neospathodus triangularis","Triassospathodus sosioensis","Chiosella gondolelloides","Chiosella timorensis","Anisian","Middle Triassic","Siberlingites mulleri","Neogondolella? regalis","Lenotropites caurus","Anagymnotoceras varium","Paragondolella excelsa","Eogymnotoceras deleeni","Frechites chischa","Eoprotrachyceras matutinum","Ladinian","Budurovignathus truempyi","Tuchodiceras poseidon","Budurovignathus hungaricus","Meginoceras meginae","Macl. maclearni","Paragondolella inclinata","Frankites sutherlandi","Metapolygnathus intermedius","Daxatina canadensis","Julian","Carnian","Late Triassic","Metapolygnathus tadpole","Tachyceras desatoyense","Austrotrachyceras obesum","Sirenites nanseni","Metapolygnathus polygnathiformis","Tropites dilleri","Metapolygnathus carpathicus","Tropites welleri","Metapolygnathus nodosus","Tuvalian","Klamathites macrolobatus","Metapolygnathus primitius","Lacian","Stikinoceras kerri","Norian","Epigondolella quadrata","Malayites dawsoni","Epigondolella triangularis","Juvavites magnus","Drepanites rutherfordi","Cypriodella postera","Alaunian","Mesohimavatites columbianus","Cypriodella spiculata","Cypriodella postera","Sevatian","Cypriodella serrulata","Cypriodella bidentata","Gnomohalorites cordilleranus","Cypriodella mosheri","Chochloceras amoenum","Rhaetian","Norigondolella sp.","Choristoceras crickmayi","Misikella posternstenini","Psiloceras planorbis","Hettangian","NJ1","Early Jurassic","Jurassic","Jurassic-Cretaceous","Alsatites liasicus","Schlotheimia angulata","Arietites bucklandi","NJ2","Sinemurian","Arnioceras semicostatum","Caenisites turneri","Asteroceras obtusum","Oxynoticeras oxynotum","NJ3","Echinoceras raricostatum","Uptonia jamesoni","Pliensbachian","Tragophylloceras ibex","NJ4","Prodactylioceras davoei","Amaltheus margaritatus","NJ5","Pleuroceras spinatum","NJ6","Dactylioceras tenuicostatum","Toarcian","NJ7","Harpoceras falciferum","Hildoceras bifrons","Haugia variabilis","NJ8","Grammoceras thouarsense","Dumortieria levesquei","Pleydellia aalensis","Leioceras opalinum","Aalenian","Middle Jurassic","Ludwigia murchisonae","Brasilia bradfordensis","NJ9","Graphoceras concavum","Hyperlioceras discites","Bajocian","NJ10","Witchellia laeviuscula","Stephanoceras humphriesianum","Strenoceras niortense","Garantiana garantiana","Parkinsonia parkinsoni","Zigzagiceras zigzag","Bathonian","NJ11","Procerites progracilis","Tulites subcontractus","Morrisiceras morrisi","Procerites hodsoni","Oxycerites orbis","Clydoniceras discus","Macrocephalites herveyi","Callovian","Prolanulites koenigi","NJ12","Sigloceras calloviense","Kosmoceras jason","Erymnoceras coronatum","Peltoceras athleta","Quenstedtoceras lamberti","Quenstedtoceras mariae","Oxfordian","Late Jurassic","NJ13","Cardioceras cordatum","NJ14","Perisphinctes plicatilis","Perisphinctes pumilus","Perisphinctes cautisnigrae","NJ15","Ringsteadia pseudocordata","Pictonia baylei","Kimmeridgian","Rasenia cymodoce","Aulacostephanus mutabilis","Aulacostephanus eudoxus","Aulacostephanus autissiodorensis","Tithonian","Pectinatites elegans","Pectinatites scitulus","NJ16","Pectinatites wheatleyensis","Pectinatites hudlestoni","Pavlovia pallasioides","NJ17","Pavlovia rotunda","Virgatopavlovia fittoni","Progalbanites albani","Glaucolithites glaucolithus","Galbanites okusensis","Galbanites kerberus","Titanites anguiformes","Paracraspedites oppressus","Subcraspedites primitivus","Subcraspedites preplicomphalus","Portlandian","Berriasian","Early Cretaceous","Cretaceous","Subcraspedites lamplughi","CC1","Runctonia runctoni");
interval <- c(interval,"CC2","Hectoroceras kochi","Surites icenii","Surites stenomphalus","Peregrinoceras albidum","CC3","Valanginian","Paratollia/Platylenticeras","Polyptychites","CC4","Prodichotomites","Dichotomites","Stolcoceras tuberulatum","Eleniceras paucinodum","Endemoceras amblygonium","Endemoceras noricum","Endemoceras regale","Speetoniceras inversum","Milanowskia speetonensis","Hauterivian","Craspedodiscus gottschei","CC5","Simbirskites marginatus","Simbirskites variabilis","Paracrioceras rarocinctum","Barremian","Haplocrioceras fissicostatum","CC6","Paracrioceras elegans","Paracrioceras denckmanni","Ancyloceras inexum/S. pingue","Simanocyloceras stolleyi","Parancyloceras bidentatum/Parancyloceras scalare","CC7","Deshayesites forbesi","Aptian","Deshayesites deshayesi","Tropaeum bowerbanki","Epicheloniceras martinoides","Parahoplites nutfieldiensis","Korangan","Hypacanthoplites jacobi","Leymeriella schrammeni","Albian","CC8","Douvilleiceras mammillatum","Hoplites dentatus","Euhoplites loricatus","Urutawan","Euhoplites lautus","Diploceras cristatum","Mortoniceras pricei","Mortoniceras inflatum","Mortoniceras fallax","Motuan","CC9","Mortoniceras rostratum","Mortoniceras perinflatum","Arrhaphoceras briacensis","Cenomanian","Late Cretaceous","Cretaceous-Paleogene","Neogastroplites haasi","Mantelliceras mantelli","Ngaterian","Neogastroplites cornutus","Neogastroplites muelleri","Neogastroplites americanus","Mantelliceras dixoni","Neogastroplites maclearni","Conlinoceras tarrantense - Conlinoceras gilberti","Acanthoceras granerosense","Acanthoceras muldoonense","Acanthoceras bellense","Acanthoceras amphibolum","Pleisacanthoceras wyomingense","CC10","Dunveganoceras pondi","Dunveganoceras problematicum","Dunveganoceras albertense","Arowhanan","Dunveganoceras conditum","Sciponoceras gracile (Vasoceras diartianum)","Sciponoceras gracile (Euomphaloceras septemseriatum)","Burroceras clydense","Neocardioceras juddii","Nigericeras scotti","Watinoceras devonense","Pseudaspidoceras flexuosum","CC11","Turonian","Vascoceras birchbyi","Mammites nodosoides","Collingnoniceras woollgari","Collingnoniceras praecox","Mangaotanean","Prionocyclus hyatti","Prionocyclus macombi","Scaphites warreni","Scaphites whitfieldi","CC12","Scaphites nigricollensis","Prionocyclus germari","Scaphites mariasensis","Scaphites preventricosus","CC13","Coniacian","Scaphites ventricosus","Teratan","Scaphites depressus","CC14","Clioscaphites saxitonianus","Piripauan","CC15","Santonian","Clioscaphites vermiformis","Clioscaphites choteauensis","Desmoscaphites erdmanni","CC16","Desmoscaphites bassleri","CC17","Scaphites leei","Haumurian","CC18","Scaphites hippocrepis","Campanian","Scaphites hippocrepis II","Scaphites hippocrepis III","Baculites sp. (smooth)","Baculites sp. (weak flank ribs)","Baculites obtusus","CC19","Baculites maclearni","Baculites asperiformis","Baculites sp. (smooth)","Baculites perplexus","Baculites gregoryensis","CC20","Baculites reduncus","Baculites scotti","CC21","Didymoceras nebrascense","CC22","Didymoceras stevensoni","CC23","Exiteloceras jenneyi","Didymoceras cheyennense","Baculites compressus","Baculites cuneatus","Baculites reesidei","Baculites jenseni","Baculites eliasi","Baculites baculus","CC24","Maastrichtian","Baculites clinolobatus","Hoploscaphites birkelundae","CC25","Hoploscaphites nicolleti","Jeletzkytes nebrascensis","CC26","Puercan","NP1","Danian","Early Paleocene","Paleocene","Teurian","Early Tertiary","Paleogene","Tertiary","Cenozoic","Torrejonian","NP2","NP3","Tiffanian","NP4","Selandian","Middle Paleocene","Selandian-Thanetian","NP5","Thanetian","Late Paleocene","NP6","NP7","Clarkforkian","NP8","NP9","Wasatchian","Ypresian","Early Eocene","Eocene","Waipawan","NP10","NP11","Mangaorapan","NP12","Bridgerian","NP13","NP14","Heretaungan","Uintan","Lutetian","Middle Eocene","NP15","Porangan","NP16","Bortonian","Bartonian","NP17","Duchesnean","Jacksonian","Priabonian","Chadronian","Late Eocene","NP18","Kaiatan","Runangan","NP19-20","Whaingaroan","NP21","Orellan","Rupelian","Early Oligocene");
interval <- c(interval,"Oligocene","Middle Tertiary","NP22","NP23","Whitneyan","Geringian","NP24","Arikareean","Chattian","Late Oligocene","NP25","Duntroonian","Monroecreekian","Waitakian","Harrisonian","NN1","Aquitanian","Early Miocene","Miocene","Neogene","Late Tertiary","NN2","Otaian","Burdigalian","NN3","Altonian","Hemingfordian","NN4","Barstovian","Langhian","Middle Miocene","Clifdenian","NN5","Lillburnian","Serravallian","NN6","Clarendonian","Waiauan","NN7","Tortonian","Late Miocene","Tongaporutuan","NN8","NN9","Hemphillian","NN10","NN11","Messinian","Kapitean","NN12","Zanclean","Early Pliocene","Pliocene","Opoitian","NN13","Blancan","NN14","NN15","NN16","Waipipian","Piacenzian","Late Pliocene","Mangapanian","NN17","Gelasian","Early Pleistocene","Pleistocene","Quaternary","NN18","Nukumaruan","NN19","Calabrian","Irvingtonian","Castlecliffian","Middle Pleistocene","Ionian","NN20","Haweran","Rancholabrean","NN21","Wisconsinan","Late Pleistocene","Holocene");
interval <- c(interval,"Oandu","Rakvere","Nabala","Vormsi","Pirgu","Tremadoc");
interval <- c(interval,"Sa1","Sa2","Ka1","Ka2-3","Ka4","Stage 10");

color <- c("#DA037F","#F0047F","#F04370","#D6D6D6","#EBEBEB","#FFFFFF","#F444A9","#F768A9","#F99BC1","#F99BC1","#F74F7C","#F74370","#F74370","#F73563","#C7C7C7","#B5B5B5","#F75B89","#F76898","#F875A7","#F875A7","#FDC07A","#FDC07A","#FDB462","#F8C682","#F3CC8A","#FED99A","#FED99A","#FEBF4E","#FEB342","#A6A6A6","#B5B5B5","#FECC5C","#FECC5C");
color <- c(color,"#FED96A","#99B575","#8CB06C","#8CB06C","#8CB06C","#8CB06C","#7FA056","#409963","#99C08D","#9AD9DD","#A6BA80","#A6BA80","#A6C583","#A6C583","#99C078","#99C078","#A6C583","#A6C583","#A6C583","#99C078","#99C078","#B3CA8E","#BFD99D","#BFD99D","#B3CA8E","#BFD99D","#B3CA8E","#B3D492","#B3D492","#B3D492","#A6CF86","#A6CF86","#B3D492","#B3D492","#B3D492","#B3D492","#BFD99D","#BFD99D","#A6CF86","#D3E5B2","#B3E095","#CCDFAA","#CCDFAA","#CCDFAA","#CCEBAE","#CCEBAE","#B3E095","#B3E095","#D9EABA","#B3E095","#D9F0BB","#D9F0BB","#D9F0BB","#B3E095","#D9F0BB","#E0F0C1","#D9F0BB","#B3E095","#E6F5C9","#E6F5C9","#E6F5C9","#E6F5C9","#E6F5C9","#99C08D","#E6F5C9","#33A97E","#33A97E","#33A97E","#1A9D6F","#009270","#86CDA5","#33A97E","#33A97E","#33A97E","#33A97E","#1A9D6F","#1A9D6F","#41B087","#139B77","#41B087","#41B087","#41B087","#009270","#66C092","#4DB47E","#27A37F","#66C092","#4DB47E","#74C69C","#3AAC86","#74C69C","#4DB58D","#74C69C","#74C69C","#61BD95","#009270","#8ED195","#8CD094","#8FD297","#7FCA93","#8CD094","#93D39A","#94D49B","#91D298","#7FCA93","#8CD094","#96D59C","#8CD094","#97D59E","#8CD094","#9DD8A3","#99D69F","#99D69F","#99D69F","#9BD7A1","#99D69F","#A0D9A5","#A4DAA9","#99D69F","#99D69F","#A2D9A7","#A4DAA9","#99D69F","#7FCA93","#A6DBAB","#A6DBAB","#A6DCB5","#A6DCB5","#A6DCB5","#99D7B3","#99D7B3","#B3E1B6","#BFB777","#99D7B3","#B3E1C2","#B3E1C2","#99D7B3","#BFE6CF","#BFE6CF","#BFE6CF","#BFE6CF","#BFE6CF","#BFE6CF","#BFE6CF","#BFE6C3","#BFE6C3","#B3E1C2","#BFE6C3","#BFE6C3","#BFE6C3","#BFE6C3","#CCEBD1","#CCEBD1","#CCEBD1","#CCEBD1","#CCEBD1","#CCECDD","#CCECDD","#BFE6CF","#E6F5FF","#CCECDD","#D9F0DF","#D9F0DF","#E0F3E0","#D9F0DF","#D9F0DF","#D9F0DF","#E6F5E1","#E6F5E1","#E6F5E1","#E5B75A","#E5B75A","#E5AC4D","#E5AC4D","#E5AC4D","#CB8C37","#998E4F","#E5B75A","#E5B75A","#E5B75A","#E5B75A","#E5B75A","#E5B75A","#E5B75A","#E5B75A","#E5AC4D","#E5C468","#E5C468","#E5C468","#E5C468","#E5D075","#E5D075","#E5D075","#E5D075","#E5D075","#E5D075","#E5D075","#E5D075","#E5D075","#E5D075","#F1D576","#F1D576","#F1D576","#F1C868","#F1C868","#F1C868","#F1C868","#F1D576","#F1D576","#F1D576","#F1E185","#F1E185","#F1E185","#F1E185","#F1E185","#F1CE78","#F1E185","#CB8C37","#F1D487","#F2EDAD","#F1E19D","#F2EDAD","#F2EDAD","#F2EDAD","#F2DB97","#F2EDAD","#F2EDAD","#F2EDC5","#F2E1A6","#F2EDC5","#F2EDC5","#F2EDC5","#F2EDC5","#F2EDC5","#F2E7B6","#F2EDC5","#F2EDC5","#F2EDC5","#F2EDC5","#8CB06C","#96B46C","#8CB06C","#678F66","#67A599","#8CB06C","#8CB06C","#8CB06C","#8CB06C","#A0B76C","#8CB06C","#8CB06C","#A6B96C","#A6B96C","#A6B96C","#ABBB6B","#A6B96C","#A6B96C","#A6B96C","#B5BE6B","#A6B96C","#BFC26B","#BFC26B","#A6C093","#BFC26B","#BFC26B","#BFC26B","#99C2B5","#9FC4B7","#99C2B5","#99C2B5","#C5816F","#99C2B5","#99C2B5","#99C2B5","#ACC9BC","#99C2B5","#99C2B5","#C7CBB9","#C7CBB9","#A6C7BA","#C7CBB9","#C7CBB9","#B3CBBE","#C7CBB9","#C7CBB9","#BFD0C5","#BFD0C5","#BFD0C3","#B9CDC0","#BFD0C5","#BFD0C5","#BFD0C5","#BFD0C5","#CCD4C7","#CCD4C7","#C6D2C5","#CCD4C7","#CCD4C7","#CCD4C7","#CCD4C7","#E36350","#E36350","#E36956","#EE5845","#EF5845","#F04028","#E36350","#E36350","#E36350","#E36350","#E36F5C","#E36F5C","#E36F5C","#E36F5C","#E37B68","#E37B68","#E37B68","#E37B68","#E3816F","#E38776","#E38776","#E38776","#E38776","#E38776","#FB8069","#F57F71","#FB8069","#EE7E79","#FB745C","#FB8D76","#FB8D76","#FB9A85","#FB9A85","#FB9A85","#FB9A85","#FB9A85","#FB9A85","#FB9A85","#E27C88","#FCB4A2","#E87D80","#FCB4A2","#E27C88","#FBA794","#FBA794","#FCB4A2","#FCB4A2","#FCB4A2","#FCB4A2","#FCB4A2","#FCB4A2","#FCB4A2","#FCC0B2","#FCC0B2","#FCC0B2","#FCC0B2","#FCC0B2","#FCC0B2","#FCC0B2","#A4469F","#A4469F","#DC7B90","#A4469F","#D57998","#983999","#812B92","#5B6FAE","#67C5CA","#A4469F","#A4469F","#A4469F","#983999","#A4469F","#A4469F");
color <- c(color,"#CF78A0","#A4469F","#A4469F","#B051A5","#B051A5","#C977A7","#B051A5","#B051A5","#B051A5","#B051A5","#B051A5","#C276AF","#B051A5","#B051A5","#B051A5","#B051A5","#B051A5","#B051A5","#BC75B7","#BC75B7","#B168B1","#BC75B7","#BC75B7","#BC75B7","#BC75B7","#B168B1","#BC75B7","#BC75B7","#C983BF","#C983BF","#C983BF","#C983BF","#C983BF","#C983BF","#C983BF","#C983BF","#C983BF","#C99BCB","#C99BCB","#C99BCB","#C99BCB","#BD8CC3","#C99BCB","#C99BCB","#C99BCB","#C99BCB","#C99BCB","#C99BCB","#C99BCB","#C99BCB","#C99BCB","#C99BCB","#C99BCB","#BD8CC3","#D6AAD3","#D6AAD3","#D6AAD3","#D6AAD3","#D6AAD3","#D6AAD3","#D6AAD3","#D6AAD3","#D6AAD3","#D6AAD3","#D6AAD3","#D6AAD3","#D6AAD3","#D6AAD3","#D6AAD3","#D6AAD3","#D6AAD3","#E3B9DB","#E3B9DB","#E3B9DB","#E3B9DB","#E3B9DB","#E3B9DB","#4EB3D3","#4EB3D3","#56B6D5","#42AED0","#34B2C9","#5ABC8C","#4EB3D3","#4EB3D3","#67BCD8","#5FB9D6","#67BCD8","#67BCD8","#67BCD8","#67BCD8","#67BCD8","#74C1DB","#67BCD8","#80C5DD","#80C5DD","#80C5DD","#88C8DF","#80C5DD","#80C5DD","#91CBE1","#80C5DD","#99D1E2","#99CEE3","#99CEE3","#9AD4E0","#99CEE3","#99CEE3","#99CEE3","#9AD6DF","#99CEE3","#99CEE3","#9AD9DD","#9AD9DD","#9AD9DD","#80CFD8","#9AD9DD","#9AD9DD","#A9DEE1","#9AD9DD","#A6DDE0","#A6DDE0","#ADE0E2","#A6DDE0","#A6DDE0","#A6DDE0","#A6DDE0","#A6DDE0","#B3E2E3","#B3E2E3","#B0E1E2","#B3E2E3","#B3E2E3","#B3E2E3","#B3E2E3","#B3E2E3","#B3E2E3","#BFE7E5","#BFE7E5","#BFE7E5","#B9E5E4","#BFE7E5","#BFE7E5","#BFE7E5","#BFE7E5","#BFE7E5","#BFE7F1","#BFE7F1","#B3E3EE","#BFE7E9","#BFE7F1","#BFE7ED","#BFE7F1","#BFE7F1","#BFE7F1","#C6EAF3","#BFE7F1","#CCECF4","#CCECF4","#CCECF4","#CCECF4","#CCECF4","#CCECF4","#D9F1F7","#D9F1F7","#D9F1F7","#D3EFF6","#D9F1F7","#D9F1F7","#D9F1F7","#D9F1EE","#D9F1F7","#D9F1F7","#D9F1F7","#D9F1F7","#D9F1F7","#D9F1F7","#D9F1F7","#D9F1F7","#D9F1F7","#D9F1F7","#8CCD60","#8CCD60","#8CCD57","#7FC64E","#8CCD60","#90CF63","#8CCD60","#95D167","#8CCD60","#8CCD60","#8CCD60","#99D36A","#9DD56E","#99D36A","#99D36A","#99D36A","#A2D771","#99D36A","#99D36A","#99D36A","#99D36A","#99D36A","#99D36A","#99D36A","#99D36A","#8CCD57","#A6D975","#A6D975","#ADDC7A","#A6D975","#A6D975","#B3DF7F","#B3DF7F","#B3DF7F","#B7E183","#B3DF7F","#B3DF7F","#B3DF7F","#B3DF7F","#B3DF7F","#BBE286","#BFE48A","#BFE48A","#BFE48A","#BFE48A","#BFE48A","#BFE48A","#70C189","#BFE48A","#CCEA97","#CCEA97","#C6E791","#CCEA97","#CCEA97","#CCEA97","#77C48E","#CCEA97","#CCEA97","#CCEA97","#CCEA97","#CCEA97","#7DC693","#C0E475","#CCEA97","#CCEA97","#B3DE55","#B3DE53","#A6D84A","#D2C055","#B3DE53","#B3DE53","#84C998","#B3DE53","#B3DE53","#B3DE53","#B3DE53","#B3DE53","#B3DE53","#B3DE53","#B3DE53","#B3DE53","#B3DE53","#B3DE53","#B7E056","#B3DE53","#B3DE53","#B3DE53","#91CEA3","#B3DE53","#B3DE53","#B3DE53","#B3DE53","#B3DE53","#B3DE53","#B3DE53","#BFE35D","#BBE15A","#BFE35D","#BFE35D","#BFE35D","#BFE35D","#BFE35D","#98D1A8","#BFE35D","#BFE35D","#BFE35D","#BFE35D","#C3E561","#BFE35D","#BFE35D","#BFE35D","#A6D84A","#C8E764","#CCE968","#CCE968","#9ED3AE","#CCE968","#D0EB6C","#A6D84A","#AEDABA","#D5ED70","#D9EF74","#D9EF74","#D9EF74","#D9EF74","#DCF077","#D9EF74","#E0F27A","#D9EF74","#B8DEC2","#E3F37C","#A6D84A","#E6F47F","#E6F47F","#E6F47F","#E6F47F","#E6F47F","#E6F47F","#E8F581","#E6F47F","#E6F47F","#E6F47F","#E6F47F","#E6F47F","#E9F683","#E6F47F","#E6F47F","#EBF785","#E6F47F","#EDF786","#E6F47F","#EFF888","#E6F47F","#E6F47F","#E6F47F","#E6F47F","#E6F47F","#E6F47F","#E6F47F","#A6D84A","#F0F98A","#F2FA8C","#F2FA8C","#F2FA8C","#F2FA81","#F2FA8C","#F2FA8C","#F2FA76","#FDB469","#FDB663","#FDB462","#FDB462","#FDA75F","#F4AA70","#FDB46C","#FD9A52","#F4B470","#F2F91D","#FEBA64","#FDB863","#FEBB64","#FEBF6A","#FEBD64","#FEBF65","#FEBF65","#FEBF6A","#FEBF6A","#FDBF6F","#FDBF6F","#FDBA70","#FDB571","#FDB371","#FCB171","#FCAC72","#FCAB78","#FCA773","#FCA773","#FDB46C","#EFA76E","#FCA976","#FCAB78","#EAA46D","#FCAE7B","#FCB07D","#FCB07D","#FCB280","#E4A06B","#FCB887","#FCB482","#FCB482","#FCB887","#DF9D69","#FDBC8C","#E3A570","#FDC091","#FDC799","#FDBC8C","#FDCDA1","#FDCDA1","#FDC799","#FED99A","#FDD09F","#EBB07A","#F4BA83");
color <- c(color,"#FED39E","#FCCA96","#FED69C","#FED39E","#FED99A","#FED99A","#FDC07A","#FDC07A","#FEDC9E","#FEE0A2","#FEE39F","#FFECA5","#FEE3A6","#FEE9B2","#FEE6AA","#FEE6AA","#FEE6AA","#FDCFA0","#FEEE82","#FDD5AA","#FFF75B","#FFF75B","#FFFF33","#FFFF33","#FFFF00","#FFE619","#FFFF00","#FFFF3A","#FEDAB4","#FFFF41","#FFFF45","#FDD6A7","#FFFF3A","#FFFF49","#FFFF47","#FFFF4D","#FFFF4D","#FDCD90","#FFFF53","#FDCA7B","#FFFF59","#FFFF5D","#FFFF60","#FDD07C","#FFFF62","#FFFF66","#FFFF66","#FED67E","#FFFF69","#FFFF6B","#FFFF6D","#FFFF6E","#FFFF70","#FFFF73","#FEE582","#FFFF93","#FFFFB3","#FFFFB3","#FFFF99","#FFF78A","#FFFFB5","#FFFFB9","#FFFFB8","#FFFFBA","#FFFFBD","#FFF88E","#FFFFBF","#FFFFBF","#FFF891","#FFF6B9","#FFEDB3","#FFEDB3","#FFF2AE","#F9F97F","#FFEFBA","#FFF994","#FFF0C0","#FFF2C7","#FFF0BD","#FFF997","#FFF2C7","#FFF2D3","#FFF2D8","#FFFA9B","#FFF2E1","#FFF2DC","#FFF2D3","#FFF2C7","#FEF2E0");
color <- c(color,"#7FCA93","#97D59E","#99D69F","#99D69F","#A2D9A7","#33A97E");
color <- c(color,"#8ED195","#96D59C","#99D69F","#A0D9A5","#A2D9A7","#E6F5C9");

onset <- c(4000,4000,4000,4000,4000,4000,3600,3200,2800,2800,2500,2500,2500,2500,2500,2500,2300,2050,1800,1800,1600,1600,1600,1400,1400,1300,1200,1000,1000,1000,1000,850,850);
onset <- c(onset,635,541,541,541,541,541,541,541,541,541,529,524,521,521,521,521,520,520,518.5,515.5,515.5,514,513,513,511.1,511,510.5,509,509,509,509,509,507.3,506.5,506.5,505.25,504.5,504.5,504.5,501,501,500.5,499.5,498.5,497,497,497,497,496.8,495.2,494,493.9,493,493,493,492.5,492.3,490.9,489.5,489.5,487.8,487.1,486.7,486.7,486.2,485.4,485.4,485.4,485.4,485.4,485.4,484.4,480.3,480,480,479.8,479.5,477.7,477.7,477.5,475,473.7,470.3,470.0,470.0,470.0,468.7,467.4,467.3,465.5,464.4,464,463.7,462.4,461.4,461.2,458.4,458.4,458.4,458.4,456.6,456.5,456.5,456.5,456.5,455.8,455.5,454.7,454.5,453.7,453.5,453,453,452.3,451,450,450,450,449.7,449,449,449,447,445.5,445.2,445.2,443.8,443.8,443.8,443.8,443.8,443.8,443.8,441.1,440.8,440.5,439.4,438.5,438.5,437.5,436.4,435.8,435.5,435.2,433.4,433.4,433.4,432.6,432.3,431.8,431.2,430.5,430.5,430.5,429.4,428.3,427.4,427.4,427.4,427.4,427,425.6,425.6,425.6,425.2,424.8,423.6,423,423,421.2,419.2,419.2,419.2,419.2,419.2,419.2,419.2,419,417.2,415.2,414.5,414,413.5,412.5,412,411,410.8,410.8,410,408.5,407.6,407.6,407.6,407.6,406,404.2,403.5,401,397.5,395,393.3,393.3,393.3,393.3,393.3,393.3,393.3,391.5,389.2,388,387.7,387.7,387.7,387.5,385,385,384,383,382.7,382.7,382.7,381.2,380.2,379,379,376.5,373,372.2,372.2,372.2,370,370,367.5,366.5,365,364,363,362.2,360,358.9,358.9,358.9,358.9,358.9,358.4,357.4,355,353,351.5,349,348,346.7,346.7,344,343.5,340,336.8,335,335,332,330.9,330.9,330.9,329.5,327.2,325.5,323.2,323.2,323.2,323.2,323.2,322.5,320.8,320,319,318.6,317,315.2,315.2,315,314.5,312.5,312,310.8,309,307,307,307,306,305.7,305.2,304.5,304,303.7,303.7,303.7,302.7,302,301,300,298.9,298.9,298.9,298.9,298.9,298.9,298,297.5,297,296.5,295.5,295.5,294.1,293.5,290.1,290.1,285,282,282,279.3,279.3,275.5,274.7,273.7,272.3,272.3,272.3,272.3,272.3,268.8,268.8,265.1,265.1,264.5,264,263.7,261.8,260.3,259.9,259.9,259.9,259.9,259.9,259.9,259.9,259,258.5,257.7,257.2,256.5,255.3,254.5,254.14,254.14,253.5,253.1,252.8,252.5,252.4,252.17,252.17,252.17,252.17,252.17,252.17,252.17,252.17,252.17,251.9,251.9,251.8,251.7,251.6,251.6,251.6,251.4,251.3,251.2,251.2,251.2,251.2,249.3,249.3,248.55,248.5,248.5,248.4,248.3,248,247.6,247.4,247.3,247.2,247.2,247.2,246.8,246.7,246.5,246.3,244.3,244,243.5,242,242,241,240.4,240.3,239,238.3,238,238,237,237,237,237,237,236.5,236.3,234.6,234,233.5,233.5,233,232.6,231.4,230.5,229.7,229.5,228.4,228,228,226.5,224.5,221.4,218.3,217.5,217.4,217.4,217,216.9,216.5,215.3,214.8,214,214,208.5,208.5,208.5,202.3,202.2,201.9,201.3,201.3,201.3,201.3,201.3,201.3,201.1,200.2,199.3,199.3,199.3,197.8,196.4,195.4,193.8,193.3,192.8,190.8,190.8,189.5,189,188.5,187.5,185.5,184.2,182.8,182.7,182.7,181,180.7,180.4,178.4,177.2,176.4,175,174.4,174.1,174.1,174.1,172.2,171.3,171.2,170.9,170.3,170.3,170.3,170,169.8,169.5,169,168.6,168.3,168.3,168,167.3,167,166.8,166.6,166.4,166.2,166.1,166.1,165.5,165.4,165,164.6,164.5,164.2,163.8,163.5,163.5,163.5,163.4,161.4,161.4,160.8,160.1,159.7,159.7,159,157.3,157.3,156,154.7,153.7,152.4,152.1,152,151.5,151.2,151,150.7,150.2,149.2,149,148.3,148,147.6,147.4,147.2,147,146.7,146.4,146,146,145,145,145,144.6,144,144,142.4,141.6,141,140.2,139.8,139.8,139.8,138.5,137.6,137.4,136.8,135.6,135,134.3,133.95,133.85,133.5,133.4,133,132.9,132.8,132.7,132.5,132.3,129.4,129.4,129,129,128.8,128.2,127.8,127.3,125.5,125.4,125,125,124.5,124,123,118,117.5,116.8,113,113,112.2,111.3,110.3,109.5,108.4,108.3,107.5,107,104.9,104.1,103.3,103,101.7,101.3,100.8,100.5,100.5,100.5,100.25,100.2,100.2,99.81,99.17,98.75,98.5,98.19,97.26,96.24,96.08,95.98,95.9,95.81,95.7,95.67,95.47,95.24,95.2,95.01,94.78,94.57,94.39,94.27,94.15,93.98,93.9,93.9,93.9,93.55,93.45,93.35,92.9,92.1,92.08,91.6,91.41,91.34,91,90.65,90.24,89.98,89.87,89.8,89.8,89.77,89.1,88.77,88,87.86,86.5,86.3,86.3,86.26,85.56,85.23,85,84.52,84.2,84.08,84,83.8,83.64,83.6,82.7,82,81.53,81.28,81.13,81,80.97,80.67,80.21,79.64,79.01,79,78.34,77.63,77.6,76.94,76.8,76.27,76,75.64,75.08,74.6,74.21)
onset <- c(onset,73.91,73.63,73.27,72.74,72.1,72.1,72.05,70.44,70.2,69.91,69.3,67.8,66,66,66,66,66,66,66,66,66,66,64.75,64.5,63.8,62.25,62.2,61.6,61.6,61.6,59.7,59.2,59.2,58.4,57.5,57.5,57.3,56.2,56,56,56,56,55.5,55,53.61,53,52.85,52,50.6,49.7,49.5,47.9,47.8,47.8,47.3,46.2,43.4,43,41.3,41.3,39.9,38,37.8,37.8,37.8,37,37,36,36,34.3,34.2,33.9,33.9,33.9,33.9,33.9,32.8,32.3,32.1,30.8,29.9,29.8,28.1,28.1,27.5,27.3,26.3,25.2,24.8,23.9,23.03,23.03,23.03,23.03,23.03,22.2,21.7,20.44,19,19,18.6,18.3,16.3,15.97,15.97,15.9,15.6,15.1,13.82,13.6,13.6,12.7,11.8,11.62,11.62,10.92,10.9,10.7,10.3,9.4,8.6,7.246,6.5,5.6,5.333,5.333,5.333,5.28,5,4.9,4.15,4,3.75,3.6,3.6,3.6,3,2.6,2.588,2.588,2.588,2.588,2.5,2.4,2,1.806,1.806,1.63,0.781,0.781,0.5,0.34,0.3,0.3,0.15,0.126,0.0117)
onset <- c(onset,453.0,452.0,450.4,449.4,448.5,485.4)
onset <- c(onset,458.4,456.5,453.0,450.3,448.3,489.5)

end <- c(3600,2500,541,485.4,252.17,0,3200,2800,2500,2500,2300,1600,1600,541,485.4,252.17,2050,1800,1600,1600,1400,1300,1000,1300,1200,850,1000,850,541,485.4,252.17,635)
end <- c(end,541,541,529,521,520,513,513,485.4,443.8,252.17,0,521,521,520,514,509,509,518.5,515.5,515.5,511.1,511,509,501,501,510.5,506.5,509,507.3,504.5,504.5,497,497,506.5,505.25,504.5,504.5,500.5,499.5,497,496.8,485.4,497,498.5,497,495.2,494,493,485.4,492.5,493.9,489.5,493,492.3,486.7,488.27,485.4,490.9,487.8,485.4,485.4,487.1,486.7,486.2,480,485.4,484.4,484.4,477.7,470.0,443.8,419.2,480.3,480,479.8,479.5,477.5,475,470.0,465.5,473.7,470.0,470.3,468.7,467.3,458.4,456.5,467.4,464.4,458.4,461.9,463.7,458.4,462.4,461.2,458.4,456.6,456.5,453,449,443.8,455.8,457,453.5,455.5,451,454.7,454.5,453.7,451,453,449,452.3,445.2,450,450,449.7,449,443.8,449,447,445.2,443.8,445.5,443.8,443.8,443.8,441.1,440.8,438.5,433.4,427.4,419.2,358.9,440.5,438.5,439.4,437.5,433.4,425.6,436.4,435.8,435.5,435.2,433.4,432.6,430.5,427.4,432.3,431.8,431.2,430.5,429.4,427.4,427.4,428.3,427.4,427,425.6,423,419.2,425.6,425.2,423,419.2,424.8,423.6,423,421.2,419.2,419,410.8,410.8,393.3,393.3,393.3,358.9,323.2,417.2,415.2,414.5,414,413.5,412.5,412,411,410,407.6,407.6,408.5,407.6,406,393.3,393.3,393.3,404.2,403.5,401,397.5,395,393.3,391.5,387.7,387.7,385,382.7,382.7,382.7,389.2,388,387.7,387.5,382.7,382.7,385,384,370,383,381.2,379,372.2,358.9,380.2,379,376.5,372.2,373,372.2,370,365,358.9,367.5,358.9,366.5,364,358.9,363,362.2,360,358.9,358.4,351.5,346.7,323.2,298.9,357.4,355,353,349,343.5,348,346.7,344,330.9,340,335,336.8,335,332,323.2,330.9,329.5,323.2,315,327.2,325.5,323.2,322.5,319,315.2,298.9,252.17,320.8,320,318.6,312,317,315.2,314.5,307,307,312.5,310.8,306,309,307,305.7,303.7,298.9,303.7,305.2,304.5,304,303.7,302.7,298.9,298.9,302,301,300,298.9,298,295.5,282,272.3,272.3,252.17,297.5,297,296.5,295.5,294.1,290.1,293.5,290.1,285,279.3,282,279.3,272.3,275.5,272.3,274.7,273.7,272.3,268.8,268.8,268.8,259.9,259.9,265.1,265.1,264.5,259.9,264,263.7,261.8,260.3,259.9,259,259,254.14,254.14,252.17,252.17,252.17,258.5,257.7,257.2,256.5,255.3,254.5,254.14,253.5,252.17,253.1,252.8,252.5,252.17,252.17,251.9,251.9,251.6,251.2,247.2,247.2,201.3,145,66,251.8,251.7,251.6,251.6,251.4,251.3,251.2,251.2,251.2,249.3,249.3,248.5,247.2,248.55,248.4,248.5,248,247.2,248.3,247.6,247.7,247.4,247.3,247.2,246.7,242,237,246.5,244.3,246.3,244,241,243.5,242,240.4,237,240.3,239,238,238.3,238,237,237,236.5,236.3,230.5,228,201.3,233.5,234.6,234,233.5,233,232.6,231.4,229.7,229.5,228,228,226.5,217.4,224.5,208.5,221.4,218.3,217.4,217.5,217,216.9,215.3,214,216.5,214.8,208.5,214,208.5,208.5,202.3,202.2,201.3,201.9,201.3,201.3,201.1,199.3,199.3,174.1,145,66,200.2,199.3,197.8,193.3,190.8,196.4,195.4,193.8,192.8,189,190.8,189.5,182.7,188.5,185.5,187.5,184.2,182.8,182.7,181,180.7,174.1,177.2,180.4,178.4,176.4,171.2,175,174.4,174.1,172.2,170.3,163.5,171.3,170.9,170.3,170.3,170,168.3,168,169.8,169.5,169,168.6,168.3,167.3,166.1,165.4,167,166.8,166.6,166.4,166.2,166.1,165.5,163.5,165,163.4,164.6,164.5,164.2,163.8,163.5,161.4,157.3,145,161.4,160.8,159.7,160.1,159.7,159,151.2,157.3,156,152.1,154.7,153.7,152.4,152,145,151.5,151,149.2,150.7,150.2,149,144,148.3,148,147.6,147.4,147.2,147,146.7,146.4,146,144.6,142,139.8,100.5,66,144,142.4,141.6,139.8,141,140.2,139.8,138.5,137.4,132.9,137.6,136.8,132.7,135.6,135,134.3,133.95,133.85,133.5,133.4,133,132.8,129.4,132.5,129,132.3,129.4,129,125,128.8,125.4,128.2,127.8,127.3,125.5,125,112.2,124.5,113,124,123,118,116.8,108.4,113,111.3,100.5,103,110.3,109.5,108.3,103.3,107.5,107,104.9,104.1,101.7,100.2,95.7,101.3,100.8,100.2,93.9,66,56,99.81,98.5,95.2,99.17,98.75,98.19,97.5,97.76,96.24,96.08,95.98,95.9,95.81,95.67,93.9,95.47,95.24,95.01,92.1,94.78,94.57,94.39,94.27,94.15,93.98,93.9,93.55,91,89.8,93.45,93.35,92.9,92.08,89.1,91.6,91.41,91.34,90.65,89.8,90.24,89.98,89.87,89.77,88,86.3,88.77,86.5,87.86,86.3,86.26,84,85,83.6,85.56,85.23,84.52,84.2,84.08,83.8,83.64,66,81,82.7,72.1,82,81.53,81.28,81.13,80.97,79,80.67,80.21,79.64,79.01,78.34,77.6,77.63,76.94,76.8,76.27,76,75.64,72.1,75.08,74.6,74.21,73.91,73.63,73.27);
end <- c(end,72.74,72.05,70.2,66,70.44,69.91,67.8,69.3,68.69,66,64.75,64.5,61.6,61.6,56,55.5,33.9,23.03,2.588,0,62.25,63.8,62.2,57.5,59.7,59.2,59.2,56,58.4,56,56,57.5,57.3,56,56.2,55,52,47.8,47.8,33.9,53,53.61,52.85,49.5,50.6,47.9,49.7,47.3,46.2,39.9,41.3,38,43.4,43,41.3,37,37.8,37,37.8,33.9,33.9,33.9,33.9,36,36,34.3,34.2,27.3,32.8,32.1,28.1,28.1,23.03,23.03,32.3,29.9,29.8,26.3,27.5,18.6,23.03,23.03,23.9,25.2,24.8,21.7,20.6,22.2,20.44,15.97,5.333,2.588,2.588,19,19,15.97,18.3,15.9,16.3,15.6,13.6,13.82,11.62,15.1,13.6,12.7,11.62,11.8,10.3,10.92,10.9,7.246,5.333,6.5,10.7,9.4,4.9,8.6,5.6,5.333,5.28,5,3.6,3.6,2.588,3.6,4.15,1.806,4,3.75,2.6,3,2.588,2.588,2.4,2.5,1.806,1.806,0.0117,0,2,1.63,0.5,0.781,0.3,0.34,0.126,0.0117,0.3,0,0.0114,0,0.05,0.0117,0);
end <- c(end,452.0,450.4,449.4,448.5,445.2,477.7);
end <- c(end,456.5,453.0,450.3,448.3,445.2,485.4);

output <- data.frame(cbind(interval,onset,end,color));
return(output)
}

accio_stage_colors <- function()	{
interval <- c("Ediacaran","Fortunian","Stage 2","Stage 3","Stage 4","Wuliuan","Drumian","Guzhangian","Paibian","Jiangshanian","Stage 10","Tremadoc","Floian","Dapingian","Darriwilian","Sandbian","Katian","Hirnantian","Rhuddanian","Aeronian","Telychian","Sheinwoodian","Homerian","Gorstian","Ludfordian","Pridoli","Lochkovian","Pragian","Emsian","Eifelian","Givetian","Frasnian","Famennian","Tournaisian","Visean","Serpukhovian","Bashkirian","Asselian","Sakmarian","Artinskian","Kungurian","Roadian","Wordian","Capitanian","Wuchiapingian","Changhsingian");
onset <- c(635.0,541.0,529.0,521.0,514.0,509.0,504.5,500.5,497.0,494.0,489.5,485.4,477.7,470.0,467.3,458.4,453.0,445.2,443.8,440.8,438.5,433.4,430.5,427.4,425.6,423.0,419.2,410.8,407.6,393.3,387.7,382.7,372.2,358.9,346.7,330.9,323.2,298.9,295.5,290.1,279.3,272.3,268.8,265.1,259.9,254.1);
end <- c(542.0,529.0,521.0,514.0,509.0,504.5,500.5,497.0,494.0,489.5,485.4,477.7,470.0,467.3,458.4,453.0,445.2,443.8,440.8,438.5,433.4,430.5,427.4,425.6,423.0,419.2,410.8,407.6,393.3,387.7,382.7,372.2,358.9,346.7,330.9,323.2,315.5,295.5,290.1,279.3,272.3,268.8,265.1,259.9,254.1,252.2);
color <- c("#FED96A","#99B575","#A6BA80","#A6C583","#B3CA8E","#B3D492","#BFD99D","#CCDFAA","#CCEBAE","#D9F0BB","#E6F5C9","#33A97E","#41B087","#66C092","#74C69C","#8CD094","#99D69F","#A6DBAB","#A6DCB5","#B3E1C2","#BFE6CF","#BFE6C3","#CCEBD1","#CCECDD","#D9F0DF","#E6F5E1","#E5B75A","#E5C468","#E5D075","#F1D576","#F1E185","#F2EDAD","#F2EDC5","#8CB06C","#A6B96C","#BFC26B","#99C2B5","#E36350","#E36F5C","#E37B68","#E38776","#FB8069","#FB8D76","#FB9A85","#FCB4A2","#FCC0B2");
return(data.frame(cbind(interval,onset,end,color)));
}

hermione_interval_color <- function(onset,end)	{
color_scale <- accio_time_scale_color();
onsets <- as.numeric(as.character(color_scale$onset));
ends <- as.numeric(as.character(color_scale$end));
start_diff <- abs(abs(onsets)-abs(onset))
end_diff <- abs(abs(ends)-abs(end))
ttl_diff <- start_diff+end_diff
color <- color_scale$color[match(min(ttl_diff),ttl_diff)]
return(color)
}

hermione_stage_color_given_age <- function(ma)	{
get_colors <- accio_time_scale_color()
after <- as.numeric(as.character(get_colors$onset))-abs(ma);
before <- abs(ma)-as.numeric(as.character(get_colors$end));
get_colors <- cbind(get_colors,after,before);
get_colors <- subset(get_colors,get_colors$after>=0);
get_colors <- subset(get_colors,get_colors$before>=0);
span <- get_colors$before+get_colors$after;
best_bet <- match(min(span),span)
#ma_color <- c(xx[best_bet,4],xx[best_bet,1])
ma_color <- get_colors$color[best_bet];
return(ma_color)
}

hermione_stage_color_given_stage_onset_and_end <- function(age_range)	{
print(age_range);
get_colors <- accio_time_scale_color()
after <- as.numeric(as.character(get_colors$onset))-abs(age_range[1]);
before <- abs(age_range[2])-as.numeric(as.character(get_colors$end));
get_colors <- cbind(get_colors,after,before);
get_colors <- subset(get_colors,get_colors$after>=0);
get_colors <- subset(get_colors,get_colors$before>=0);
span <- get_colors$before+get_colors$after;
best_bet <- match(min(span),span)
#ma_color <- c(xx[best_bet,4],xx[best_bet,1])
ma_color <- get_colors$color[best_bet];
return(ma_color)
}

hermione_stage_color_given_strat_unit <- function(strat_unit)	{
get_colors <- accio_time_scale_color()
su <- match(strat_unit,get_colors$interval)
return(get_colors[su,4])
}

pentagon_symbol <- function(x,y,abc,ord,size)	{
pent <- matrix(0,2,5)
an <- 18
for (r in 1:5)  {
	#	pent[1,r] <- x+cos(pi*an/180)*(abc*size)
	#	pent[2,r] <- y+sin(pi*an/180)*(ord*size)
	pent[1,r] <- cos(pi*an/180)
	pent[2,r] <- sin(pi*an/180)
	an <- an+72
	}
  # center and rescale
pent <- size*pent
pent[1,] <- abc*pent[1,]
pent[2,] <- ord*pent[2,]
pent[1,] <- x+pent[1,]
pent[2,] <- y+pent[2,]	
  
return(pent)
}

bug_symbol <- function(x,y,abc,ord,size)	{
bug <- matrix(0,2,6)
bug[1,1] <- x-(0.75*abc*size*(sqrt(2)/2))
bug[2,1] <- y+(0.75*ord*size*(sqrt(2)/2))
bug[1,2] <- x
bug[2,2] <- y+(ord*size)
bug[1,3] <- x+(0.75*abc*size)*(sqrt(2)/2)
bug[2,3] <- y+(0.75*ord*size*(sqrt(2)/2))
bug[1,4] <- x+(0.75*abc*size)*(sqrt(2)/2)
bug[2,4] <- y-(0.75*ord*size*(sqrt(2)/2))
bug[1,5] <- x
bug[2,5] <- y-(ord*size)
bug[1,6] <- x-(0.75*abc*size*(sqrt(2)/2))
bug[2,6] <- y-(0.75*ord*size*(sqrt(2)/2))
return(bug)
}

star_symbol <- function(x,y,abc,ord,size)	{
# x: x-coordinate
# y: y-coordinate
# abc: relative size of X-axis
# ord: relative size of Y-axis
# size: size of shell
star <- matrix(0,2,10)
an <- 18
for (r in 1:10)  {
	if (r%%2==1)	{
      #		star[1,r] <- x+cos(pi*an/180)*(abc*size)
      #		star[2,r] <- y+sin(pi*an/180)*(ord*size)
		star[1,r] <- cos(pi*an/180)/2
		star[2,r] <- sin(pi*an/180)/2
    	}
    if (r%%2==0)	{
      #		star[1,r] <- x+cos(pi*an/180)*(size*abc/2)		
      #		star[2,r] <- y+sin(pi*an/180)*(size*ord/2)		
		star[1,r] <- cos(pi*an/180)/4		
		star[2,r] <- sin(pi*an/180)/4
    	}
	an <- an+36
	}
 # rescale
star <- size*star
star[1,] <- abc*star[1,]
star[2,] <- ord*star[2,]
star[1,] <- x+star[1,]
star[2,] <- y+star[2,]
return(star)
}

flower_symbol <- function(x,y,abc,ord,size)	{
flower <- matrix(0,2,12)
an <- 0
for (r in 1:12)	{
	if (r%%2==1)	{
		flower[1,r] <- cos(pi*an/180)/2
		flower[2,r] <- sin(pi*an/180)/2
		}
	if (r%%2==0)	{
		flower[1,r] <- cos(pi*an/180)/4		
		flower[2,r] <- sin(pi*an/180)/4
	  }
	an <- an+30
	}
  # rescale
flower <- size*flower
flower[1,] <- abc*flower[1,]
flower[2,] <- ord*flower[2,]
flower[1,] <- x+flower[1,]
flower[2,] <- y+flower[2,]
return(flower)
}

mollusc_symbol <- function(x,y,abc,ord,size,whorls,W)	{
# x: x-coordinate
# y: y-coordinate
# abc: relative size of X-axis
# ord: relative size of Y-axis
# size: size of shell
# whorls: number of whorls to draw
# W: Raup's W.
stops <- (72*whorls)+1+(72*(whorls-1)+1)
snail <- matrix(0,2,stops)
  #r <- (size*abc)*/(W^whorls)
r <- 1/(W^whorls)
  #snail[1,1] <- x-r
snail[1,1] <- -r
  #snail[2,1] <- y
snail[2,1] <- 0
ang <- 5
for (i in 2:((72*whorls)+1))	{
	ri <- r*W^(5*(i-1)/360)
    #	snail[1,i] <- x+ri*cos(pi*((180-ang)/180))
    #	snail[2,i] <- y-(ord/abc)*ri*sin(pi*((180-ang)/180))
	snail[1,i] <- ri*cos(pi*((180-ang)/180))
	snail[2,i] <- -ri*sin(pi*((180-ang)/180))
	ang <- ang+5
	}
a <- 72*(whorls-1)+1
i <- (72*whorls)+2
for (i in ((72*whorls)+2):stops)	{
	snail[1,i] <- snail[1,a]
	snail[2,i] <- snail[2,a]
	a <- a-1
	}
mn <- (max(snail[1,])-min(snail[1,])+max(snail[2,])-min(snail[2,]))/2
snail <- size*snail/mn
snail[1,] <- abc*snail[1,]
snail[2,] <- ord*snail[2,]
snail[1,] <- x+snail[1,]
snail[2,] <- y+snail[2,]
return (snail)	
}

brachiopod_symbol <- function(x,y,abc,ord,pedicle_width=5,sulcus_width=10,sulcus_depth=0.1,shell_width=2,size)	{
# have the angle fo the hinge start high and decrease
#seq(sin(pedicle_angle),0,by=-sin(pedicle_angle)/60)
#yy <- seq(sin(pedicle_angle),0,by=-sin(pedicle_angle)/60)/60
#for (i in 2:length(yy))	yy[i] <- yy[i]+yy[i-1]
#plot(xx,yy);
#pa <- pedicle_angle

## draw the "base" (hinge) of the shell
xx <- seq(0,shell_width/2,by=(shell_width/2)/59)
yy <- c()
pa <- 0.01
for (i in 1:pedicle_width)	{
	yy <- c(yy,sin(pa));
	pa <- 1.3*pa;
#	i <- i+1;
	}
for (i in (pedicle_width+1):60)	{
	yy <- c(yy,sin(pa));
	pa <- 0.95*pa;
#	i <- i+1;
	}
for (i in 2:length(yy))	yy[i] <- yy[i]+yy[i-1];
yy <- 0.2*yy
#plot(xx,yy,xlim=c(-1,1),ylim=c(-1,1),type="l");
radius <- ((xx[length(xx)]^2)+(yy[length(yy)])^2)^0.5;
base <- tan(yy[length(yy)]/xx[length(xx)]);

### draw the outer margin up to the sulcus
sul_wdth <- pi*sulcus_width/180
top <- (pi/2)-sul_wdth;
edges <- seq(base+((top-base)/108),top,(top-base)/108);
yyy <- radius*sin(edges);
xxx <- radius*cos(edges);
#lines(xxx,yyy);

if (sulcus_width > 0)	{
	### draw the sulcus
	edges2 <- seq(0+(pi/20),pi/2,by=pi/20);
	xxx <- c(xxx,0,xxx[length(xxx)])
	yyy <- c(yyy,0,yyy[length(yyy)])
	xxxx <- xxx[length(xxx)]*cos(edges2);
	yyyy <- max(yy,yyy)-xxx[length(xxx)]*sin(edges2);
	brachiopod_x <- c(xx,xxx,xxxx);
	brachiopod_y <- c(yy,yyy,yyyy);
	} else	{
	# if there is no sulcus, then we skip that part
	brachiopod_x <- c(xx,xxx);
	brachiopod_y <- c(yy,yyy);
	}

# draw the left side of the shell
brachiopod_x <- c(brachiopod_x,-brachiopod_x[length(brachiopod_x):1]);
brachiopod_y <- c(brachiopod_y,brachiopod_y[length(brachiopod_y):1]);

#recenter and rescale
brachiopod_y <- brachiopod_y-mean(brachiopod_y);
brachiopod <- rbind(brachiopod_x,brachiopod_y);
brachiopod <- size*brachiopod;
brachiopod[1,] <- abc*brachiopod[1,];
brachiopod[2,] <- ord*brachiopod[2,];
brachiopod[1,] <- x+brachiopod[1,];
brachiopod[2,] <- y+brachiopod[2,];
return(brachiopod);
}

trilobite_symbol <- function(x,y,abc,ord,size,cph,tho,thow,pyg,pygwu,pygwb)	{
# x, y: x & y coabcissa (center of symbol)
# abc, ord: length of x & y axes (to make it "square")
# cph: relative height of cephalon (head)
# tho: relative height of thorax
# thow: relative width of thorax
# pyg: relative height of pygidium (tail)
# pygwu: relative width of pygidium at top
# pygwb: relative width of pygidium at bottom
l <- cph+tho+pyg	# total length
cph <- cph/l
tho <- tho/l
pyg <- pyg/l
hby <- (0.5-cph)				#where body starts relative to y
bty <- (0.5-(cph+tho))	#where body ends relative to y
a <- 0.5
b <- cph
xxx <- vector(length=34)
yyy <- vector(length=34)
trilo <- matrix(0,2,62)
j <- 1
trilo[1,j+0] <- -(thow*a)	#head/body connction
trilo[1,j+1] <- -(((1+thow)/2)*a)	#head/body connction
trilo[2,j+3] <- trilo[2,j+1] <- trilo[2,j+0] <- hby
trilo[1,j+3] <- trilo[1,j+2] <- -a		#tip of spine
trilo[2,j+2] <- -(hby-(bty/2))
ang <- 175
j <- 4
for (i in 1:35)	{
	# go in 5?? increments
	#	ddd[i] <- (0.25*(abs(ang-90)/90))
	trilo[1,i+j] <- (a*cos(ang*pi/180))
	xxx[i] <- a*cos(ang*pi/180)
	yyy[i] <- sqrt(((b^2)-((b^2)*xxx[i]^2)/(a^2)))
	trilo[2,i+j] <- trilo[2,1]+yyy[i]
	ang <- ang-5	
	}
  j <- j+35
trilo[1,j+3] <- (thow*a)	#head/body connction
trilo[1,j+2] <- (((1+thow)/2)*a)	#head/body connction
trilo[2,j] <- trilo[2,j+2] <- trilo[2,j+3] <- hby
trilo[1,j] <- trilo[1,j+1] <- a		#tip of spine
trilo[2,j+1] <- -(hby-(bty/2))
  
j <- j+4
trilo[1,j] <- trilo[1,j+2] <- trilo[1,j+4] <- trilo[1,j+6] <- ((7/6)*thow*a)
trilo[1,j+1] <- trilo[1,j+3] <- trilo[1,j+5] <- trilo[1,j+7] <- (thow*a)
tip <- (hby-bty)/8
for (i in 1:8)	{
	trilo[2,j] <- hby-(i*tip)
	j <- j+1
	}
rev <- j-1
# do pygidium
trilo[1,j] <- pygwu*a
trilo[2,j] <- bty
j <- j+1
trilo[1,j] <- pygwb*a
trilo[2,j] <- -0.5
j <- j+1
trilo[1,j] <- -pygwb*a
trilo[2,j] <- -0.5
j <- j+1
trilo[1,j] <- -pygwu*a
trilo[2,j] <- bty
j <- j+1
  
# do rest of thorax
for (i in 1:8)	{
	k <- i-1
	trilo[1,j] <- -trilo[1,(rev-k)]
	trilo[2,j] <- trilo[2,(rev-k)]
	j <- j+1
	}
  
# center and rescale
trilo <- size*trilo
trilo[1,] <- abc*trilo[1,]
trilo[2,] <- ord*trilo[2,]
trilo[1,] <- x+trilo[1,]
trilo[2,] <- y+trilo[2,]
  
#polygon(trilo[1,],trilo[2,],col="orange",lwd=0.25)
return(trilo)
}

fish_symbol <- function(x,y,abc,ord,size,a,b)	{
  #	a <- 0.5
  #	b <- 0.25
  spin <- 25
  rev <- 3+(2*spin)
  xxx <- vector(length=52)
  yyy <- vector(length=52)
  #ddd <- vector(length=34)
  fish <- matrix(0,2,rev)
  ang <- 160
  for (i in 1:spin)	{
    # go in 5?? increments
    #	ddd[i] <- (0.25*(abs(ang-90)/90))
    fish[1,i] <- fish[1,rev-i] <- (a*cos(ang*pi/180))
    xxx[i] <- a*cos(ang*pi/180)
    yyy[i] <- sqrt(((b^2)-((b^2)*xxx[i]^2)/(a^2)))
    fish[2,i] <- yyy[i]
    fish[2,rev-i] <- yyy[i]*-1
    ang <- ang-5
  }
  fish[1,i+1] <- fish[1,i+2] <- 0.5
  fish[2,i+1] <- 1.5*max(fish[2,])
  fish[2,i+2] <- fish[2,i+1]*-1
  fish[1,rev] <- -0.30
  fish[2,rev] <- 0
  
  # center and rescale fish
  scx <- max(fish[1,])-min(fish[1,])
  fish <- fish/scx
  scxx <- max(fish[1,])-0.5
  fish[1,] <- fish[1,]-scxx
  
  # set the size so that the area is about equal to a circle with r=a
  fish <- (a/b)*size*fish
  fish[1,] <- abc*fish[1,]
  fish[2,] <- ord*fish[2,]
  fish[1,] <- x+fish[1,]
  fish[2,] <- y+fish[2,]
  return (fish)
}

log_axes <- function(axe,min_ax,max_ax,increment,numbers,linewd=4/3,orient)	{
axis(axe,at=seq(min_ax,max_ax,by=(max_ax-min_ax)),tcl=0,labels=NA,lwd.ticks=NA,lwd=linewd)
ticks <- length(numbers)
#if ((numbers[2]/numbers[1])>=increment)	{
for (i in 1:ticks)	{
	l <- numbers[i]
	axis(2,at=log(l),tcl=-.3,labels=numbers[i],lwd=0,lwd.ticks=linewd,las=orient)
	}
for (i in 1:(ticks-1))	{
	l <- numbers[i]+increment
	axis(2,at=log(l),tcl=-.15,labels=NA,lwd=0,lwd.ticks=linewd)
	}
}
#log10_axes(axe=2,min_ax=mny,max_ax=mxy,numbers=10^(mny:mxy),linewd=4/3,orient=2)
#	produces y-axis with 10^-7 to 1.0 with increments from 2:9 & labels on 10^x for x=integer

log10_axes <- function(axe,min_ax,max_ax,numbers,linewd=4/3,orient)	{
#log10_axes(axe=2,min_ax=mny,max_ax=mxy,numbers,linewd=1.5,orient=2)
# axe: 1 for x; 2 for y
# min_ax: minimum, already log10 transformed
# max_ax: maximum, already log10 transformed
# numbers: array of numbers
# linewd: line width
# orient: orientation of text, with 2 making y-axis the way I like it
axis(axe,at=seq(min_ax,max_ax,by=(max_ax-min_ax)),tcl=0,labels=NA,lwd.ticks=NA, lwd=linewd, las=orient)
#axis(axe,at=seq(min_ax,max_ax,by=(max_ax-min_ax)),tcl=0.0,labels=FALSE,lwd=linewd,lwd.ticks=0.0,las=orient)

#if ((numbers[2]/numbers[1])>=increment)	{
#mnn <- ceiling(log10(min(numbers)))
#mxn <- floor(log10(max(numbers)))
#for (i in 1:length(numbers))	{
#	l <- numbers[i]
#	axis(axe,at=log10(l),tcl=-.3,labels=numbers[i],lwd=0,lwd.ticks=linewd,las=orient)
#	}
mnn <- ceiling(min_ax)
mxn <- floor(max_ax)
for (i in mnn:mxn)	{
#	l <- numbers[i]
	axis(axe,at=i,tcl=-.3,labels=FALSE,lwd=0,lwd.ticks=linewd)
	}
for (i in 1:length(numbers))	{
	l <- numbers[i]
	axis(axe,at=log10(l),tcl=-.3,labels=numbers[i],lwd=0,lwd.ticks=0,las=orient)
	}
if (ceiling(min_ax)>min_ax)	{
	strt <- match(min_ax,ceiling(min_ax)-(1-log10(1:9)))
	if (is.na(strt))	{
		strt <- 1+sum((ceiling(min_ax)-(1-log10(1:9)))<min_ax)
		}
	ticks <- (strt:10)*10^floor(min_ax)
	} else	{
	ticks <- c()
	}
if (ceiling(min_ax)<floor(max_ax))	{
	for (i in ceiling(min_ax):(floor(max_ax)-1))	{
		ticks <- c(ticks,(1:10)*10^i)
		}
	}
if (max_ax>floor(max_ax))	{
	add <- max(ticks)
	while (log10(max(ticks)+add)<=max_ax)	ticks <- c(ticks,max(ticks)+add)
#	end <- match(round(max_ax,4),round(floor(max_ax)+log10(1:9),4))
#	ticks <- c(ticks,(1:end)*10^floor(max_ax))
	}
ticks <- ticks[!ticks %in% 10^(mnn:mxn)]
#axis(axe,at=log10(0.2),tcl=-.15,labels=FALSE,lwd=0,lwd.ticks=linewd)
#axis(axe,at=log10(0.4),tcl=-.15,labels=FALSE,lwd=0,lwd.ticks=linewd)
axis(axe,at=log10(ticks),tcl=-.15,labels=FALSE,lwd=0,lwd.ticks=linewd)
}

log2_axes <- function(axe,min_ax,max_ax,numbers,linewd=4/3,orient)	{
#log2_axes(axe=2,min_ax=mny,max_ax=mxy,numbers,linewd=1.5,orient=2)
# axe: 1 for x; 2 for y
# min_ax: minimum, already log2 transformed
# max_ax: maximum, already log2 transformed
# numbers: array of numbers
# linewd: line width
# orient: orientation of text, with 2 making y-axis the way I like it
axis(axe,at=seq(min_ax,max_ax,by=(max_ax-min_ax)),tcl=0,labels=NA,lwd.ticks=NA, lwd=linewd, las=orient)
ticks <- length(numbers)
for (i in 1:ticks)	{
	if (numbers[1]==0)	{
		l <- max(1,2*numbers[i])
		}	else	l <- numbers[i]
	axis(axe,at=log2(l),tcl=-.3,labels=numbers[i],lwd=0,lwd.ticks=linewd,las=orient)
#	i <- i+1
	}
}

spindle_diagram <- function(bin_onsets,spindled_midpts,spindled_counts,bin_colors,plot_on_y=T,bar_legend=T,legend_width=1,legend_case="Unit")	{
# bin_onsets: vector giving where (on X or Y axis)
# spindled_midpts: the position where the middle of each spindle segment is plotted
# spindled_counts: matrix giving the number of times an observation is made
#	with the mean of that variable
# bin_colors: colors for separate spindle diagrams
# plot_on_y: if TRUE, then spindles go "up" Y-axis, with different ones plotted
#	on X-axis
# bar_legend: if TRUE,then print a legend for width
# legend_width: the width being plotted
# legend_case: Unit name (e.g., "case" or "clade") that one example represents
bins <- dim(spindled_counts)[1];
bin_mids <- (bin_onsets[1:bins]+bin_onsets[2:(bins+1)])/2;
bin_widths <- abs(bin_onsets[1:bins]-bin_onsets[2:(bins+1)]);
bin_maxs <- vector(length=bins)
for (i in 1:bins)	bin_maxs[i] <- max(spindled_counts[i,]);
width_of_one <- min(bin_widths/bin_maxs);
y_hts <- abs(spindled_midpts[1]-spindled_midpts[2])/2
#cbind(bin_maxs,bin_widths,bin_widths/bin_maxs)
#which(spindled_counts==mxwd,arr.ind=TRUE)
print(c(width_of_one,legend_width));
for (s in 1:bins)	{
	toplot <- spindled_midpts[spindled_counts[s,]>0]
	tocount <- spindled_counts[s,spindled_counts[s,]>0]
	for (i in 1:length(toplot))	{
		if (plot_on_y)	{
			x1 <- as.numeric(bin_mids[s]-(width_of_one/2)*tocount[i])
			x2 <- as.numeric(bin_mids[s]+(width_of_one/2)*tocount[i])
			y1 <- as.numeric(toplot[i])-y_hts
			y2 <- as.numeric(toplot[i])+y_hts
			}	else	{
			y1 <- as.numeric(bin_mids[s]-(width_of_one/2)*tocount[i])
			y2 <- as.numeric(bin_mids[s]+(width_of_one/2)*tocount[i])
			x1 <- as.numeric(toplot[i])-y_hts
			x2 <- as.numeric(toplot[i])+y_hts
			}
		rect(x1,y1,x2,y2,col=bin_colors[s],lwd=0.5)
		}
	}
# do legend if requested
if (bar_legend)	{
	xleg <- max(bin_onsets)-(0.2*(max(bin_onsets)-min(bin_onsets)));
	if (plot_on_y)	{
		x1 <- xleg-((width_of_one/2)*legend_width);
		x2 <- xleg+((width_of_one/2)*legend_width);
		x3 <- xleg+((width_of_one/2)*legend_width);
		y1 <- max(spindled_midpts)+y_hts;
		y2 <- max(spindled_midpts)-y_hts;
		y3 <- (y1+y2)/2;
		}	else	{
		y1 <- xleg-((width_of_one/2)*legend_width)
		y2 <- xleg+((width_of_one/2)*legend_width)
		y3 <- xleg+((width_of_one/2)*legend_width)
		x1 <- max(spindled_midpts)+y_hts
		x2 <- max(spindled_midpts)-y_hts
		x3 <- (x1+x2)/2
		}
	rect(x1,y1,x2,y2,col="gray50",lwd=0.5)
	legend_text <- paste(":",legend_width,legend_case,sep=" ")
	text(x3,y3,legend_text,cex=0.75,pos=4)
	}
}

single_spindle <- function(axe=2,midpoint,spindle,max_width,min_axe,max_axe,spindle_color="white",spindle_lwdth=4/3,spindle_linecol="black")	{
# axe: 1 for x-axis (going left-right), 2 for y-axis (going up-down)
# midpoint: midpoint on the other axis
# spindle: width of spindel at each point along axis (= height of
#	histogram for the same data)
# max_width: maximum width on the other axis
# spindle_color: color of spindel bars
# spindle_lwdth: width of lines around spindel
# spindel_linecol: color of lines
spinds <- length(spindle)
axe_incr <- seq(min_axe,max_axe,by=(max_axe-min_axe)/spinds)
scaled_spindle <- max_width*spindle
for (sp in 1:spinds)	{
	if (scaled_spindle[sp]>0 && (axe==2 || axe==4))	{
		rect(midpoint-scaled_spindle[sp]/2,axe_incr[sp],midpoint+scaled_spindle[sp]/2,axe_incr[sp+1],col=spindle_color,lwd=4/3,border=spindle_color)
		} else if (scaled_spindle[sp]>0) {
		rect(axe_incr[sp],midpoint-scaled_spindle[sp]/2,axe_incr[sp+1],midpoint+scaled_spindle[sp]/2,col=spindle_color,lwd=4/3,border=spindle_color)
		}
	}
if (axe==2 || axe==4)	{
	for (sp in 1:spinds)	{
		# draw bottom if there is no segment beneath it
		if (scaled_spindle[sp]>0)	{
			if (sp==1 || (sp>1 && scaled_spindle[sp-1]==0))
				segments(midpoint-scaled_spindle[sp]/2,axe_incr[sp],midpoint+scaled_spindle[sp]/2,axe_incr[sp],lwd=spindle_lwdth,col=spindle_linecol)
			# draw sides
			segments(midpoint-scaled_spindle[sp]/2,axe_incr[sp],midpoint-scaled_spindle[sp]/2,axe_incr[sp+1],lwd=spindle_lwdth,col=spindle_linecol)
			segments(midpoint+scaled_spindle[sp]/2,axe_incr[sp],midpoint+scaled_spindle[sp]/2,axe_incr[sp+1],lwd=spindle_lwdth,col=spindle_linecol)
			# draw top if necessary
			if (sp==spinds || (sp<spinds && scaled_spindle[sp+1]==0))	{
				segments(midpoint-scaled_spindle[sp]/2,axe_incr[sp+1],midpoint+scaled_spindle[sp]/2,axe_incr[sp+1],lwd=spindle_lwdth,col=spindle_linecol)
				} else if (scaled_spindle[sp]>scaled_spindle[sp+1])	{
#				dev <- (scaled_spindle[sp]-scaled_spindle[sp+1])/2
				segments(midpoint-scaled_spindle[sp]/2,axe_incr[sp+1],midpoint-scaled_spindle[sp+1]/2,axe_incr[sp+1],lwd=spindle_lwdth,col=spindle_linecol)
				segments(midpoint+scaled_spindle[sp]/2,axe_incr[sp+1],midpoint+scaled_spindle[sp+1]/2,axe_incr[sp+1],lwd=spindle_lwdth,col=spindle_linecol)
				} else if (scaled_spindle[sp]<scaled_spindle[sp+1])	{
				segments(midpoint-scaled_spindle[sp]/2,axe_incr[sp+1],midpoint-scaled_spindle[sp+1]/2,axe_incr[sp+1],lwd=spindle_lwdth,col=spindle_linecol)
				segments(midpoint+scaled_spindle[sp]/2,axe_incr[sp+1],midpoint+scaled_spindle[sp+1]/2,axe_incr[sp+1],lwd=spindle_lwdth,col=spindle_linecol)
				}
			}
		}
#	top <- 1+max(seq(1:spinds)*(scaled_spindle>0))
#	segments(midpoint-scaled_spindle[top]/2,min_axe,midpoint+scaled_spindle[1]/2,min_axe,lwd=spindle_lwdth,col=spindle_linecol)
	}	else	{
	for (sp in 1:spinds)	{
		# draw bottom if there is no segment beneath it
		if (scaled_spindle[sp]>0)	{
			if (sp==1 || (sp>1 && scaled_spindle[sp-1]==0))
				segments(axe_incr[sp],midpoint-scaled_spindle[sp]/2,axe_incr[sp],midpoint+scaled_spindle[sp]/2,lwd=spindle_lwdth,col=spindle_linecol)
			# draw sides
			segments(axe_incr[sp],midpoint-scaled_spindle[sp]/2,axe_incr[sp+1],midpoint-scaled_spindle[sp]/2,lwd=spindle_lwdth,col=spindle_linecol)
			segments(midpoint+scaled_spindle[sp]/2,axe_incr[sp],midpoint+scaled_spindle[sp]/2,axe_incr[sp+1],lwd=spindle_lwdth,col=spindle_linecol)
			# draw top if necessary
			if (sp==spinds || (sp<spinds && scaled_spindle[sp+1]==0))	{
				segments(axe_incr[sp+1],midpoint-scaled_spindle[sp]/2,axe_incr[sp+1],midpoint+scaled_spindle[sp]/2,lwd=spindle_lwdth,col=spindle_linecol)
				} else if (scaled_spindle[sp]>scaled_spindle[sp+1])	{
#				dev <- (scaled_spindle[sp]-scaled_spindle[sp+1])/2
				segments(axe_incr[sp+1],midpoint-scaled_spindle[sp]/2,axe_incr[sp+1],midpoint-scaled_spindle[sp+1]/2,lwd=spindle_lwdth,col=spindle_linecol)
				segments(axe_incr[sp+1],midpoint+scaled_spindle[sp]/2,axe_incr[sp+1],midpoint+scaled_spindle[sp+1]/2,lwd=spindle_lwdth,col=spindle_linecol)
				} else if (scaled_spindle[sp]<scaled_spindle[sp+1])	{
				segments(axe_incr[sp+1],midpoint-scaled_spindle[sp]/2,axe_incr[sp+1],midpoint-scaled_spindle[sp+1]/2,lwd=spindle_lwdth,col=spindle_linecol)
				segments(axe_incr[sp+1],midpoint+scaled_spindle[sp]/2,axe_incr[sp+1],midpoint+scaled_spindle[sp+1]/2,lwd=spindle_lwdth,col=spindle_linecol)
				}
			}
		}
#	top <- 1+max(seq(1:spinds)*(scaled_spindle>0))
#	segments(midpoint-scaled_spindle[top]/2,min_axe,midpoint+scaled_spindle[1]/2,min_axe,lwd=spindle_lwdth,col=spindle_linecol)
	}
}

single_spindle_set_increments <- function(axe=2,midpoint,axe_incr,spindle,max_width,spindle_color="white",spindle_lwdth=4/3,spindle_linecol="black")	{
# axe: 1 for x-axis (going left-right), 2 for y-axis (going up-down)
# midpoint: midpoint on the other axis
# spindle: width of spindel at each point along axis (= height of
#	histogram for the same data)
# max_width: maximum width on the other axis
# spindle_color: color of spindel bars
# spindle_lwdth: width of lines around spindel
# spindel_linecol: color of lines
spindle <- spindle/max(spindle)
spinds <- length(spindle)
scaled_spindle <- max_width*spindle
for (sp in 1:spinds)	{
	if (scaled_spindle[sp]>0 && (axe==2 || axe==4))	{
		rect(midpoint-scaled_spindle[sp]/2,axe_incr[sp],midpoint+scaled_spindle[sp]/2,axe_incr[sp+1],col=spindle_color,lwd=4/3,border=spindle_color)
		} else if (scaled_spindle[sp]>0) {
		rect(axe_incr[sp],midpoint-scaled_spindle[sp]/2,axe_incr[sp+1],midpoint+scaled_spindle[sp]/2,col=spindle_color,lwd=4/3,border=spindle_color)
		}
	}
if (axe==2 || axe==4)	{
	for (sp in 1:spinds)	{
		# draw bottom if there is no segment beneath it
		if (scaled_spindle[sp]>0)	{
			if (sp==1 || (sp>1 && scaled_spindle[sp-1]==0))
				segments(midpoint-scaled_spindle[sp]/2,axe_incr[sp],midpoint+scaled_spindle[sp]/2,axe_incr[sp],lwd=spindle_lwdth,col=spindle_linecol)
			# draw sides
			segments(midpoint-scaled_spindle[sp]/2,axe_incr[sp],midpoint-scaled_spindle[sp]/2,axe_incr[sp+1],lwd=spindle_lwdth,col=spindle_linecol)
			segments(midpoint+scaled_spindle[sp]/2,axe_incr[sp],midpoint+scaled_spindle[sp]/2,axe_incr[sp+1],lwd=spindle_lwdth,col=spindle_linecol)
			# draw top if necessary
			if (sp==spinds || (sp<spinds && scaled_spindle[sp+1]==0))	{
				segments(midpoint-scaled_spindle[sp]/2,axe_incr[sp+1],midpoint+scaled_spindle[sp]/2,axe_incr[sp+1],lwd=spindle_lwdth,col=spindle_linecol)
				} else if (scaled_spindle[sp]>scaled_spindle[sp+1])	{
#				dev <- (scaled_spindle[sp]-scaled_spindle[sp+1])/2
				segments(midpoint-scaled_spindle[sp]/2,axe_incr[sp+1],midpoint-scaled_spindle[sp+1]/2,axe_incr[sp+1],lwd=spindle_lwdth,col=spindle_linecol)
				segments(midpoint+scaled_spindle[sp]/2,axe_incr[sp+1],midpoint+scaled_spindle[sp+1]/2,axe_incr[sp+1],lwd=spindle_lwdth,col=spindle_linecol)
				} else if (scaled_spindle[sp]<scaled_spindle[sp+1])	{
				segments(midpoint-scaled_spindle[sp]/2,axe_incr[sp+1],midpoint-scaled_spindle[sp+1]/2,axe_incr[sp+1],lwd=spindle_lwdth,col=spindle_linecol)
				segments(midpoint+scaled_spindle[sp]/2,axe_incr[sp+1],midpoint+scaled_spindle[sp+1]/2,axe_incr[sp+1],lwd=spindle_lwdth,col=spindle_linecol)
				}
			}
		}
#	top <- 1+max(seq(1:spinds)*(scaled_spindle>0))
#	segments(midpoint-scaled_spindle[top]/2,min_axe,midpoint+scaled_spindle[1]/2,min_axe,lwd=spindle_lwdth,col=spindle_linecol)
	}	else	{
	for (sp in 1:spinds)	{
		# draw bottom if there is no segment beneath it
		if (scaled_spindle[sp]>0)	{
			if (sp==1 || (sp>1 && scaled_spindle[sp-1]==0))
				segments(axe_incr[sp],midpoint-scaled_spindle[sp]/2,axe_incr[sp],midpoint+scaled_spindle[sp]/2,lwd=spindle_lwdth,col=spindle_linecol)
			# draw sides
			segments(axe_incr[sp],midpoint-scaled_spindle[sp]/2,axe_incr[sp+1],midpoint-scaled_spindle[sp]/2,lwd=spindle_lwdth,col=spindle_linecol)
			segments(midpoint+scaled_spindle[sp]/2,axe_incr[sp],midpoint+scaled_spindle[sp]/2,axe_incr[sp+1],lwd=spindle_lwdth,col=spindle_linecol)
			# draw top if necessary
			if (sp==spinds || (sp<spinds && scaled_spindle[sp+1]==0))	{
				segments(axe_incr[sp+1],midpoint-scaled_spindle[sp]/2,axe_incr[sp+1],midpoint+scaled_spindle[sp]/2,lwd=spindle_lwdth,col=spindle_linecol)
				} else if (scaled_spindle[sp]>scaled_spindle[sp+1])	{
#				dev <- (scaled_spindle[sp]-scaled_spindle[sp+1])/2
				segments(axe_incr[sp+1],midpoint-scaled_spindle[sp]/2,axe_incr[sp+1],midpoint-scaled_spindle[sp+1]/2,lwd=spindle_lwdth,col=spindle_linecol)
				segments(axe_incr[sp+1],midpoint+scaled_spindle[sp]/2,axe_incr[sp+1],midpoint+scaled_spindle[sp+1]/2,lwd=spindle_lwdth,col=spindle_linecol)
				} else if (scaled_spindle[sp]<scaled_spindle[sp+1])	{
				segments(axe_incr[sp+1],midpoint-scaled_spindle[sp]/2,axe_incr[sp+1],midpoint-scaled_spindle[sp+1]/2,lwd=spindle_lwdth,col=spindle_linecol)
				segments(axe_incr[sp+1],midpoint+scaled_spindle[sp]/2,axe_incr[sp+1],midpoint+scaled_spindle[sp+1]/2,lwd=spindle_lwdth,col=spindle_linecol)
				}
			}
		}
#	top <- 1+max(seq(1:spinds)*(scaled_spindle>0))
#	segments(midpoint-scaled_spindle[top]/2,min_axe,midpoint+scaled_spindle[1]/2,min_axe,lwd=spindle_lwdth,col=spindle_linecol)
	}
}

wagner_set_axes <- function (ax_min,ax_max)	{
if ((ax_max-ax_min)<=10)	{
	min_ax <- floor(ax_min)
	max_ax <- ceiling(ax_max)
	lbl_prn <- tcs <- c()
	tcs <- rbind(tcs,seq(min_ax,max_ax,by=1))
	lbl_prn <- rbind(lbl_prn,-1*(min_ax:max_ax))
	} else if ((ax_max-ax_min)<=25)	{
	min_ax <- floor(ax_min)
	max_ax <- ceiling(ax_max)
	tcs <- seq(2*ceiling(min_ax/2),2*floor(max_ax/2),by=2)
	lbl_prn <- -1*tcs
	add <- length(tcs)-length((min_ax:max_ax)[!(min_ax:max_ax) %in% tcs])
	added <- c((min_ax:max_ax)[!(min_ax:max_ax) %in% tcs])
	if (add==0)	{
		tcs <- rbind(tcs,added)
		lbl_prn <- rbin(lbl_prn,rep("",length(added)))
		}	else if (add>0)	{
		added <- c(added,rep(0,add))
		tcs <- rbind(tcs,added)
		lbl_prn <- rbind(lbl_prn,rep("",length(added)))
		} else	{
		tcs <- rbind(c(tcs,rep(0,abs(add))),added)
		lbl_prn <- rbind(c(lbl_prn,rep("",abs(add)),rep("",length(added))))
		}
	tick_str <- c(-0.30,-0.15)
	} else if ((ax_max-ax_min)<=50)	{
	min_ax <- floor(ax_min)
	max_ax <- ceiling(max(divergences))
	tcs <- seq(5*ceiling(min_ax/5),5*floor(max_ax/5),by=5)
	lbl_prn <- -1*tcs
	add <- length(tcs)-length((min_ax:max_ax)[!(min_ax:max_ax) %in% tcs])
	added <- c((min_ax:max_ax)[!(min_ax:max_ax) %in% tcs])
	if (add==0)	{
		tcs <- rbind(tcs,added)
		lbl_prn <- rbin(lbl_prn,rep("",length(added)))
		}	else if (add>0)	{
		added <- c((min_ax:max_ax)[!(min_ax:max_ax) %in% tcs],rep(0,abs(add)))
		tcs <- rbind(tcs,added)
		lbl_prn <- rbind(c(lbl_prn,rep("",add)),rep("",length(added)))
		}	else	{
		tcs <- rbind(c(tcs,rep(0,abs(add))),added)
		lbl_prn <- rbind(c(lbl_prn,rep("",abs(add))),rep("",length(added)))
		}
	tick_str <- c(-0.30,-0.15)
#	tcs <- rbind(tcs,(min_ax:max_ax)[!(min_ax:max_ax) %in% tcs])
#	lbl_prn <- c(TRUE,FALSE)
	}	else if ((ax_max-ax_min)<=100)	{
	min_ax <- floor(ax_min)
	max_ax <- ceiling(max(divergences))
	tcs <- seq(10*ceiling(min_ax/10),10*floor(max_ax/10),by=10)
	tcs <- rbind(tcs,seq(5*ceiling(min_ax/5),5*floor(max_ax/5),by=5)[!seq(5*ceiling(min_ax/5),5*floor(max_ax/5),by=5) %in% seq(10*ceiling(min_ax/10),10*floor(max_ax/10),by=10)])
	tcs <- rbind(tcs,(min_ax:max_ax)[!(min_ax:max_ax) %in% seq(5*ceiling(min_ax/5),5*floor(max_ax/5),by=5)])
	tick_str <- c(-0.30,-0.20,-0.10)
	lbl_prn <- c(TRUE,FALSE,FALSE)
	}
output <- list(tcs,lbl_prn,tick_str)
names(output) <- c("Ticks","Labels","Tick_Strength")
return(output)
}

fitted_linear_axis <- function(axe,max_val,min_val,linewd=4/3,orient=1,decimals=TRUE)	{
if (max_val<=10)	{
	axis(axe,at=seq(min_val,max_val,by=max_val-min_val),tcl=-0.0,labels = FALSE,lwd=linewd,las=orient)
	if (decimals==TRUE)	{
		ticks <- old_ticks <- seq(min_val,max_val,by=1)
		axis(axe,at=ticks,tcl=-0.30,labels=FALSE,lwd=0.0,lwd.ticks=linewd,las=orient)
		bs <- 0.5
		ticks <- seq(min_val,max_val,by=bs)[!seq(min_val,max_val,by=bs) %in% old_ticks]
		axis(axe,at=ticks,tcl=-0.20,labels=FALSE,lwd=0.0,lwd.ticks=linewd,las=orient)
		axis(axe,at=seq(min_val,max_val,by=bs),tcl=-0.0,labels=TRUE,lwd=0.0,lwd.ticks=0.0,las=orient)
		bs <- 0.1
		ticks <- seq(min_val,max_val,by=bs)[!seq(min_val,max_val,by=bs) %in% old_ticks]
		axis(axe,at=ticks,tcl=-0.10,labels=FALSE,lwd=0.0,lwd.ticks=linewd,las=orient)
		} else	{
		axis(axe,at=seq(min_val,max_val,by=1),tcl=-0.3,labels=TRUE,lwd=0.0,lwd.ticks=linewd,las=orient)
		}
	} else if (max_val <=20)	{
	axis(axe,at=seq(min_val,max_val,by=(max_val-min_val)),tcl=0.0,labels=FALSE,lwd=linewd,las=orient)
	ticks <- seq(min_val,max_val,by=2)
	axis(axe,at=ticks,tcl=-0.30,labels=TRUE,lwd=0.0,lwd.ticks=linewd,las=orient)
	old_ticks <- ticks
	bs <- 1
	ticks <- seq(min_val,max_val,by=bs)[!seq(min_val,max_val,by=bs) %in% old_ticks]
	axis(axe,at=ticks,tcl=-0.15,labels=FALSE,lwd=0.0,lwd.ticks=linewd)
	} else if (max_val <=50)	{
	axis(axe,at=seq(min_val,max_val,by=(max_val-min_val)),tcl=0.0,labels=FALSE,lwd=linewd,las=orient)
	ticks <- seq(min_val,max_val,by=5)
	axis(axe,at=ticks,tcl=-0.30,labels=TRUE,lwd=0.0,lwd.ticks=linewd,las=orient)
	old_ticks <- ticks
	bs <- 1
	ticks <- seq(min_val,max_val,by=bs)[!seq(min_val,max_val,by=bs) %in% old_ticks]
	axis(axe,at=ticks,tcl=-0.15,labels=FALSE,lwd=0.0,lwd.ticks=linewd)
	} else if (max_val <=100)	{
	axis(axe,at=seq(min_val,max_val,by=(max_val-min_val)),tcl=0.0,labels=FALSE,lwd=linewd,las=orient)
	ticks <- seq(min_val,max_val,by=10)
	axis(axe,at=ticks,tcl=-0.30,labels=TRUE,lwd=0.0,lwd.ticks=linewd,las=orient)
	old_ticks <- ticks
	bs <- 5
	ticks <- seq(min_val,max_val,by=bs)[!seq(min_val,max_val,by=bs) %in% old_ticks]
	axis(axe,at=ticks,tcl=-0.20,labels=FALSE,lwd=0.0,lwd.ticks=linewd)
	old_ticks <- ticks
	bs <- 1
	ticks <- seq(min_val,max_val,by=bs)[!seq(min_val,max_val,by=bs) %in% old_ticks]
	axis(axe,at=ticks,tcl=-0.10,labels=FALSE,lwd=0.0,lwd.ticks=linewd)
	} else if (max_val <=250)	{
	axis(axe,at=seq(min_val,max_val,by=(max_val-min_val)),tcl=0.0,labels=FALSE,lwd=linewd,las=orient)
	ticks <- seq(min_val,max_val,by=25)
	axis(axe,at=ticks,tcl=-0.30,labels=TRUE,lwd=0.0,lwd.ticks=linewd,las=orient)
	old_ticks <- ticks
	bs <- 5
	ticks <- seq(min_val,max_val,by=bs)[!seq(min_val,max_val,by=bs) %in% old_ticks]
	axis(axe,at=ticks,tcl=-0.15,labels=FALSE,lwd=0.0,lwd.ticks=linewd)
	} else if (max_val <=500)	{
	axis(axe,at=seq(min_val,max_val,by=(max_val-min_val)),tcl=0.0,labels=FALSE,lwd=linewd,las=orient)
	ticks <- seq(min_val,max_val,by=50)
	axis(axe,at=ticks,tcl=-0.30,labels=TRUE,lwd=0.0,lwd.ticks=linewd,las=orient)
	old_ticks <- ticks
	bs <- 25
	ticks <- seq(min_val,max_val,by=bs)[!seq(min_val,max_val,by=bs) %in% old_ticks]
	axis(axe,at=ticks,tcl=-0.20,labels=FALSE,lwd=0.0,lwd.ticks=linewd)
	bs <- 5
	ticks <- seq(min_val,max_val,by=bs)[!seq(min_val,max_val,by=bs) %in% old_ticks]
	axis(axe,at=ticks,tcl=-0.10,labels=FALSE,lwd=0.0,lwd.ticks=linewd)
	} else if (max_val <=1000)	{
	axis(axe,at=seq(min_val,max_val,by=(max_val-min_val)),tcl=0.0,labels=FALSE,lwd=linewd,las=orient)
	ticks <- seq(min_val,max_val,by=100)
	axis(axe,at=ticks,tcl=-0.30,labels=TRUE,lwd=0.0,lwd.ticks=linewd,las=orient)
	old_ticks <- ticks
	bs <- 50
	ticks <- seq(min_val,max_val,by=bs)[!seq(min_val,max_val,by=bs) %in% old_ticks]
	axis(axe,at=ticks,tcl=-0.20,labels=FALSE,lwd=0.0,lwd.ticks=linewd)
	bs <- 10
	ticks <- seq(min_val,max_val,by=bs)[!seq(min_val,max_val,by=bs) %in% old_ticks]
	axis(axe,at=ticks,tcl=-0.10,labels=FALSE,lwd=0.0,lwd.ticks=linewd)
	}
}

specify_basic_plot <- function(mxx, mnx, mxy, mny, abcissa="", ordinate="", xsize=3, ysize=3)	{
par(pin=c(xsize,ysize));
plot(NA,type='n',axes=FALSE,main="",xlab=abcissa,ylab=ordinate,xlim=c(mnx,mxx),ylim=c(mny,mxy));
}

# routine to make axes as you want them.
# axe: axis # (1 = x; 2= y)
# max_val: maximum value
# min_val: minimum value
# maj_break: major (labelled) breaks
# med_break: intermediate breaks
# min_break: minor breaks.  NOTE: if med_break or min_break=0, then just two breaks
# specified_axis(axe=2,max_val=mxy,min_val=mny,maj_break=100,med_break=50,min_break=10,linewd=4/3,orient=2) gives y-axis with 0:700 labeled,
#	2ndary ticks on 50 & tertiary ticks on 10
specified_axis <- function(axe,max_val,min_val,maj_break,med_break,min_break,linewd=4/3,orient=1,print_label=TRUE)	{
axis(axe,at=seq(min_val,max_val,by=max_val-min_val),tcl=-0.0,labels = FALSE,lwd=linewd,las=orient)
if ((min_val/maj_break)-floor(min_val/maj_break)<(10^-10))	{
	# this is a kluge necessitated by tiny rounding errors....
	mnv1 <- min_val
	} else	{
	mnv1 <- (maj_break*ceiling(min_val/maj_break))
	}
ticks <- old_ticks <- seq(mnv1,max_val,by=maj_break)
axis(axe,at=ticks,tcl=-0.30,labels=print_label,lwd=0.0,lwd.ticks=linewd,las=orient)
if (med_break!=0)	{
	mnv2 <- (med_break*ceiling(min_val/med_break))
	ticks <- seq(mnv2,max_val,by=med_break)[!seq(mnv2,max_val,by=med_break) %in% old_ticks]
	if (min_break!=0)	{
		tck_sz <- -0.20
		} else	tck_sz <- -0.15
	axis(axe,at=ticks,tcl=tck_sz,labels=FALSE,lwd=0.0,lwd.ticks=linewd,las=orient)
	old_ticks <- sort(c(ticks,old_ticks))
	}
if (min_break!=0)	{
	ticks <- seq(min_val,max_val,by=min_break)[!seq(min_val,max_val,by=min_break) %in% old_ticks]
	if (med_break!=0)	{
		tck_sz <- -0.10
		} else	tck_sz <- -0.15
	axis(axe,at=ticks,tcl=tck_sz,labels=FALSE,lwd=0.0,lwd.ticks=linewd,las=orient)
	}
}

specified_axis_w_labels <- function(axe,max_val,min_val,maj_break,med_break=0,min_break=0,axis_labels,axis_label_pts=NULL,axis_label_size=1,linewd=4/3,orient=1)	{
axis(axe,at=seq(min_val,max_val,by=max_val-min_val),tcl=-0.0,labels = FALSE,lwd=linewd,las=orient);
#if ((min_val/maj_break)-floor(min_val/maj_break)<(10^-10))	{
	# this is a kluge necessitated by tiny rounding errors....
#	mnv1 <- min_val;
#	} else	{
#	mnv1 <- (maj_break*ceiling(min_val/maj_break));
#	}
ticks <- old_ticks <- seq(min_val,max_val,by=maj_break);
axis(axe,at=ticks,tcl=-0.30,labels=FALSE,lwd=0.0,lwd.ticks=linewd,las=orient);
if (med_break!=0)	{
	mnv2 <- (med_break*ceiling(min_val/med_break))
	ticks <- seq(mnv2,max_val,by=med_break)[!seq(mnv2,max_val,by=med_break) %in% old_ticks]
	if (min_break!=0)	{
		tck_sz <- -0.20
		} else	tck_sz <- -0.15
	axis(axe,at=ticks,tcl=tck_sz,labels=FALSE,lwd=0.0,lwd.ticks=linewd,las=orient)
	old_ticks <- sort(c(ticks,old_ticks))
	}
if (min_break!=0)	{
	ticks <- seq(min_val,max_val,by=min_break)[!seq(min_val,max_val,by=min_break) %in% old_ticks]
	if (med_break!=0)	{
		tck_sz <- -0.10
		} else	tck_sz <- -0.15
	axis(axe,at=ticks,tcl=tck_sz,labels=FALSE,lwd=0.0,lwd.ticks=linewd,las=orient)
	}
#axis(axe,at=axis_labels,tcl=-0.30,labels=TRUE,lwd=0.0,lwd.ticks=0.0,las=orient)
if(is.numeric(axis_labels[1]) && is.null(axis_label_pts))	{
	axis_label_pts <- axis_labels;
	}
axis(axe,at=axis_label_pts,tcl=-0.30,labels=axis_labels,lwd=0.0,lwd.ticks=0.0,las=orient,cex.axis=axis_label_size);
}

specified_right_y_axis <- function(mxy1,mny1,max_val,min_val,maj_break,med_break,min_break,linewd=4/3,orient=2,print_label=TRUE)	{
# mxy1: maximum value on main (left) y-axis
# mny1: minimum value on main (left) y-axis; we rescale (max_val-min_val)/(mxy1-mny1)
axis(side=4,at=seq(mny1,mxy1,by=(mxy-mny)),tcl=-0.0,labels = FALSE,lwd=linewd,las=orient)
if ((min_val/maj_break)-floor(min_val/maj_break)<(10^-10))	{
	# this is a kluge necessitated by tiny rounding errors....
	mnv1 <- min_val
	} else	{
	mnv1 <- (maj_break*ceiling(min_val/maj_break))
	}
labels <- seq(mnv1,max_val,by=maj_break)
ticks <- old_ticks <- mny1 + (mxy1-mny1)*labels/(max(labels)-min(labels))
#ticks <- old_ticks <- seq(mnv1,mxy1,by=(mxy1-mny1)/(max_val/maj_break))
if (print_label)	{
	axis(side=4,at=ticks,tcl=-0.30,labels=seq(mnv1,max_val,by=maj_break),lwd=0.0,lwd.ticks=linewd,las=orient)
	} else	{
	axis(side=4,at=ticks,tcl=-0.30,labels=print_label,lwd=0.0,lwd.ticks=linewd,las=orient)
	}

if (med_break!=0)	{
	ticks <- mny1 + (mxy1-mny1)*seq(min_val,max_val,by=med_break)/(max(labels)-min(labels))
	mnv2 <- (med_break*ceiling(min_val/med_break))
	ticks <- ticks[!ticks %in% old_ticks]
	if (min_break!=0 && min_break!=med_break)	{
		tck_sz <- -0.20
		} else	tck_sz <- -0.15
	axis(side=4,at=ticks,tcl=tck_sz,labels=FALSE,lwd=0.0,lwd.ticks=linewd,las=orient)
	old_ticks <- sort(c(ticks,old_ticks))
	}
if (min_break!=0 && min_break!=med_break)	{
	ticks <- mny1 + (mxy1-mny1)*seq(min_val,max_val,by=min_break)/(max(labels)-min(labels))
#	ticks <- seq(mny1,mxy1,by=(mxy1-mny1)/(max_val/min_break))[!seq(mny1,mxy1,by=(mxy1-mny1)/(max_val/min_break)) %in% old_ticks]
	ticks <- ticks[!ticks %in% old_ticks]
	if (med_break!=0)	{
		tck_sz <- -0.10
		} else	tck_sz <- -0.15
	axis(side=4,at=ticks,tcl=tck_sz,labels=FALSE,lwd=0.0,lwd.ticks=linewd,las=orient)
	}
}

#slice_coll
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
####  Phylogeny drawing routines ####
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
center_budding_phylogeny <- function(ctree,durations)	{
venntree <- accio_venn_tree_from_vector_tree(ctree);
mtree <- accio_matrix_tree_from_vector_tree(ctree);
node_richness <- accio_node_richness_from_vector_tree(vector_tree = ctree);
node_ages <- c();
for (nd in 1:nNodes)	node_ages <- c(node_ages,min(durations[venntree[nd,venntree[nd,]>0],1]));
branch_ages <- c(durations[,1],node_ages)
nNode <- nrow(mtree);
notu <- length(ctree)-nNode;
ttl_richness <- c(rep(1,notu),node_richness);
patristic_distances <- accio_patristic_distance_from_base(atree=mtree);
max_nodes <- max(patristic_distances);

last_left <- "left";	# move up the axis
last_right <- "right";	# move down the axis
accounted <- c();
nd <- 0;
phy_pos <- rep(0,nNode+notu);
for (nd in 1:nNodes)	{
#	nd <- nd+1;
	htu <- nd+notu;			# htu number of node;
	tf1 <- sum(mtree[nd,]>0);
	if (sampled_ancestors[htu]!=0)	{
		tf1 <- tf1-1;
		phy_pos[sampled_ancestors[htu]] <- phy_pos[notu+nd];
		}  
	f1 <- mtree[nd,!mtree[nd,] %in% c(sampled_ancestors[htu],0)];
	f1 <- f1[order(-ttl_richness[f1])];
	if (length(f1)>2)	{
		right <- left <- 0;
		prop_richness <- ttl_richness[f1]/sum(ttl_richness[f1]);
		f1cc <- 1;
		left <- f1[1]; 
		sum_prop <- prop_richness[1];
		while (sum_prop < 0.5)	{
			f1cc <- 1+f1cc;
			sum_prop <- sum_prop+prop_richness[f1cc];
			left <- c(left,f1);
			}
		right <- f1[!f1 %in% left];
		right <- right[order(-branch_ages[right])];
		left <- left[order(-branch_ages[left])];
		# shift rest of the tree away from ancestral node
		phy_pos[phy_pos<phy_pos[htu]] <- phy_pos[phy_pos<phy_pos[htu]]-sum(ttl_richness[right]);
		phy_pos[phy_pos>phy_pos[htu]] <- phy_pos[phy_pos>phy_pos[htu]]+sum(ttl_richness[left]);
		rr <- 1;
		phy_pos[right[rr]] <- phy_pos[htu]-ttl_richness[right[rr]];
		while (rr < length(right))	{
			rr <- rr+1;
			phy_pos[right[rr]] <- (phy_pos[right[rr-1]]-ttl_richness[right[rr-1]])-ttl_richness[right[rr]];
			}
		ll <- 1;
		phy_pos[left[ll]] <- phy_pos[htu] + ttl_richness[left[ll]];
		while (ll < length(left))	{
			ll <- ll+1;
			phy_pos[left[ll]] <- (phy_pos[left[ll-1]]+ttl_richness[left[ll-1]])+ttl_richness[left[ll]];
			}
		}	else if (length(f1)==2)	{
		f1 <- f1[order(-ttl_richness[f1])];
		if (phy_pos[htu]<phy_pos[1+notu])	{
			# going left is positive, so shift everything above this up by this amount
			if (last_right=="right")	{
				right <- f1[2];
				left <- f1[1];
				last_right <- "left"
				} else	{
				right <- f1[1];
				left <- f1[2];
				last_right <- "right"
				}
			} else	{
			# going right is negative, so shift everything below this down by this amount
			if (last_left=="left")	{
				right <- f1[1];
				left <- f1[2];
				last_left <- "right";
				} else	{
				right <- f1[2];
				left <- f1[1];
				last_left <- "left";
				}
			}
		# shift rest of the tree away from ancestral node
		phy_pos[phy_pos>phy_pos[htu]] <- phy_pos[phy_pos>phy_pos[htu]] + ttl_richness[left];
		phy_pos[phy_pos<phy_pos[htu]] <- phy_pos[phy_pos<phy_pos[htu]] - ttl_richness[right];
		phy_pos[left] <- phy_pos[htu] + ttl_richness[left];
		phy_pos[right] <- phy_pos[htu] - ttl_richness[right];
		}	else if (length(f1)==1)	{
		if (phy_pos[htu]<phy_pos[1+notu])	{
			if (last_right=="right")	{
				# going left is positive, so shift everything above this up by this amount
#				phy_y[phy_y>(phy_y[htu]+ttl_richness[f1])] <- phy_y[phy_y>phy_y[htu]] + ttl_richness[f1];
				# shift rest of the tree away from ancestral node
				phy_pos[phy_pos>phy_pos[htu]] <- phy_pos[phy_pos>phy_pos[htu]] + ttl_richness[f1];
				phy_pos[f1] <- phy_pos[htu] + ttl_richness[f1];
				last_right <- "left";
				} else {
				# going right is negative, so shift everything below this down by this amount
#				phy_y[phy_y<(phy_y[htu]-ttl_richness[f1])] <- phy_y[phy_y>ttl_richness[f1]] - ttl_richness[f1];
				# shift rest of the tree away from ancestral node
				phy_pos[phy_pos<phy_pos[htu]] <- phy_pos[phy_pos<phy_pos[htu]] - ttl_richness[f1];
				phy_pos[f1] <- phy_pos[htu] - ttl_richness[f1];
				last_right <- "right";
				}
			} else	{
			if (last_left=="right")	{
#				phy_y[phy_y>(phy_y[htu]+ttl_richness[f1])] <- phy_y[phy_y>phy_y[htu]] + ttl_richness[f1];
				# shift rest of the tree away from ancestral node
				phy_pos[phy_pos>phy_pos[htu]] <- phy_pos[phy_pos>phy_pos[htu]] + ttl_richness[f1];
				phy_pos[f1] <- phy_pos[htu] + ttl_richness[f1];
				last_left <- "left";
				} else	{
				# going right is negative, so shift everything below this down by this amount
#				phy_y[phy_y<(phy_y[htu]-ttl_richness[f1])] <- phy_y[phy_y>ttl_richness[f1]] - ttl_richness[f1];
				# shift rest of the tree away from ancestral node
				phy_pos[phy_pos<phy_pos[htu]] <- phy_pos[phy_pos<phy_pos[htu]] - ttl_richness[f1];
				phy_pos[f1] <- phy_pos[htu] - ttl_richness[f1];
				last_left <- "right";
				}
			}
		} 
	
	phy_pos <- phy_pos-phy_pos[notu+1];	# recenter around the base of the tree
#	print(c(nd,phy_pos));
	# now do species
	accounted <- c(accounted,mtree[nd,mtree[nd,]>0]);
	}
final_pos <- match(phy_pos,sort(unique(phy_pos)));
needed_edits <- hist(final_pos,breaks=0:max(final_pos),plot=F)$counts
too_many_here <- (1:length(needed_edits))[needed_edits>2];
while (length(too_many_here)>0) {
	ttl_minions <- 1:(notu+nNode);
	tmh <- 0;
	while (tmh < length(too_many_here))	{
		tmh <- tmh+1;
		problems <- ttl_minions[final_pos==too_many_here[tmh]];	# taxa overlapping each other
		these_nodes <- c();
		for (pp in 1:length(problems))
			these_nodes <- c(these_nodes,which(mtree==problems[pp],arr.ind = T)[1]);	# get the nodes containing problem cases;
		problem_ancestors <- problems[(notu+these_nodes) %in% problems];				# separate out sampled ancestors
		problem_ancestors_htu <- notu+these_nodes[match(problem_ancestors,problems)];	# keep track of the htu to which they belong, however!
		these_nodes <- these_nodes[!problems %in% problem_ancestors];
		problems <- problems[!problems %in% problem_ancestors];	# remove sampled ancestors for now
		starting_points <- final_pos[notu+these_nodes];			# positions of ancestral nodes/taxa
		adjust2 <- adjust <- starting_points-too_many_here[tmh];
		adjust2[adjust<0]<- -(length(adjust[adjust<0]):1);
		adjust2[adjust>0]<- 1:length(adjust[adjust>0]);
		final_pos[final_pos<too_many_here[tmh]] <- final_pos[final_pos<too_many_here[tmh]]+min(adjust2);
		final_pos[final_pos>too_many_here[tmh]] <- final_pos[final_pos>too_many_here[tmh]]+max(adjust2);
		phy_pos[phy_pos<phy_pos[problems[1]]] <- phy_pos[phy_pos<phy_pos[problems[1]]]+min(adjust2);
		phy_pos[phy_pos>phy_pos[problems[1]]] <- phy_pos[phy_pos>phy_pos[problems[1]]]+max(adjust2);
		final_pos[problems] <- final_pos[problems]+adjust;
		phy_pos[problems] <- phy_pos[problems]+adjust2;
		final_pos[problem_ancestors] <- final_pos[problem_ancestors_htu];
		phy_pos[problem_ancestors] <- phy_pos[problem_ancestors_htu];
		too_many_here <- too_many_here+max(adjust2);
		}
	final_pos <- match(phy_pos,sort(unique(phy_pos)));
	needed_edits <- hist(final_pos,breaks=0:max(final_pos),plot=F)$counts
	too_many_here <- (1:length(needed_edits))[needed_edits>2]
	}
return(final_pos);
}
# otu_cols <- rep("blue",notu)
# ctree <- vector tree in which each number gives the node from which a species evolved
# strat_ranges <- first and last appearance times of taxa
# durations <- originations and extinctions of taxa
# apos <- rep(1,length(divergence_times_1))
draw_calibrated_phylogeny_vertical <- function(ctree,strat_ranges,durations,apos,oldest=NA,youngest=NA,taxon_labels=NA,otu_cols,lazarus_col="gray50",branching_col="black")	{
# draws phylogeny onto an already configured plot
notu <- nrow(strat_ranges)
sctree <- cbind(rank(ctree),ctree)
otu_order <- order(ctree[1:notu])
sampled_ancestors <- accio_poss_ancestors_for_nodes(ctree,FA=strat_ranges[,1],apos)

phy_x <- vector(length=length(ctree))
if ((taxon_labels[1]!="numbers" && taxon_labels[1]!="Numbers")  && !is.na(taxon_labels[1]))	{
	y_adj <- (mxy-mny)/50
	x_adj <- -(notu+1)/37.5
	}
for (n in 1:notu)	{
	phy_x[n] <- match(n,otu_order)
	if (!is.na(match(n,sampled_ancestors)))	{
		segments(phy_x[n],durations[n,2],phy_x[n],strat_ranges[n,1],col=branching_col,lwd=4)
		segments(phy_x[n],durations[n,2],phy_x[n],strat_ranges[n,1],col=otu_cols[n],lwd=2)
		}
	segments(phy_x[n],durations[n,1],phy_x[n],durations[n,2],col=lazarus_col,lwd=4)
	if (strat_ranges[n,1]!=strat_ranges[n,2])	{
		rect((phy_x[n]-0.25),strat_ranges[n,1],(phy_x[n]+0.25),min(youngest,strat_ranges[n,2]),col=otu_cols[n])
		}	else	{
		segments((phy_x[n]-0.25),strat_ranges[n,1],(phy_x[n]+0.25),strat_ranges[n,1],lwd=4)
		segments((phy_x[n]-0.25),strat_ranges[n,1],(phy_x[n]+0.25),strat_ranges[n,1],lwd=2,col=otu_cols[n])
		}
	if (taxon_labels[1]=="numbers" || taxon_labels[1]=="Numbers")	{
		text(phy_x[n],strat_ranges[n,2],n,pos=3)
		} else if (!is.na(taxon_labels[1]))	{
		text(phy_x[n]+x_adj,strat_ranges[n,2]+y_adj,taxon_labels[n],srt=90,pos=4)
		}
	}
Nnode <- max(ctree) - notu
mtree <- accio_matrix_tree_from_vector_tree(ctree)
for (nn in Nnode:1)	{
	n <- notu+nn
	if (sampled_ancestors[n]==0)	{
		phy_x[n] <- mean(phy_x[mtree[nn,]])
		segments(min(phy_x[mtree[nn,]]),durations[n,2],max(phy_x[mtree[nn,]]),durations[n,2],lwd=1)
		segments(phy_x[n],durations[n,1],phy_x[n],durations[n,2],col=lazarus_col,lwd=4)
		}	else	{
		phy_x[n] <- phy_x[sampled_ancestors[n]]
		f1 <- mtree[nn,mtree[nn,]!=sampled_ancestors[nn]]
		for (f in 1:length(f1))	{
			segments(phy_x[f1[f]],durations[f1[f],1],phy_x[sampled_ancestors[n]],durations[f1[f],1],lwd=1)
			}
		segments(phy_x[n],durations[n,2],phy_x[n],durations[n,1],lwd=4,col=lazarus_col)
		}
#	nn <- nn-1
	}
}

#draw_calibrated_phylogeny_horizontal <- function(ctree,finds,durations,apomorphies,oldest=NA,youngest=NA,taxon_labels=NA,otu_cols,lazarus_col="gray50",branching_col="black",plot_stratigraphy="ranges",new_plot=F,xsize=4,ysize=6)	{
draw_calibrated_phylogeny_horizontal <- function(ctree,finds,durations,apomorphies,taxon_labels=NA,otu_cols,lazarus_col="gray50",branching_col="black",plot_stratigraphy="ranges",new_plot=F,xsize=4,ysize=6)	{
# draws phylogeny onto an already configured plot
notu <- match(-1,ctree)-1;
durations <- -abs(durations);
finds <- -abs(finds);
if (plot_stratigraphy=="ranges")	{
	strat_ranges <- data.frame(FAD=as.numeric(rep(0,notu)),LAD=as.numeric(rep(0,notu)),stringsAsFactors = F);
	for (n in 1:notu)	{
		if (sum(finds[n,]!=0)>0)	{
			strat_ranges$FAD[n] <- min(finds[n,]);
			strat_ranges$LAD[n] <- max(finds[n,finds[n,]!=0]);
			}
		}
	}

if (new_plot)	{
	mxx <- -abs(0.5*ceiling(max(durations)/0.5));
	mnx <- -abs(0.5*ceiling(min(durations)/0.5));
	par(pin=c(ysize,xsize));
	plot(NA,type='n',axes=FALSE,main="",xlab="",ylab="",xlim=c(mnx,mxx),ylim=c(1,notu))
	}

phy_y <- center_budding_phylogeny(ctree,durations,sampled_ancestors)
#for (n in 1:notu)	{
#	phy_y[n] <- match(n,atu_order);
#	if (!is.na(match(n,sampled_ancestors)))	{
#		segments(durations[n,2],phy_y[n],strat_ranges[n,1],phy_y[n],col=branching_col,lwd=3)
#		segments(durations[n,2],phy_y[n],strat_ranges[n,1],phy_y[n],col=otu_cols[n],lwd=1.5)
#		}
adj_y <- (mxx-mnx)*0.005;
nNode <- max(ctree) - notu;
for (nn in nNode:1)	{
	n <- notu+nn;
	f1 <- mtree[nn,!mtree[nn,] %in% c(0,sampled_ancestors[n])];
	for (f in 1:length(f1))
		segments(durations[f1[f],1]-adj_y,phy_y[f1[f]],durations[f1[f],1]-adj_y,phy_y[n],lwd=1);
	nn <- nn-1;
	}

n <- 0;
while (n < notu)	{
	n <- n+1;
	segments(durations[n,1],phy_y[n],durations[n,2],phy_y[n],col=lazarus_col,lwd=3);
	if (plot_stratigraphy=="ranges")	{
		if (strat_ranges[n,1]!=strat_ranges[n,2])	{
			rect(strat_ranges[n,1],(phy_y[n]-0.25),min(youngest,strat_ranges[n,2]),(phy_y[n]+0.25),col=otu_cols[n])
			}	else if (strat_ranges[n,1]!=0)	{
#			segments(strat_ranges[n,1],(phy_y[n]-0.25),strat_ranges[n,1],(phy_y[n]+0.25),lwd=4)
#			segments(strat_ranges[n,1],(phy_y[n]-0.25),strat_ranges[n,1],(phy_y[n]+0.25),lwd=2,col=otu_cols[n])
			points(strat_ranges[n,1],phy_y[n],pch=22,cex=1,bg=otu_cols[n]);
			} 
		} else if (plot_stratigraphy=="points")	{
		these_finds <- finds[n,finds[n,]!=0];
		tf <- 0;
		while (tf < length(these_finds))	{
			tf <- 1+tf;
			points(these_finds[tf],phy_y[n],pch=21,cex=1,bg=otu_cols[n])
			}
		}
	if (!is.na(taxon_labels))	{
		if (taxon_labels[1]=="numbers" || taxon_labels[1]=="Numbers")	{
			text(phy_y[n],durations[n,2],n,pos=3)
			} else if (!is.na(taxon_labels[1]))	{
			text(phy_y[n]+x_adj,durations[n,2]+y_adj,taxon_labels[n],srt=90,pos=4)
			}
		}
	}
return(atu_order);
}

draw_calibrated_phylogeny_flex <- function(ctree,finds,durations,apos,orientation="vertical",oldest=NA,youngest=NA,taxon_labels=NA,otu_cols,lazarus_col="gray50",branching_col="black",plot_stratigraphy="ranges",new_plot=F,xsize=4,ysize=6)	{
# draws phylogeny onto an already configured plot
notu <- match(-1,ctree)-1;
durations <- -abs(durations);
finds <- -abs(finds);
if (plot_stratigraphy=="ranges")	{
	strat_ranges <- data.frame(FAD=as.numeric(rep(0,notu)),LAD=as.numeric(rep(0,notu)),stringsAsFactors = F);
	for (n in 1:notu)	{
		if (sum(finds[n,]!=0)>0)	{
			strat_ranges$FAD[n] <- min(finds[n,]);
			strat_ranges$LAD[n] <- max(finds[n,finds[n,]!=0]);
			}
		}
	}

phy_z <- center_budding_phylogeny(ctree,durations,sampled_ancestors);
if (new_plot)	{
	mxz <- -abs(0.5*ceiling(max(durations)/0.5));
	mnz <- -abs(0.5*ceiling(min(durations)/0.5));
	
	if (orientation=="vertical")	{
		par(pin=c(min(ysize,xsize),max(ysize,xsize)));
		plot(NA,type='n',axes=FALSE,main="",xlab="",ylab="",ylim=c(mnz,mxz),xlim=c(1,max(phy_z)));
		plot(NA,type='n',axes=T,main="",xlab="",ylab="",ylim=c(mnz,mxz),xlim=c(1,max(phy_z)));
		} else	{
		par(pin=c(max(ysize,xsize),min(ysize,xsize)));
		plot(NA,type='n',axes=FALSE,main="",xlab="",ylab="",xlim=c(mnz,mxz),ylim=c(1,max(phy_z)));
		plot(NA,type='n',axes=T,main="",xlab="",ylab="",xlim=c(mnz,mxz),ylim=c(1,max(phy_z)));
		}
	}

adj_z <- (mxz-mnz)*0.0025;
nNode <- max(ctree) - notu;
for (nn in nNode:1)	{
	n <- notu+nn;
	f1 <- mtree[nn,!mtree[nn,] %in% c(0,sampled_ancestors[n])];
	for (f in 1:length(f1))	{
		if (orientation=="horizontal")	{
			segments(durations[f1[f],1]-adj_z,phy_z[f1[f]],durations[f1[f],1]-adj_z,phy_z[n],lwd=1);
			} else	{
			segments(phy_z[f1[f]],durations[f1[f],1]-adj_z,phy_z[n],durations[f1[f],1]-adj_z,lwd=1);
			}
		}
	nn <- nn-1;
	}

n <- 0;
while (n < notu)	{
	n <- n+1;
	if (orientation=="horizontal")	{
		segments(durations[n,1],phy_z[n],durations[n,2],phy_z[n],col=lazarus_col,lwd=3);
		} else	{
		segments(phy_z[n],durations[n,1],phy_z[n],durations[n,2],col=lazarus_col,lwd=3);
		}
	if (plot_stratigraphy=="ranges")	{
		if (strat_ranges[n,1]!=strat_ranges[n,2])	{
			rect(strat_ranges[n,1],(phy_z[n]-0.25),min(youngest,strat_ranges[n,2]),(phy_z[n]+0.25),col=otu_cols[n])
			}	else if (strat_ranges[n,1]!=0)	{
#			segments(strat_ranges[n,1],(phy_y[n]-0.25),strat_ranges[n,1],(phy_y[n]+0.25),lwd=4)
#			segments(strat_ranges[n,1],(phy_y[n]-0.25),strat_ranges[n,1],(phy_y[n]+0.25),lwd=2,col=otu_cols[n])
			points(strat_ranges[n,1],phy_z[n],pch=22,cex=1,bg=otu_cols[n]);
			} 
		} else if (plot_stratigraphy=="points")	{
		these_finds <- finds[n,finds[n,]!=0];
		tf <- 0;
		while (tf < length(these_finds))	{
			tf <- 1+tf;
			points(these_finds[tf],phy_z[n],pch=21,cex=1,bg=otu_cols[n])
			}
		}
	if (!is.na(taxon_labels))	{
		if (taxon_labels[1]=="numbers" || taxon_labels[1]=="Numbers")	{
			text(phy_z[n],durations[n,2],n,pos=3)
			} else if (!is.na(taxon_labels[1]))	{
			text(phy_z[n]+x_adj,durations[n,2]+y_adj,taxon_labels[n],srt=90,pos=4)
			}
		}
	}
#	if (sampled_ancestors[n]==0)	{
##		phy_y[n] <- mean(phy_y[mtree[nn,]])
#		segments(durations[n,2],min(phy_y[mtree[nn,]]),durations[n,2],max(phy_y[mtree[nn,]]),lwd=1)
#		segments(durations[n,1],phy_y[n],durations[n,2],phy_y[n],col=lazarus_col,lwd=4)
#		}	else	{
#		phy_y[n] <- phy_y[sampled_ancestors[n]];
		
#		segments(durations[n,2],phy_y[n],durations[n,1],phy_y[n],lwd=4,col=lazarus_col)
return(atu_order);
}

plot_calibrated_phylogeny <- function(ctree,strat_ranges,durations,apos,oldest=NA,youngest=NA,taxon_labels=NA,otu_cols,xsize=3.5,ysize=3.5)	{
notu <- dim(strat_ranges)[1]
sctree <- cbind(rank(ctree),ctree)
ordinate <- "Ma"
otu_order <- order(ctree[1:notu])
sampled_ancestors <- accio_poss_ancestors_for_nodes(ctree,FA=strat_ranges[,1],apos)

# set up Y-axis (time) scale
if (is.na(oldest))	{
	ax_min <- min(durations)
	}	else	{
	ax_min <- min(oldest,min(durations))
	}
if (is.na(youngest))	{
	ax_max <- max(strat_ranges)
	} else	{
	ax_max <- max(youngest,max(strat_ranges))
	}
if (!is.na(taxon_labels[1]))	{
	y_add <- TRUE
	}	else	{
	y_add <- FALSE
	}
axis_info <- wagner_set_axes(ax_min,ax_max,y_add)
tcs <- axis_info$Ticks
lbl_prn <- axis_info$Labels
tick_str <- axis_info$Tick_Strength
par(pin=c(xsize,ysize))
mny <- min(tcs)
mxy <- max(tcs[dim(tcs)[1],tcs[dim(tcs)[1],]!=0])
plot(NA,type='n',axes=FALSE,xlab="",ylab=ordinate,xlim=c(0,notu+1),ylim=c(mny,mxy))
axis(2,at=seq(mny,mxy,by=(abs(mxy-mny))),tcl=-0.00,labels=FALSE,lwd=1.1,lwd.ticks=0.0,las=2)
for (i in 1:dim(tcs)[1])	{
	up <- match(max(tcs[i,]),tcs[i,])
	axis(2,at=tcs[i,1:up],tcl=tick_str[i],labels=lbl_prn[i,1:up],lwd=0.0,lwd.ticks=1.1,las=2)
	}
phy_x <- vector(length=length(ctree))
if ((taxon_labels[1]!="numbers" && taxon_labels[1]!="Numbers")  && !is.na(taxon_labels[1]))	{
	y_adj <- (mxy-mny)/50
	x_adj <- -(notu+1)/37.5
	}
for (n in 1:notu)	{
	phy_x[n] <- match(n,otu_order)
	if (!is.na(match(n,sampled_ancestors)))	{
		segments(phy_x[n],durations[n,2],phy_x[n],strat_ranges[n,1],col="black",lwd=4)
		segments(phy_x[n],durations[n,2],phy_x[n],strat_ranges[n,1],col=otu_cols[n],lwd=2)
		}
	segments(phy_x[n],durations[n,1],phy_x[n],durations[n,2],col="gray50",lwd=4)
	if (strat_ranges[n,1]!=strat_ranges[n,2])	{
		rect((phy_x[n]-0.25),strat_ranges[n,1],(phy_x[n]+0.25),strat_ranges[n,2],col=otu_cols[n])
		}	else	{
		segments((phy_x[n]-0.25),strat_ranges[n,1],(phy_x[n]+0.25),strat_ranges[n,1],lwd=4)
		segments((phy_x[n]-0.25),strat_ranges[n,1],(phy_x[n]+0.25),strat_ranges[n,1],lwd=2,col=otu_cols[n])
		}
	if (taxon_labels[1]=="numbers" || taxon_labels[1]=="Numbers")	{
		text(phy_x[n],strat_ranges[n,2],n,pos=3)
		} else if (!is.na(taxon_labels[1]))	{
		text(phy_x[n]+x_adj,strat_ranges[n,2]+y_adj,taxon_labels[n],srt=90,pos=4)
		}
	}
Nnode <- max(ctree) - notu
mtree <- accio_matrix_tree_from_vector_tree(ctree)
for (nn in Nnode:1)	{
	n <- notu+nn
	if (sampled_ancestors[n]==0)	{
		phy_x[n] <- mean(phy_x[mtree[nn,]])
		segments(min(phy_x[mtree[nn,]]),durations[n,2],max(phy_x[mtree[nn,]]),durations[n,2],lwd=1)
		segments(phy_x[n],durations[n,1],phy_x[n],durations[n,2],col="gray50",lwd=4)
		}	else	{
		phy_x[n] <- phy_x[sampled_ancestors[n]]
		f1 <- mtree[nn,mtree[nn,]!=sampled_ancestors[nn]]
		for (f in 1:length(f1))	{
			segments(phy_x[f1[f]],durations[f1[f],1],phy_x[sampled_ancestors[n]],durations[f1[f],1],lwd=1)
			}
		segments(phy_x[n],durations[n,2],phy_x[n],durations[n,1],lwd=4,col="gray50")
		}
#	nn <- nn-1
	}
}

draw_cladogram_from_newick_string <- function(newick_string)	{
### NOTE: This assumes that each numbered taxon is in the cladogram
###	  IF you have only SOME taxa in the string, then this needs too be rewritten
otu_order_string <- simplify2array(strsplit(newick_string,split=",")[1])[,1]
otu_order_string <- gsub("\\(","",otu_order_string);
otu_order_string <- gsub("\\)","",otu_order_string);
otu_order_string <- as.numeric(gsub(";","",otu_order_string));
otu_order <- match(1:max(otu_order_string),otu_order_string);

cladogram <- read_newick_string(newick_string);
notu <- length(need_to_be_drawn);
tttu <- length(cladogram);
nhtu <- tttu-notu;
otu_order <- c(otu_order,rep(0,nhtu));
undrawn <- 1:nhtu;

mxy <- mxx <- notu+1;
mnx <- 0;
mny <- -mxy;
par(pin=c(4.5,4.5));

lineage_tops <- array(0,dim=tttu);

node_span <- array(0,dim=c(tttu,2));
node_span[1:notu,1] <- node_span[1:notu,2] <- otu_order[1:notu];
plot(NA,type='n',axes=FALSE,main="",xlab="",ylab="",xlim=c(mnx,mxx),ylim=c(mny,mxy));
need_to_be_drawn <- sort(otu_order_string);
for (nn in nhtu:1)	{
	node <- nn+notu;
	descendants <- (1:tttu)[cladogram==(notu+nn)];
	if (sum(descendants %in% need_to_be_drawn)==length(descendants))	{
		node_span[node,1] <- min(node_span[descendants,1])
		node_span[node,2] <- max(node_span[descendants,2])
		otu_order[node] <- mean(node_span[node,]);
		lineage_tops[node] <- -(node_span[node,2]-node_span[node,1]);
		for (dd in 1:length(descendants))	{
			segments(otu_order[descendants[dd]],lineage_tops[descendants[dd]],otu_order[node],lineage_tops[node],lwd=4);
			}
		need_to_be_drawn <- need_to_be_drawn[!need_to_be_drawn %in% descendants];
		need_to_be_drawn <- sort(c(need_to_be_drawn,node));
		}
	nn <- nn-1;

	}
}