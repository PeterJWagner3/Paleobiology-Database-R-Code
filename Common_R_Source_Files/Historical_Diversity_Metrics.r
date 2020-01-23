MAXNO <- 1.797693e+308;
# Routines to estimate different aspects of historical diversity patterns of the sort described
#	in Paleobiology all of the time!

#source('~/Documents/R_Projects/Common_R_Source_Files/Wagner_Stats_and_Probability_101.r');
#source('~/Documents/R_Projects/Common_R_Source_Files/Wagner_kluges.r')  #

# calculate center of gravity for continuous data (e.g., disparity)
Center_of_Gravity_continuous <- function(metric, start, end)	{
# metric: vector giving the data at any point in time
CS <- vector(length=4)
#SS <- vector(length=3)
c0 <- c1 <- c2 <- c3 <- sm <- 0
sm <- sum(metric)
SS <- c()
for (tt in start:end)	{
	ttt <- 1+tt-start
	if (tt==start)	{
#		SS[tt] <- metric[tt]/sm
		SS <- metric[tt]/sm
		}	else	{
#		SS[tt] <- SS[tt-1]+metric[tt]/sm
		SS <- c(SS,SS[ttt-1]+metric[tt]/sm)
		}
	c0 <- c0+metric[tt]
	c1 <- c1+(metric[tt]*(ttt-0.5))
	c2 <- c2+(metric[tt]*(ttt-0.5)*(ttt-0.5))
	c3 <- c3+(metric[tt]*(ttt-0.5)*(ttt-0.5)*(ttt-0.5))
		#	c4 <- c4+(metric[chrons]*(chrons+0.5)*(chrons+0.5)*(chrons+0.5)*(chrons+0.5))
	}
CS[1] <- (c1/c0)/(1+end-start)	# center of gravity	*/
A <- ((c0*c3)-(3*c1*c2)+(2*c1*c1*c1/c0))
z <- c0*c0*(((c2-(c1*c1/c0))/c0)^1.5)	# asymmetry	*/
	#CS[1] <- A <- ((c0*c3)-(3*c1*c2) + (2*c1*c1*c1)/c0)/(c0*c0)*(((c2-((c1*c1)/c0))/c0)^1.5)
if (A>0)	CS[2] <- A/z
	
ttt <- 1
while (SS[ttt]<0.5)	ttt <- ttt+1
if (ttt>1)	{
	x <- CS[2] <- (0.5-SS[ttt-1])
	y <- (SS[ttt]-SS[ttt-1])
	x <- x/y
	x <- x+(ttt-1)
	}	else	{
	x <- 0.5/SS[1]
	}
CS[3] <- x/(1+end-start)		# where midpoint is reached	*/
CS[4] <- sm/(1+end-start)
names(CS) <- c("Center_of_Gravity","Asymmetry","Midpoint","Weight")
return(CS)
}

convert_time_scale_to_interval_durations <- function(timescale)	{
bins <- length(timescale)-1
dt <- vector(length=bins)
for (i in 1:bins)	dt[i] <- abs(timescale[i]-timescale[i+1])
return(dt)
}

freqrat <- function(taxon_ranges)	{
if (!is.array(taxon_ranges))	{
	counts <- hist(plot=FALSE,taxon_ranges,breaks=0:max(taxon_ranges))$counts
	} else	{
	tx_ranges <- 1+abs(taxon_ranges[,2]-taxon_ranges[,1])
	}
return((counts[2]^2)/(counts[1]*counts[3]))
}

ML_FreqRat_Rate_Support_Bars <- function(ages, end, final, support_width=1.0, precision=0.01)	{
# span will be the oldest considered.  For survivorship, that will be the recent.  
#	All survivors should be put there so the exponential is thrown off.
#	For prenascents, that will be the first interval. Instead of assuming 
#	everything originates there, start from there, stop one before.  This means no
#	origination rates before the 2nd interval
if (final || final==1)	{
	span <- end
	}	else span <- (end-1)
rates <- (1/precision)-1
meth <- 1+max(ages)
#brks <- vector(length=meth)
#for (i in 1:meth)	brks[i] <- i-1
#lifespans <- vector(length=bins)
bb <- count(ages)
lives <- vector(length=max(bb[,1]))
for (i in 1:max(bb[,1]))	lives[bb[i,1]] <- bb[i,2]
#lives <- c(240,60,30,15,8,4,2,1)
freqratsson <- matrix(0,rates,rates)
for (chrons in 1:rates)	{
	exp_ranges <- vector(length=bins)
	turn <- chrons*precision
	# get expected lifetimes
	for (b in 1:(bins-1))	{
		exp_ranges[b] <- turn*((1-turn)^(b-1))
		}
	exp_ranges[bins] <- 1-sum(exp_ranges[1:(bins-1)])
	if (exp_ranges[bins]==0)	exp_ranges[bins] <- turn*((1-turn)^(bins-1))
	
	for (f in 1:rates)	{
		fr <- f*precision
		# get expected ranges given chrons & fr
		exp_finds <- vector(length=bins)
		for (a in 1:bins)	{
			# multiply proportion of taxa living that long times proportion expected to be sampled for certain ranges
			for (r in 1:a)	{
				if (r==1)	{
					# probablity of living 1 stage & being found in it
					if (a==1)	{
						exp_finds[r] <- exp_ranges[a]*fr
						}	else {
						# probablity of living 2+ stages & being found in only in the first
						exp_finds[r] <- exp_finds[r]+exp_ranges[a]*fr*((1-fr)^(a-r))
						}
					}	else {
					# probablity of being found in this stage and one subsequent one & none afterwards
					exp_finds[r] <- exp_finds[r]+exp_ranges[a]*(fr^2)*((1-fr)^(a-r))
					}
				}	
			}
		exfn <- exp_finds[1:(meth-1)]/sum(exp_finds)
		freqratsson[chrons,f] <- sum(lives*log(exfn))
		}
	}
freqratsson <- freqratsson-max(freqratsson)
turn <- vector(length=rates)
pres <- vector(length=rates)
for (r in 1:rates)	{
	turn[r] <- sum(exp(freqratsson[r,]))
	pres[r] <- sum(exp(freqratsson[,r]))
	}
turn <- turn/sum(turn)
pres <- pres/sum(pres)
turn <- log(turn)
pres <- log(pres)
turn <- max(turn)-turn
pres <- max(pres)-pres
output <- matrix(0,2,3)
rownames(output) <- c("turnover","preservation")
colnames(output) <- c("ml","lb","ub")
output[1,1] <- precision*match(0,turn)
output[2,1] <- precision*match(0,pres)
output[,3] <- 0.995
for (r in 2:rates)	{
	if (turn[r-1]>support_width && turn[r]<support_width)	{
		output[1,2] <- precision*r
		} else if (turn[r-1]<support_width && turn[r]>support_width)	{
		output[1,3] <- precision*r
		}
	if (pres[r-1]>support_width && pres[r]<support_width)	{
		output[2,2] <- precision*r
		} else if (pres[r-1]<support_width && pres[r]>support_width)	{
		output[2,3] <- precision*r
		}
	}
return(output)
}

ML_FreqRat_with_ranges <- function(lives)	{
# span will be the oldest considered.  For survivorship, that will be the recent.  
#	All survivors should be put there so the exponential is thrown off.
#	For prenascents, that will be the first interval. Instead of assuming 
#	everything originates there, start from there, stop one before.  This means no
#	origination rates before the 2nd interval
#if (final || final==1)	{
#	span <- bins
#	}	else spans <- bins-1
#lives <- c(240,60,30,15,8,4,2,1)
bins <- length(lives)
freqratsson <- matrix(0,99,99)
for (chrons in 1:99)	{
	exp_ranges <- vector(length=bins)
	turn <- chrons/100
	# get expected lifetimes
	for (b in 1:(bins-1))	{
		exp_ranges[b] <- turn*((1-turn)^(b-1))
		}
	exp_ranges[bins] <- 1-sum(exp_ranges[1:(bins-1)])
	if (exp_ranges[bins]==0)	exp_ranges[bins] <- turn*((1-turn)^(bins-1))
	
	for (f in 1:99)	{
		fr <- f/100
		# get expected ranges given chrons & fr
		exp_finds <- vector(length=bins)
		for (a in 1:bins)	{
			# multiply proportion of taxa living that long times proportion expected to be sampled for certain ranges
			for (r in 1:a)	{
				if (r==1)	{
					# probablity of living 1 stage & being found in it
					if (a==1)	{
						exp_finds[r] <- exp_ranges[a]*fr
						}	else {
						# probablity of living 2+ stages & being found in only in the first
						exp_finds[r] <- exp_finds[r]+exp_ranges[a]*fr*((1-fr)^(a-r))
						}
					}	else {
					# probablity of being found in this stage and one subsequent one & none afterwards
					exp_finds[r] <- exp_finds[r]+exp_ranges[a]*(fr^2)*((1-fr)^(a-r))
					}
				}	
			}
		exfn <- exp_finds[1:(meth-1)]/sum(exp_finds)
		freqratsson[chrons,f] <- sum(lives*log(exfn))
		}
	}
freqratsson <- freqratsson-max(freqratsson)
return(freqratsson)
}

ML_FreqRat <- function(synoptic, interval_scale)	{
# span will be the oldest considered.  For survivorship, that will be the recent.  
#	All survivors should be put there so the exponential is thrown off.
#	For prenascents, that will be the first interval. Instead of assuming 
#	everything originates there, start from there, stop one before.  This means no
#	origination rates before the 2nd interval
#if (final || final==1)	{
#	span <- bins
#	}	else spans <- bins-1
bins <- length(interval_scale)
notu <- dim(synoptic)[1]
ages <- vector(length=notu)
interval_scale_up <- c(interval_scale[2:length(interval_scale)],0)
for (n in 1:notu)	{
	ages[n] <- sum((interval_scale>synoptic[n,2])*(interval_scale_up<=synoptic[n,1]))
	}
meth <- 1+max(ages)
brks <- 0:max(ages)
#lifespans <- vector(length=ages)
lives <- hist(ages,breaks=brks)$counts
#lives <- c(240,60,30,15,8,4,2,1)
freqratsson <- matrix(0,99,99)
for (chrons in 1:99)	{
	exp_ranges <- vector(length=bins)
	turn <- chrons/100
	# get expected lifetimes
	for (b in 1:(bins-1))	{
		exp_ranges[b] <- turn*((1-turn)^(b-1))
		}
	exp_ranges[bins] <- 1-sum(exp_ranges[1:(bins-1)])
	if (exp_ranges[bins]==0)	exp_ranges[bins] <- turn*((1-turn)^(bins-1))
	
	for (f in 1:99)	{
		fr <- f/100
		# get expected ranges given chrons & fr
		exp_finds <- vector(length=bins)
		for (a in 1:bins)	{
			# multiply proportion of taxa living that long times proportion expected to be sampled for certain ranges
			for (r in 1:a)	{
				if (r==1)	{
					# probablity of living 1 stage & being found in it
					if (a==1)	{
						exp_finds[r] <- exp_ranges[a]*fr
						}	else {
						# probablity of living 2+ stages & being found in only in the first
						exp_finds[r] <- exp_finds[r]+exp_ranges[a]*fr*((1-fr)^(a-r))
						}
					}	else {
					# probablity of being found in this stage and one subsequent one & none afterwards
					exp_finds[r] <- exp_finds[r]+exp_ranges[a]*(fr^2)*((1-fr)^(a-r))
					}
				}	
			}
		exfn <- exp_finds[1:(meth-1)]/sum(exp_finds)
		freqratsson[chrons,f] <- sum(lives*log(exfn))
		}
	}
freqratsson <- freqratsson-max(freqratsson)
return(freqratsson)
}

expected_origination_given_logistic_constant_extinction <- function(mu,R,K,S=1)	{
# mu: extinction
# R: intrinsic rate of diversification (drops as S increases)
# K: "carrying capacity" richness.
# S: initial richness (1 unless otherwise specified)
pp <- vector(length=1)
dS <- R*S*(1-(S/K))
pp[1] <- mu+(log(S+dS)-log(S))
while (min(pp)>mu)	{
#	S <- S*exp((pp[chrons]-mu)*dt)
	S <- min(K,(S+(R*S*(1-(S/K)))))
	pp <- c(pp,mu+log(S+(R*S*(1-(S/K))))-log(S))
#	print(c(length(pp),S,pp[length(pp)]))
	}
return(pp)
}

Per_Bin_Sampling_Given_Range_Throughs <- function(ranges, finds)	{
genera <- dim(ranges)[1]
bins <- dim(finds)[1]
total_range_throughs <- vector(length=bins)
unsampled_range_throughs <- vector(length=bins)
for (g in 1:genera) {
    if (ranges[g,1]<(ranges[g,2]-1) && ranges[g,1]!=(bins+1))  {
        for (b in (ranges[g,1]+1):(ranges[g,2]-1))  {
            total_range_throughs[b] <- total_range_throughs[b]+1
            if (finds[b,g]==0)
                unsampled_range_throughs[b] <- unsampled_range_throughs[b]+1
            }
        }
    }
preservation <- vector(length=bins)
for (b in 1:bins)   {
    if (total_range_throughs[b]>0)  preservation[b] <- 1-(unsampled_range_throughs[b]/total_range_throughs[b])
    }
}

Per_Taxon_Sampling_Given_Range_Throughs <- function(ranges, finds)	{
genera <- dim(ranges)[1]
bins <- dim(finds)[1]
total_range_throughs <- vector(length=bins)
unsampled_range_throughs <- vector(length=bins)
for (g in 1:genera) {
    if (ranges[g,1]<(ranges[g,2]-1) && ranges[g,1]!=(bins+1))  {
        for (b in (ranges[g,1]+1):(ranges[g,2]-1))  {
            total_range_throughs[b] <- total_range_throughs[b]+1
            if (finds[b,g]==0)
                unsampled_range_throughs[b] <- unsampled_range_throughs[b]+1
            }
        }
    }

preservation <- 1-(sum(unsampled_range_throughs)/sum(total_range_throughs))
}

### Routines for Phylogeny Prior Probabilities
prob_sampling_clade_bapst <- function(p,q,r)	{
# modified from Bapst (2013) equation 2
# p: origination (lambda in Bapst 2013)
# q: extinction  (mu in Bapst 2013)
# r: sampling (psi in Bapst 2013)
s <- 1		# richness of clade (K in Bapst's paper)
pmiss <- mxp <- 0
num <- den <- 1
while ((mxp/100000)<(num/den))	{
	num <- (p^(s-1))*(q^s)*choose((2*s)-2,s-1)
	den <- s*((p+q+r)^((2*s)-1))
	pmiss <- pmiss+(num/den)
	if (mxp<(num/den))	mxp <- (num/den)
#	print(c(s,num,den,num/den,pmiss))
	s <- s+1
	}
return(1-pmiss)
}

prob_missing_clade_bapst <- function(p,q,r)	{
# modified from Bapst (2013) equation 2
# p: origination (lambda in Bapst 2013)
# q: extinction  (mu in Bapst 2013)
# r: sampling (psi in Bapst 2013)
s <- 1		# richness of clade (K in Bapst's paper)
pmiss <- mxp <- 0
num <- den <- 1
while ((mxp/100000)<(num/den))	{
	num <- (p^(s-1))*(q^s)*choose((2*s)-2,s-1)
	den <- s*((p+q+r)^((2*s)-1))
	pmiss <- pmiss+(num/den)
	if (mxp<(num/den))	mxp <- (num/den)
#	print(c(s,num,den,num/den,pmiss))
	s <- s+1
	pmiss
	}
return(pmiss)
}

Prob_Gap_Given_Parameters <- function(gap,origination,extinction,sampling)	{
prob_samp <- prob_sampling_clade_bapst(origination,extinction,sampling)
return (((sampling+(prob_samp*origination))^3)*(gap^2)*exp((sampling+(prob_samp*origination))*gap*-1)/2)
}

prob_missing_ancestor_from_node <- function(extinction,sampling)	{
return(dpois(0,0.5*sampling/extinction))
}

expected_gap_to_node <- function(origination,extinction,sampling)	{
pmiss <- 1-prob_sampling_clade_bapst(origination,extinction,sampling)
amiss <- prob_missing_ancestor_from_node(extinction,sampling)
gap <- origination*(1-(pmiss*amiss))+sampling
return(gap)
}

prob_ghost_lineage_duration <- function(p,q,r,chrons)	{
# THIS IS THE GOOD ONE
pfind <- prob_sampling_clade_bapst(p,q,r)	# P clades is found
pgl <- exp(-(r+(pfind*p))*chrons)		# P no finds and no sampled sister taxa
	#pgl <- dpois(1,(pfind*p)*2*chrons)*dpois(0,r*chrons)
return(pgl)
}

prob_ghost_taxon_duration <- function(p,q,r,chrons)	{
# THIS IS THE GOOD ONE
pfind <- prob_sampling_clade_bapst(p,q,r)	# P clades is found
#pnd <- 1-prob_amalgamate_at_node(p,q,r,chrons)		# P[have not encountered another node]
psab <- exp(-p*chrons)								# P[still direct ancestor at td]
psap <- prob_sampling_taxon(q,r)				# P[sample ancestor after
pgp <- exp(-r*chrons)								# P[unsampled over tg=td]
punsmp <- 1-(psab*psap)						# P[sampled not filling gap with ancestor]
pgt <- punsmp*exp(-(r+(pfind*p))*chrons)
#pgt <- (pgp*pfind)*punsmp
return(pgt)
}

assign_interval_rates_to_finer_time_scale <- function(bin_rates,bin_spans,prec)	{
#slices <- 1+(max(bin_spans)-min(bin_spans))/prec
bins <- length(bin_rates)
fine_rates <- c()
#db_cells <- c()
for (b in 1:bins)	{
	cells <- round(bin_spans[b]/prec,0)
	fine_rates <- c(fine_rates,rep(bin_rates[b],cells))
#	db_cells <- c(db_cells,cells)
#	added_rates <- rep(bin_rates[b],cells)
#	if (b==1)	{
#		fine_rates <- added_rates
#		}	else {
#		fine_rates <- c(fine_rates,added_rates)
#		}
	}
#slices <- (sum(bin_spans)/prec)
#tdz <- base::cumsum(bin_spans)/prec
#tda <- c(1,tdz[1:(bins-1)]+1)
#fine_rates <- vector(length=slices)
# for unknown reasons, this does not work!  it leaves 0's at #'s divisible by 256....
#for (b in 1:bins)	for (c in tda[b]:tdz[b])	fine_rates[c] <- bin_rates[b]
return(fine_rates)
}

# written 2017-01-06
# assigns all expected origination to first interval or all expected
#	extinction to last interval
make_rates_pulsed_at_finer_scale <- function(bin_rates,bin_spans,pulse_type,prec)	{
# bin_rates: expected events per bin
# bin_spans: span of bin in time
# pulse_type: 1 for origination, 2 for extinction
# prec: how finely time is divided; if 0.1, then 10 cells for every 1 in bin_spans	
bins <- length(bin_rates)
for (b in 1:bins)	{
	pulse <- bin_rates[b]/prec
	cells <- round(bin_spans[b]/prec,0)-1
	if (pulse_type==1)	{
		if (b==1)	{
			fine_rates <- c(pulse,rep(0,cells))
			}	else	{
			fine_rates <- c(fine_rates,c(pulse,rep(0,cells)))
			}
		} else	{
		if (b==1)	{
			fine_rates <- c(rep(0,cells),pulse)
			}	else	{
			fine_rates <- c(fine_rates,c(rep(0,cells),pulse))
			}
		}
	}
return(fine_rates)
}

# written 2017-01-01
accio_prob_sampling_paraclade_of_unknown_size_from_interval <- function(p,q,r,chrons,continuous=TRUE)	{
# p: origination rate during interval
# q: extinction rate during interval
# r: sampling rate for interval
# chrons: duration of interval
# routine integrates over Raup's equations for P[richness=N | pp,qq,dur]
	# this gives expected amount of time with N species
# ∑P[1+ find | r*N*chrons] x P[N | chrons,p,q] over all values of N creating meaningful numbers
mx_exp_n <- exp_n <- integrate(prob_paraclade_extinction,lower=0,upper=chrons,a=1,p=p,q=q)$value
n<-1
while (exp_n >= (mx_exp_n/10000000))	{
	exp_n <- integrate(prob_n_species_at_time_tm,lower=0,upper=chrons,n=n,a=1,p=p,q=q)$value
	if (mx_exp_n < exp_n)	mx_exp_n <- exp_n
	if (n==1)	{
		exp_outs <- exp_n
		}	else	{
		exp_outs <- c(exp_outs,exp_n)
		}
	n <- n+1
	}
exp_outs_cond <- exp_outs * (1:length(exp_outs))
p_samp <- 1-dpois(0,r*sum(exp_outs_cond))
return(p_samp)
#prec <- 0.01
#chrons<-seq(prec,3,by=prec)
#pns <- sapply(chrons,prob_n_species_at_time_tm,n=0,a=1,p=p,q=q)
#sum(pns)
#(integrate(prob_n_at_ti_extinct_at_end,lower=0,upper=ti,n=n,p=p,q=q,chrons=chrons)$value)/pdodo	
}

# routine written during strophomenoid project 2017-06-14 to figure out how the hell I found phi
accio_paraclade_sampling_per_time_slice <- function(p,q,r,ts,prec) {
# routine to get phis at different points in time
# returns a vector of length = length of 10xtimescale duration
# ts: time scale giving onset of each interval.  If there are 12 intervals, then there should
#	be 13 intervals.
#bins <- length(p)
bs <- c()
for (i in 1:(length(ts)-1)) bs <- c(bs,abs (ts[i+1]-ts[i]))
slices <- sum(bs)/prec
pp <- assign_interval_rates_to_finer_time_scale(bin_rates=p,bs,prec)
qq <- assign_interval_rates_to_finer_time_scale(bin_rates=q,bs,prec)
rr <- assign_interval_rates_to_finer_time_scale(bin_rates=r,bs,prec)
dt <- rep(prec,slices)

# phi for individual time slice
phi1s <- accio_p_sample_lineage_present_over_interval_continuous(pp,qq,rr,dt)
return(phi1s)
}
# pp: per-interval origination rates
# qq: per-interval extinction rates
# rr: per-interval sampling rates
# ts: time-scale, with first interval the oldest and in negative numbers
# prec: how finely you want the timescale, with 1 = 1 million years

# rewritten 2017-06-14 after talking with Mark & Curtis.  I rock.
prob_sampling_paraclade_unknown_size_per_time_step <- function(p,q,r,bs,prec,per_ma=T) {
# p: interval origination rate
# q: interval extinction rate
# r: interval sampling rate
# bs: span of each interval (bs[i] corresponding to p[i],q[i],r[i])
# routine to get phis at different points in time
# returns a vector of length = length of 10xtimescale duration
bins <- length(p);
#slices <- 1+(max(ts)-min(ts))/prec
slices <- 1+(sum(bs)/prec);
if (per_ma)	{
	pp <- assign_interval_rates_to_finer_time_scale(bin_rates=p,bin_spans=bs,prec)
	qq <- assign_interval_rates_to_finer_time_scale(bin_rates=q,bin_spans=bs,prec)
	rr <- assign_interval_rates_to_finer_time_scale(bin_rates=r,bin_spans=bs,prec)
	}	else	{
	# 2018-01-05: modified this so that the bin rate is a per time unit
	pp <- assign_interval_rates_to_finer_time_scale(bin_rates=p/bs,bin_spans=bs,prec)
	qq <- assign_interval_rates_to_finer_time_scale(bin_rates=q/bs,bin_spans=bs,prec)
	rr <- assign_interval_rates_to_finer_time_scale(bin_rates=r/bs,bin_spans=bs,prec)
	}
#dtt <- create_vector_of_single_value(prec,slices)
#dt <- rep(prec,slices)
#### 2018-01-05: I think that this needs to be modified to allow for bs!=1
dt <- rep(prec,length(pp))
#bbss <- assign_interval_rates_to_finer_time_scale(bin_rates=bs,bin_spans=bs,prec)
#dt <- dt/bbss
# debug this!!!!!
#pp[1:52]
#length(pp)
#length(qq)
#length(rr)
#length(dt)
phi_sl <- accio_p_sample_lineage_present_over_interval_continuous(pp,qq,rr,dt)
#sum(abs(phi_sl-foote_99))
#sum(abs(phi_sl_f-foote_99_ts))

return(phi_sl)
}

# rewritten 2018-01-12 to do this more simply
prob_sampling_clade_unknown_size_per_time_step <- function(p,q,r,bs,prec,per_ma=TRUE) {
# p: interval origination rate
# q: interval extinction rate
# r: interval sampling rate
# bs: span of each interval (bs[i] corresponding to p[i],q[i],r[i])
# routine to get phis at different points in time
# returns a vector of length = length of 10xtimescale duration
bins <- length(p)
#slices <- 1+(max(ts)-min(ts))/prec
slices <- 1+(sum(bs)/prec)
if (per_ma)	{
	pp <- assign_interval_rates_to_finer_time_scale(bin_rates=p,bin_spans=bs,prec)
	qq <- assign_interval_rates_to_finer_time_scale(bin_rates=q,bin_spans=bs,prec)
	rr <- assign_interval_rates_to_finer_time_scale(bin_rates=r,bin_spans=bs,prec)
	}	else	{
	# 2018-01-05: modified this so that the bin rate is a per time unit
	pp <- assign_interval_rates_to_finer_time_scale(bin_rates=p/bs,bin_spans=bs,prec)
	qq <- assign_interval_rates_to_finer_time_scale(bin_rates=q/bs,bin_spans=bs,prec)
	rr <- assign_interval_rates_to_finer_time_scale(bin_rates=r/bs,bin_spans=bs,prec)
	}

#### 2018-01-05: I think that this needs to be modified to allow for bs!=1
mxS <- 1
tstf <- 1
mxf <- 0
ttt <- median(bs)
#expected_births <- bs*p
test_bin <- match(max(bs*p),bs*p)
while (tstf > mxf/(10^10))	{
	tstf <- prob_n_species_at_time_tm(n=mxS,a=1,p=pp[test_bin],q=qq[test_bin],chrons=bs[test_bin])
	if (mxf < tstf)	mxf <- tstf
	mxS <- mxS+1
	}

phi_sl <- survivorship <- c()
tslices <- length(pp)
for (sl in 1:tslices)	{
	homer <- fitness <- c()
	for (i in 0:mxS)	{
	#	integrate(prob_n_species_at_time_tm,lower=0,upper=chrons,n=0,a=1,p=p,q=q)$value/chrons
		homer <- c(homer,integrate(prob_n_species_at_time_tm,lower=0,upper=prec,n=i,a=1,p=pp[sl],q=qq[sl])$value)
		fitness <- c(fitness,prob_n_species_at_time_tm(n=i,a=1,p=pp[sl],q=qq[sl],chrons=prec))
		}
	S_homer <- (0:mxS)*homer
	phi_sl <- c(phi_sl,sum(1-dpois(0,S_homer*rr[sl])))
	survivorship <- rbind(survivorship,fitness)
	}

phi <- phi_sl
for (sl in (tslices-1):1)	{
	phi[sl] <- phi[sl]+((1-phi[sl])*sum(survivorship[sl,]*(1-(1-phi[sl+1])^(0:mxS))))
	}
return(phi)
}

# written 2017-06-14 after talking with Mark & Curtis.  I probably have this somewhere already
prob_sampled_branching_per_time_step <- function(p,q,r,bs,prec)	{
phi_sl <- prob_sampling_paraclade_unknown_size_per_time_step(p,q,r,bs,prec);
pp <- assign_interval_rates_to_finer_time_scale(bin_rates=p,bs,prec);
ppp <- 1-dpois(0,prec*pp)	# p branching over time=prec give expectation of p over time=1
return(phi_sl*ppp)	# probability of branching x probability that branching is sampled
}

calculate_paraclade_sampling_per_time_step <- function(pp,qq,rr,dt,ts,prec) {
# routine to get phis at different points in time
# returns a vector of length = length of 10xtimescale duration
# ts: time scale giving onset of each interval.  If there are 12 intervals, then there should
#	be 13 intervals.
bins <- length(pp)
slices <- 1+(max(ts)-min(ts))/prec
ppp <- qqq <- rrr <- dtt <- vector(length=slices)
b <- 1
chrons <- ts[b]
for (i in 1:slices)	{
	dtt[i] <- prec
	if (chrons>=ts[b+1])	b <- b+1
	ppp[i] <- pp[b]
	qqq[i] <- qq[b]
	rrr[i] <- rr[b]
	chrons <- chrons+prec
	}
phi1s <- accio_p_sample_lineage_present_over_interval_continuous(ppp,qqq,rrr,dtt)
prob_fs <- accio_p_sample_lineage_present_over_interval_continuous(pp,qq,rr,dt)
b <- 1
onset <- chrons <- ts[b]
while (chrons<ts[bins+1])	{
	mdpt <- (dt[b]-abs(chrons-onset))/dt[b]
	pmd <- mid_bin_probs_subsequent_sampling_continuous(stg=b,mdpt=mdpt,pp,qq,rr,dt,prob_fs)
	if (chrons==ts[1])	{
		phi <- pmd
		}	else phi <- c(phi,pmd)
	chrons <- chrons+prec
	if (chrons>=ts[b+1])	{
		b <- b+1
		onset <- ts[b]
		}
	}
return(phi)
}

probability_ind_gaps_varying_rates <- function(pp,qq,rr,ts,phis,prec)	{
# was likelihood_ind_gaps_varying_rates for some reason
# determines the likelihood of 0 finds on phylogeny over particular gaps
pgaps <- vector(length=1)
b <- bins <- length(pp)
bs <- bb <- abs(ts[bins+1]-ts[1])/prec
rpl <- rr[b]+(phis[bb]*pp[b])
pgaps <- c(exp(-rpl*prec),1.0)
while (bb>1)	{
	bb <- bb-1
	myr <- ts[bins+1]-(prec*(bs-bb))
#	print(c(b,bb,round(myr,1),round(ts[b],1)))
	if (round(abs(myr),1)==round(abs(ts[b]),1))	b <- b-1
	rpl <- rr[b]+(phis[bb]*pp[b])
	pgaps <- c(exp(-rpl*prec),pgaps)
	}
return(pgaps)
}

prob_sampling_taxon <- function(q,r)	{
# 2016-02-29: simple probability of sampling a taxon before it becomes extinct
# NOTE: Use this only if q & r are static....
return(r/(r+q))
}

individual_lineage_myr <- function(q,chrons)	{
# get the expected lineage myrs for a single lineage taking into account probability of extinction
return(exp(-q*chrons))
}

prob_ever_sampling_lineage_per_interval <- function(qq,rr,ts)	{
bins <- length(pp)
p_samp_ever <- vector(length=bins)
psurv <- 1
for (b in 1:bins)	{
	chrons <- abs(ts[b+1]-ts[b])
	lmyr <- lmyr <- integrate(individual_lineage_myr,lower=0,upper=chrons,q=qq[b])$value
	p_samp_ever[b] <- exp(-rr[b]*lmyr)
	psmp <- 1-exp(-rr[b]*lmyr)
	psurv <- vector(length=b)
	psurv[b] <- 1
	# go through remaining bins.  Get P[find lineage] x P[lineage extant]
	if (b<bins)	{
		for (sb in (b+1):bins)	{
			psv <- psurv[sb-1]*exp(-qq[sb-1]*chrons)
			psurv <- c(psurv,psv)
			chrons <- abs(ts[sb+1]-ts[sb])
			lmyr <- lmyr <- integrate(individual_lineage_myr,lower=0,upper=chrons,q=qq[sb])$value
#			p_samp_ever[b] <- p_samp_ever[b]+psv*exp(-rr[sb]*lmyr)
			psmp <- c(psmp,psv*(1-exp(-rr[sb]*lmyr)))
			}  # HERE!
		}
	pmiss <- 1-psmp
	p_samp_ever[b] <- 1-prod(pmiss)
	}
names(p_samp_ever) <- names(pp)
return(p_samp_ever)
}

prob_sampling_taxon_variable_rates <- function(start,qq,rr,ts)	{
# 2016-03-27
# replaces r/(r+q) when r and/or q vary over time
bins <- length(qq)
sb <- 1
while(round(abs(start),2)<=round(abs(ts[sb+1]),2) && sb<(bins+1))	sb <- sb+1	# get starting inteval
pfind <- 0
if (sb<=bins)	{
	for (b in sb:bins)	{
			# get probability that lineage made it to this stage
		dur <- abs(abs(ts[b+1])-min(abs(start),abs(ts[b])))
		psv <- exp(-qq[b]*dur)	# probability of surviving
		lmyr <- integrate(individual_lineage_myr,lower=0,upper=dur,q=qq[b])$value
		psm <- 1-exp(-rr[b]*lmyr)# probability of being sampled given sampling & chance of extinction
		if (b==sb)	{
			psurv <- c(1,psv)
			psamp <- psm
			} else {
			psurv <- c(psurv,psv*psurv[1+b-sb])
			psamp <- c(psamp,psm)
			}
		}
	psamp <- c(psamp,0)
	if (length(psamp)==length(psurv))	{
		pssm <- psamp*psurv
		}	else {
		pssm <- vector(length=length(psamp))
		for (g in 1:length(psamp))	pssm[g] <- psamp[g]*psurv[g]
		}
	pmiss <- 1-pssm
	pfind <- 1-prod(pmiss)
	}
return(pfind)
}

likelihood_ghost_lineage_duration_varying_rates <- function(start,pgaps,ts,prec)	{
b <- 1
d <- ceiling(abs(log10(prec)))
while(round(abs(start),d)<round(abs(ts[b+1]),d))	b <- b+1	# get starting inteval
bs <- bb <- (start-ts[1])/prec
lgl <- vector(length=1)
lgl[1] <- pgaps[bb]

while (bb>1)	{
	bb <- bb-1
	lgl <- c(lgl,prod(pgaps[bb:bs]))
	}
return(lgl)
}

likelihood_ghost_taxon_duration_varying_rates <- function(start,pp,qq,rr,ts,pgaps,phis,prec)	{
b <- 1
d <- ceiling(abs(log10(prec)))
# if outset of interval, then we need rates from prior interval
while(round(abs(start),d)<round(abs(ts[b+1]),d))	b <- b+1	# get starting inteval
sb <- b
if (start<max(ts))	{
	psap <- prob_sampling_taxon_variable_rates(start,qq,rr,ts)
	} else psap <- 0.0
bs <- bb <- (start-ts[1])/prec
psab <- exp(-pp[b]*prec)						# P[still direct ancestor at td]
punsmp <- 1-(psab*psap)						# P[sampled not filling gap with ancestor]
lgl <- lgt <- vector(length=1)
pg <- pgaps[bb]
lgt[1] <- punsmp*pg
while (bb>1)	{
	bb <- bb-1
	tach <- bs-bb
	myr <- start-(tach*prec)
	if (round(abs(myr),1)==round(abs(ts[b]),1))	b <- b-1
	psab <- psab*exp(-pp[b]*prec)		# keep this at prec
	punsmp <- 1-(psab*psap)						# P[sampled not filling gap with ancestor]
#	rpl <- rr[b]+(phis[bb]*pp[b])
#	lgl <- c(lgl,lgl[tach-1]*exp(-rpl*prec))
	pg <- prod(pgaps[bb:bs])
#	pgt <- punsmp*lgl[tach]
	lgt <- c(lgt,punsmp*pg)
	}
return(lgt)
}

quantiles_ghost_lineage_duration_varying_rates <- function(start,pp,qq,rr,ts,pgaps,phis,quants,prec)	{
#lgl <- likelihood_ghost_lineage_duration_varying_rates(start,pp,qq,rr,ts,phis,prec)
lgl <- likelihood_ghost_lineage_duration_varying_rates(start,pgaps,ts,prec)
quants <- sort(quants,decreasing=FALSE)
qu <- length(quants)
ret_quants <- vector(length=qu)
pgl <- lgl/sum(lgl)
for (i in 2:length(pgl))	pgl[i] <- pgl[i-1]+pgl[i]
j <- i <- 1
while (j<=qu)	{
	if (pgl[i]>=quants[j])	{
		ret_quants[j] <- i*prec
		j <- j+1
		}
	i <- i+1
	}
return(ret_quants)
}

quantiles_ghost_taxon_duration_varying_rates <- function(start,pp,qq,rr,ts,pgaps,phis,quants,prec)	{
lgt <- likelihood_ghost_taxon_duration_varying_rates(start,pp,qq,rr,ts,pgaps,phis,prec)
pgt <- lgt/sum(lgt)
subbins <- length(pgt)
for (bb in 2:subbins)	pgt[bb] <- pgt[bb]+pgt[bb-1]
quants <- sort(quants,decreasing=FALSE)
qu <- length(quants)
j <- i <- 1
ret_quants <- vector(length=qu)
while (j<=qu)	{
	if (pgt[i]>=quants[j])	{
		ret_quants[j] <- i*prec
		j <- j+1
		}
	i <- i+1
	}
return(ret_quants)
}

mean_ghost_lineage_duration <- function(p,q,r)	{
pfind <- prob_sampling_clade_bapst(p,q,r)	# P clades is found
egl <- 1/(r+(pfind*p))		# P no finds and no sampled sister taxa
#pgl <- dpois(1,(pfind*p)*2*chrons)*dpois(0,r*chrons)
return(egl)
}

mean_ghost_taxon_duration <- function(p,q,r)	{
pgt <- vector(length=10000)
pgta <- vector(length=10000)
for (i in 1:10000)	{
	it <- i/100
	pgt[i] <- prob_ghost_taxon_duration(p,q,r,it)
	if (i==1)	{
		pgta[i] <- pgt[i]
		}	else {
		pgta[i] <- pgta[i-1]+pgt[i]
		}
	}
ar <- pgta[10000]
pgta <- pgta/ar
pgt <- pgt/ar

pgexp <- 0
for (i in 1:10000)	pgexp <- pgexp+pgt[i]*i/100
#pglexp <- 0
#for (i in 1:10000)	pglexp <- pglexp+pgl[i]*i/100
#1/(f)
return(pgexp)
}

median_ghost_lineage_duration <- function(p,q,r)	{
pfind <- prob_sampling_clade_bapst(p,q,r)	# P clades is found
mgl <- log(2)/(r+(pfind*p))		# P no finds and no sampled sister taxa
#pgl <- dpois(1,(pfind*p)*2*chrons)*dpois(0,r*chrons)
return(mgl)
}

median_ghost_lineage_duration_alt <- function(p,q,r)	{
pfind <- prob_sampling_clade_bapst(p,q,r)	# P clades is found
chrons <- it <- 0.01
pgla <- vector(length=1)
pgltst <- mxpgl <- pgla[1] <- pgl <- exp(-(r+(pfind*p))*chrons)		# P no finds and no sampled sister taxa
a <- 1
while (pgl>(mxpgl/10000))	{
	chrons <- chrons+it
	pgl <- exp(-(r+(pfind*p))*chrons)
	pgltst <- c(pgltst,pgl)
	pgla <- c(pgla,pgla[a]+pgl)
	if (mxpgl<pgl)	mxpgl <- pgl
	a <- a+1
	}
ar <- pgla[a]
pgla <- pgla/ar
i <- 1
while (pgla[i]<0.5)	i <- i+1
mgl <- it*((i-1)+((0.5-pgla[i-1])/((pgla[i]-0.5)+(0.5-pgla[i-1]))))
	
return(mgl)
}

median_ghost_taxon_duration <- function(p,q,r)	{
chrons <- it <- 0.01
pgta <- vector(length=1)
pgttst <- mxpgt <- pgta[1] <- pgt <- prob_ghost_taxon_duration(p,q,r,chrons)
a <- 1
while(pgt>=(mxpgt/100000))	{
	chrons <- chrons+it
	pgt <- prob_ghost_taxon_duration(p,q,r,chrons)
	pgttst <- c(pgttst,pgt)
	pgta <- c(pgta,pgta[a]+pgt)
	if (mxpgt<pgt)	mxpgt <- pgt
	a <- a+1
	}

ar <- pgta[a]
pgta <- pgta/ar
i <- 1
while (pgta[i]<0.5)	i <- i+1
mgt <- it*((i-1)+((0.5-pgta[i-1])/((pgta[i]-0.5)+(0.5-pgta[i-1]))))
	
return(mgt)
}

median_ghost_taxon_duration_old <- function(p,q,r)	{
pgla <- vector(length=10000)
for (i in 1:10000)	{
	it <- i/100
	pgt <- prob_ghost_taxon_duration(p,q,r,it)
	if (i==1)	{
		pgla[i] <- pgt
		}	else {
		pgla[i] <- pgla[i-1]+pgt
		}
	}
ar <- pgla[10000]
pgla <- pgla/ar
i <- 1
while (pgla[i]<0.5)	i <- i+1
mgt <- ((i-1)+((0.5-pgla[i-1])/((pgla[i]-0.5)+(0.5-pgla[i-1]))))/100

return(mgt)
}

prob_miss_ancestor_before <- function(p,r,chrons)	{
pma <- exp(-(r+p)*chrons)		# P no finds and no sampled sister taxa
#pgl <- dpois(1,(pfind*p)*2*chrons)*dpois(0,r*chrons)
return(pma)
}

expected_ghost_lineage_duration <- function(p,q,r)	{
pfind <- prob_sampling_clade_bapst(p,q,r)	# P clades is found
mgl <- log(2)/(r+(pfind*p))
return(mgl)
}

expected_ghost_taxon_duration <- function(p,q,r)	{
chrons <- inc <- 0.1
ppg <- go <- 1
mxpg <- 0
while (go==1)	{
	ppg <- prob_ghost_taxon_duration(p,q,r,chrons)
	if (mxpg==0)	{
		track <- c(chrons,ppg)
		}	else {
		track <- rbind(track,c(chrons,ppg))
		}
	if (ppg>mxpg)	{
		mxpg <- ppg
		} else if (ppg<(mxpg/1000000))	{
		go <- 0
		}
#	print(c(chrons,ppg,mxpg))
	chrons <- chrons+inc
	}
alpha <- track[,2]/sum(track[,2])
for (i in 2:length(alpha))	{
	alpha[i] <- alpha[i]+alpha[i-1]
	if (alpha[i]>=0.5 && alpha[i-1]<0.5)	dd <- i
	}
xpct <- inc*(dd+(1-(alpha[dd]-0.5)/(alpha[dd]-alpha[dd-1])))
return(xpct)
}

# estimate probability of paraclade going extinct by time chrons
prob_paraclade_extinction <- function(a=1,p,q,chrons)	{
# from Raup 1985
if (p!=q)	{
	ps <- q*(exp((p-q)*chrons)-1)/((p*exp((p-q)*chrons))-q)	# Raup 1985 Eq. A13
	}	else
	ps <- (q*chrons)/(1+(q*chrons))							# Raup 1985 Eq. A12
return(ps^a)
}

# estimate probability of paraclade being still extant at time chrons
prob_paraclade_survival <- function(a=1,p,q,chrons)	{
# from Raup 1985
if (p!=q)	{
	ps <- q*(exp((p-q)*chrons)-1)/((p*exp((p-q)*chrons))-q)	# Raup 1985 Eq. A13
	}	else
	ps <- (q*chrons)/(1+(q*chrons))							# Raup 1985 Eq. A12
return(1-(ps^a))
}

# written 2017-01-05
# modified from Feller's 1966 equations
# get probability of S=n species starting from S=a species
prob_n_species_at_time_tm_pure_birth <- function(n,a=1,p,chrons)	{
if (n==0)	{
	return(0)
	}	else if (n==1)	{
	return(dpois(0,p*chrons))
	}	else	{
	return(choose(n-a,n-1)*exp(-p*a*chrons)*((1-exp(-p*chrons))^(n-a)))
	}
}

# get probability of S=n species starting from S=a species
prob_n_species_at_time_tm <- function(n,a,p,q,chrons)  {
# probability that a paraclade starting with a species has n species
#	after time chrons; from Raup 1985
#print(c(n,a))
if (p==0 && q==0)	{
	pn <- 1.0;
	}	else if (n>0 && q==0)	{
	pn <- prob_n_species_at_time_tm_pure_birth(n,a,p,chrons);
	} else if (n==0)	{
	pn <- prob_paraclade_extinction(a,p,q,chrons)
	} else {
	if (a==1)	{
		if (p!=q)	{
			# Raup 1985 Eq. A13
			alpha <- prob_paraclade_extinction(1,p,q,chrons)
			beta <- alpha*(p/q)
			# Raup 1985 Eq. A17
			pn <- (1-alpha)*(1-beta)*beta^(n-1)
			}	else {
			# Raup 1985 Eq. A15
			# prob of n species given 1 initial species & p=q
			pn <- ((p*chrons)^(n-1))/((1+(p*chrons))^(n+1));
			}	# end case where turnover rates are unequal and we start with 1.  This should be the most commonly used
		}	else {
		pt2 <- pt <- 0
		pts_p <- pts <- pt2s <- c();
		if (p!=q)	{
			# Raup 1985 Eq. A13
			alpha <- prob_paraclade_extinction(1,p,q,chrons);
			beta <- alpha*(p/q);
			# Raup 1985 Eq. A16 pt 2
			for (j in 0:min(a,n))	{
				pt2s <- c(pt2s,(choose_with_logs(a,j)*choose_with_logs((a+n-j-1),(a-1))*(alpha^(a-j))*(beta^(n-j))*((1-alpha-beta)^j)));
				pt2 <- pt2+    (choose_with_logs(a,j)*choose_with_logs((a+n-j-1),(a-1))*(alpha^(a-j))*(beta^(n-j))*((1-alpha-beta)^j));
				pt3 <- log(choose_with_logs(a,j))+log(choose_with_logs(n-1,j-1))+(a-j)*log(alpha)+(n-j)*log(beta)+j*log((1-alpha)*(1-beta));
				pts <- c(pts,(choose(a,j)*choose((a+n-j-1),(a-1))*(alpha^(a-j))*(beta^(n-j))*((1-alpha-beta)^j)));
				pt <- pt+    (choose(a,j)*choose((a+n-j-1),(a-1))*(alpha^(a-j))*(beta^(n-j))*((1-alpha-beta)^j));
				pts_p <- c(pts_p,pt);
				pts_p3 <- c(pts_p3,pt3);
				}
			pn <- pt
			}	else {
			# case with equal extinction & origination rates
			# NOTE: There are errors in Raup's 1985 equations. I have Stanley et al.'s 1981 Equation 3 instead.
			if (n==0)	{
				pn <- ((p*chrons)/(1+(p*chrons)))^a		# prob of extinction
				}	else	{
				pn <- 0
				for (j in 1:a)
					pn <- pn+choose(a,j)*choose(n-1,j-1)*((p*chrons)^(a+n-(2*j)))
				pn <- pn*((1+(p*chrons))^-(a+n))
				}
			}	# end case where we start with 2+ taxa and turnover rates are equal
		}	# end case where we start with 2+ taxa
	}	# end case of non-extinction
return(pn)
}

# written 2017-01-05
# given probabilities of standing richness, what is the probability of
# 	z standing survivors given extinction rate?
conditional_prob_survivors_given_pulsed_extinction <- function(Z,fq)	{
# Z: vector giving probability of 1…z species existing
# fq: probability of a lineage surviving.
#	NOTE: if you start with Poisson expectations, then convert it to
	#	probability using Poisson_rate_to_probability above.
mx_z <- length(Z)
survivor_matrix <- matrix(0,mx_z,(mx_z+1))
for (z in 1:mx_z)	{
	survivor_matrix[z,(1:(z+1))] <- Z[z]*sapply((z:0),dbinom,z,fq)
	}
return(colSums(survivor_matrix))
}

prob_richnesses_at_time_tm <- function(a,p,q,chrons,cutoff=100000000)  {
# prob of going from a to 0…n species at the end of time chrons given origination p & extinction q
if (p<=0 && q<=0)	{
	prich <- c(0,1)
	} else	{
	n <- 0
	prich <- pn <- mxpn <- prob_paraclade_extinction(a,p,q,chrons)
	mxn <- mxprich <- 0
	while (pn>=(mxpn/cutoff))	{
		n <- n+1
		pn <- prob_n_species_at_time_tm(n,a,p,q,chrons)
		prich <- c(prich,pn)
		if (mxpn<pn)	mxpn <- pn
		}
	}
return(prich)
}

prob_final_S_multiple_bins <- function(pp,qq,dt)	{
# 2016-03-22
# create matrix giving P[1…N species | p,q,chrons,a=1]
# note: row 1 is extinction
bins <- length(pp)
mxs <- 1
p_final_S <- matrix(0,mxs,bins)
for (b in 1:bins)	{
	pfs <- prob_richnesses_at_time_tm(1,pp[b],qq[b],dt[b])
	s <- length(pfs)
	if (s>mxs)	{
		dummy <- matrix(0,(s-mxs),bins)
		p_final_S <- rbind(p_final_S,dummy)
		mxs <- s
		}	else if (s<mxs)	{
		dummy <- vector(length=(mxs-s))
		pfs <- c(pfs,dummy)
		}
	p_final_S[,b] <- pfs
	}
rich <- vector(length=mxs)
for (i in 1:mxs)	rich[i] <- i-1
rownames(p_final_S) <- rich
return(p_final_S)
}

prob_n_at_ti_extinct_at_end <- function(n,p,q,ti,chrons)	{
# routine to get the probability of n species at point ti given extinction by point chrons
pn <- prob_n_species_at_time_tm(n,1,p,q,ti)
pex <- prob_paraclade_extinction(n,p,q,(chrons-ti))
return(pn*pex)
}

prob_n_at_ti_extant_at_end <- function(n,p,q,ti,chrons)	{
# routine to get the probability of n species at point ti given extant at point chrons
pn <- prob_n_species_at_time_tm(n,1,p,q,ti)
psr <- 1-prob_paraclade_extinction(n,p,q,(chrons-ti))
return(pn*psr)
}

prob_n_at_ti_extant_at_end_diversifies <- function(n,p,q,ti,chrons)	{
# routine to get the probability of n species at point ti given speciation and net survival
# NOTE: this is scaled to P[survive]-P[stasis]!!!
pn <- prob_n_species_at_time_tm(n,1,p,q,ti)
pstasis <- 0
if (n==1)	pstasis <- exp(-(p+q)*chrons)	# probability that one lineage survives entire interval with no daughters
pdivrs <- 1-(prob_paraclade_extinction(n,p,q,(chrons-ti)))
return((pn*pdivrs)-pstasis)
}

integrate_lineage_time_given_extinction <- function(p,q,chrons,cutoff=100000000)	{
n <- 1
#integrate over P[n | p, q, ti] from ti=0…chrons given n==0 at chrons
# routine integrate over density function for P[n given extinction by chrons]
# routine continues until the area under DF's is trivially small
# ev time = n x DF: if n = 2, then there are two such curves
ev_time <- vector(length=n)
ti <- chrons
if (q>0)	{
	pdodo <- prob_paraclade_extinction(1,p,q,chrons)
	mxt <- ev_time[n] <- (integrate(prob_n_at_ti_extinct_at_end,lower=0,upper=ti,n=n,p=p,q=q,chrons=chrons)$value)/pdodo
	while (ev_time[n]>(mxt/cutoff))	{
		n <- n+1
		ar <- (integrate(prob_n_at_ti_extinct_at_end,lower=0,upper=ti,n=n,p=p,q=q,chrons=chrons)$value)/pdodo
		ev_time <- c(ev_time,n*ar)
		if (mxt<ev_time[n])	mxt <- ev_time[n]
		}
	return(sum(ev_time))
	} else	{
	return(0)
	}
}

integrate_lineage_time_given_survival <- function(p,q,chrons,cutoff=100000000)	{
#integrate over P[n | p, q, ti] from ti=0…chrons given n≥1 at chrons
n <- 1
ev_time <- vector(length=n)
ti <- chrons
proach <- 1-prob_paraclade_extinction(1,p,q,chrons)
mxt <- ev_time[n] <- (integrate(prob_n_at_ti_extant_at_end,lower=0,upper=ti,n=n,p=p,q=q,chrons=chrons)$value)/proach
while (ev_time[n]>(mxt/cutoff))	{
	n <- n+1
	ar <- (integrate(prob_n_at_ti_extant_at_end,lower=0,upper=chrons,n=n,p=p,q=q,chrons=chrons)$value)/proach
	ev_time <- c(ev_time,n*ar)
	if (mxt<ev_time[n])	mxt <- ev_time[n]
	}
return(sum(ev_time))
}

integrate_lineage_time_given_diversification <- function(p,q,chrons,cutoff=100000000)	{
#integrate over P[n | p, q, ti] from ti=0…chrons given n≥1 at chrons
n <- 1	# number of lineages
ev_time <- vector(length=n)
ti <- chrons
pdull <- prob_single_lineage_throughout(p,q,chrons)
psucc <- 1-(prob_paraclade_extinction(1,p,q,chrons)+pdull)
mxt <- ev_time[n] <- (integrate(prob_n_at_ti_extant_at_end,lower=0,upper=ti,n=n,p=p,q=q,chrons=chrons)$value-(pdull*chrons))/psucc
if (ev_time[n]==-Inf)	mxt <- ev_time[n] <- 0
while (ev_time[n]>(mxt/cutoff))	{
	n <- n+1
	ar <- (integrate(prob_n_at_ti_extant_at_end,lower=0,upper=chrons,n=n,p=p,q=q,chrons=chrons)$value)/psucc
	ev_time <- c(ev_time,n*ar)
	if (mxt<ev_time[n])	mxt <- ev_time[n]
	}
return(sum(ev_time))
}

prob_single_lineage_throughout <- function(p,q,chrons)	{
psurv <- exp(-q*chrons)
pnof1 <- exp(-p*chrons)
return(psurv*pnof1)
}

prob_sampling_branching_integrate <- function(p,q,r,chrons)	{
# 2016-03-21: estimate probability of sampling a lineage starting at time 0 through time chrons given separate 
#	expectations for net extinction and net survivorship.
exp_mean_dead <- integrate_lineage_time_given_extinction(p,q,chrons)
#exp_mean_surv <- integrate_lineage_time_given_survival(p,q,chrons)
exp_mean_divers <- integrate_lineage_time_given_diversification(p,q,chrons)
exp_prop_stasis <- prob_single_lineage_throughout(p,q,chrons)
exp_prop_dead <- prob_paraclade_extinction(1,p,q,chrons)
exp_prop_divers <- 1-(exp_prop_dead+exp_prop_stasis)
#pmiss <- (exp_prop_dead*exp(-r*exp_mean_dead))+(exp_prop_divers*exp(-r*exp_mean_divers))+(exp_prop_stasis*exp(-r*chrons))
pfind <- (exp_prop_dead*(1-exp(-r*exp_mean_dead)))+(exp_prop_divers*(1-exp(-r*exp_mean_divers)))+(exp_prop_stasis*(1-exp(-r*chrons)))
#pfind <- 1-pmiss
return(pfind)
}

prob_missing_after_bin <- function(pS,pmiss_one)	{
possS <- length(pS)
pfind <- 0
#pn <- vector(length=possS)
for (s in 1:possS)	{
#	pn[s] <- pS[s]*(1-pmiss_one^s)
	pfind <- pfind+(pS[s]*(1-pmiss_one^s))
	}
pmiss_after <- 1-pfind
return(pmiss_after)
}

#pp <- c(1.2,1.0,0.8,1.1)
#qq <- c(0.8,0.8,1.2,0.7)
#rr <- c(0.25,0.15,0.30,0.25)
#dt <- c(1,1,1,1)
accio_p_sample_lineage_present_over_interval_continuous <- function(pp,qq,rr,dt)	{
# function to get the probability that a lineage present at the outset of B bins is ever sampled
# dt: length of bins, not absolute dates
bins <- length(dt);
p_sample_from_bin <- p_miss_in_bin <- vector(length=bins);
#sapply((1:bins),prob_sampling_branching_integrate,pp,qq,rr,dt)
# first estimate the probability of missing in each time slice
rr_tiles <- length(rr)/length(pp);
for (b in 1:bins)	{
	if (rr_tiles==1)	{
		# if using single rates, then get average pfind
		if (pp[b]<=0 && qq[b]<=0)	{
			pfind <- 1-dpois(0,rr[b]*dt[b])
			} else	{
			pfind <- prob_sampling_branching_integrate(p=pp[b],q=qq[b],r=rr[b],chrons=dt[b])
			}
		}	else {
		# if using distributed sampling rates, then get average pfind
		pfind <- 0
		for (rt in 1:rr_tiles)	{
			if (pp[b]<=0 && qq[b]<=0)	{
				pfind <- pfind + 1-dpois(0,rr[b,rt]*dt[b])
				} else	{
				pfind <- pfind + prob_sampling_branching_integrate(pp[b],qq[b],rr[b,rt],dt[b])
				}
			}
		pfind <- pfind/rr_tiles
		}
	p_miss_in_bin[b] <- 1-pfind
 	#	p_sample_in_bin[b] <- prob_sampling_branching_integrate(pp[b],qq[b],rr[b],dt[b])
	}
	
# now estimate the probability of missing survivors of time slice X in slices X+1…N
for (b in bins:1)	{
	# pout gives probability of outcomes given starting at 1 linege
	pout <- prob_richnesses_at_time_tm(1,p=pp[b],q=qq[b],chrons=dt[b],cutoff=10^20);
	possS <- length(pout);
	if (b<bins)	{
		pS <- pout[2:possS]
		# probability of missing one lineage entering the next interval
		pmiss_one <- 1-p_sample_from_bin[b+1]
		# probability of missing one lineage entering the next interval
		pmiss_after <- prob_missing_after_bin(pS,pmiss_one)
		}	else pmiss_after <- 1
	p_sample_from_bin[b] <- 1-(p_miss_in_bin[b]*pmiss_after)
	#	print(b)
	#	b <- b-1
	}
return(p_sample_from_bin)
}

mid_bin_probs_subsequent_sampling_continuous <- function(stg,mdpt,pp,qq,rr,dt,prob_fa)	{
#### NEW IDEA!  2016-03-15
# For b=bin:1
#   get P[0,1,2,3,4,5] surivving from b-1
#	get P[sample | 1] in bin b
#   get ∑P[n surviving] x (1-P[sample])^n
#   This gives the probability of one lineage starting the interval being sampled
#	Keep going downwards
bins <- length(pp)
# get probabilities of finding a lineage entering an interval or any of its descendants in that interval
prob_fd <- prob_sampling_branching_integrate(p=pp[stg],q=qq[stg],r=rr[stg],chrons=mdpt*dt[stg])
#test_fd <- prob_sampling_branching_integrate(p=pp[stg],q=qq[stg],r=rr[stg],chrons=dt[stg])
pout <- prob_richnesses_at_time_tm(1,p=pp[stg],q=qq[stg],chrons=mdpt*dt[stg])
#test_out <- prob_richnesses_at_time_tm(1,p=pp[stg],q=qq[stg],chrons=dt[stg])
pS <- pout[2:length(pout)]
#test_pS <- test_out[2:length(test_out)]
pmiss_one_next <- 1-prob_fa[stg+1]
if (stg<bins)	{
	prob_mp <- prob_missing_after_bin(pS,pmiss_one_next)
#	test_prob_mp <- prob_missing_after_bin(test_pS,pmiss_one_next)
	} else prob_mp <- 1.0
#test_prob_ft <- 1-((1-test_fd)*test_prob_mp)
prob_ft <- 1-((1-prob_fd)*prob_mp)
return(prob_ft)
}

prob_all_outcomes_given_starting_species <- function(n_mn,n_mx,a,p,q,chrons)	{
# n_mn: minimum richness
# n_mx: maximum richness
span <- 1+(n_mx-n_mn)
pr_out <- vector(length=span)
n <- n_mn
for (po in 1:span)	{
	if (n==0)	{
		pn <- prob_paraclade_extinction(a,p,q,chrons)
		} else {
		if (a==1)	{
			if (p!=q)	{
				alpha <- prob_paraclade_extinction(1,p,q,chrons)		# Raup 1985 Eq. A13
#				alpha <- (q*(exp((p-q)*chrons)-1))/(p*exp((p-q)*chrons)-q)		# Raup 1985 Eq. A13
				beta <- alpha*(p/q)
				pn <- (1-alpha)*(1-beta)*beta^(n-1)					# Raup 1985 Eq. A17
				}	else {
				pn <- ((p*chrons)^(n-1))/((1+(p*chrons))^(n+1))					# Raup 1985 Eq. A15
				}	# end case where turnover rates are unequal and we start with 1.  This should be the most commonly used
			}	else {
				pt <- 0
				if (p!=q)	{
					alpha <- prob_paraclade_extinction(1,p,q,chrons)		# Raup 1985 Eq. A13
#					alpha <- (q*(exp((p-q)*chrons)-1))/(p*exp((p-q)*chrons)-q)		# Raup 1985 Eq. A13
					beta <- alpha*(p/q)
					for (j in 0:min(a,n))	{
						pt <- pt+(choose(a,j)*choose((a+n-j-1),(a-1))*(alpha^(a-j))*(beta^(n-j))*((1-alpha-beta)^j))
						}												# Raup 1985 Eq. A16 pt 2
					pn <- pt
					}	else {
					po <- ((p*chrons)/(1+(p*chrons))^(a+n))							# Raup 1985 Eq. A16 pt 1
					for (j in 1:min(a,n))	{
						pt <- pt+(choose(a,j)*choose(n-1,j-1)*((p*chrons)^(-2*j)))
						}												# Raup 1985 Eq. A16 pt 2
					pn <- po*pt
				}	# end case where we start with 2+ taxa and turnover rates are equal
			}	# end case where we start with 2+ taxa
		}	# end case of non-extinction
	pr_out[po] <- pn
	n <- n+1
	}
return(pr_out)
}

expected_evolutionary_time_given_survival <- function(p,q,chrons,prec)  {
# get probability  of n species
postpmat <- prob_richnesses_at_time_tm_given_survival(p,q,chrons,prec)
mxn <- dim(postpmat)[1]
s_probs <- vector(length=mxn)
for (s in 1:mxn)	s_probs[s] <- s*mean(postpmat[s,])
return(chrons*sum(s_probs))
}

expected_evolutionary_time_given_extinction <- function(p,q,chrons,prec)  {
# 2016-03-21: get probability  of n species
# REMEMBER: Proportion of extinct taxa with richness S in any time slice must 
# include proportion with zero already!!!! 
# simple solution: 
#	1) construct matrix of probability of 0…n lineages at time dt
#	2) construct matrix of probability of extinction given 0…n extinctions given
#		remaining time rt & richness
#	3) product of these is P[extinction | S] x P[S | turnover]
postpmat <- prob_richnesses_at_time_tm_given_extinction(p,q,chrons,prec)
mxn <- dim(postpmat)[1]
s_probs <- vector(length=mxn)
for (s in 1:mxn)	s_probs[s] <- s*mean(postpmat[s,])
return(chrons*sum(s_probs))
}

prob_richnesses_at_time_tm_given_survival <- function(p,q,chrons,prec)  {
# 2016-03-21: get probability  of n species
#dt <- rt <- chrons/2
#pext1 <- prob_paraclade_extinction(1,p,q,rt)			# probability of a lineage being dead without issue by the end
mxn <- prich <- mxprich <- 0
while (prich>=(mxprich/1000000))	{
	mxn <- mxn+1
	prich <- prob_n_species_at_time_tm(mxn,1,p,q,chrons)	# prob that we have mxn species
	if (mxprich<prich)	mxprich <- prich
	}
postpmat <- matrix(0,mxn,prec)
#while (pn>=(mxpn/100000))	{
for (i in 1:prec)	{
	dt <- i*chrons/prec
	rt <- chrons-dt
	pext1 <- prob_paraclade_extinction(1,p,q,rt)		# probability of a lineage being dead without issue by the end
	pcladeextn <- prob_paraclade_extinction(1,p,q,dt)		# probability clade is already extinct
	pstillalive <- 1-pcladeextn							# probability clade is still extant
	pcladesurv <- 1-pext1								# probability that clade is going to survive at this point
	prich <- mxprich <- 0
#	while (prich>=(mxprich/1000000))	{
	for (n in 1:mxn)	{
		psurv <- 1-(pext1^n)					# prob that these species or descendants will survive to end
		pn <- prob_n_species_at_time_tm(n,1,p,q,dt)	# prob that we have n species
#		pnc <- pn/pstillalive					# prob that we have n species given clade is alive
		# ((p[n taxa extant] x p[any survive])/p[still alive then])/p[clade going to survive]
		postpmat[n,i] <- prich <- (pn*psurv)	#prop cases of survivors with n species
		if (mxprich<prich)	mxprich <- prich
#		sum(postpmat[,i])
		}
	}
return(postpmat)
}

prob_richnesses_at_time_tm_given_extinction <- function(p,q,chrons,prec)  {
# get probability  of n species
# REMEMBER: Proportion of extinct taxa with richness S in any time slice must 
# include proportion with zero already!!!! 
# simple solution: 
#	1) construct matrix of probability of 0…n lineages at time dt
#	2) construct matrix of probability of extinction given 0…n extinctions given
#		remaining time rt & richness
#	3) product of these is P[extinction | S] x P[S | turnover]
dt <- chrons/2
rt <- chrons/2
pext1 <- prob_paraclade_extinction(1,p,q,rt)			# probability of a lineage being dead without issue by the end
mxn <- prich <- mxprich <- 0
while (prich>=(mxprich/1000000))	{
	mxn <- mxn+1
	prich <- prob_n_species_at_time_tm(mxn,1,p,q,chrons)	# prob that we have mxn species
	if (mxprich<prich)	mxprich <- prich
	}

#while (prich>=(mxprich/1000000))	{
#	mxn <- mxn+1
#	pextn <- pext1^mxn								# prob that these species and descendants will be extinct at end
#	pn <- prob_n_species_at_time_tm(mxn,1,p,q,dt)	# prob that we have mxn species
#	prich <- (pn*pextn) 			#prop cases of survivors with n species
#	if (mxprich<prich)	mxprich <- prich
#	}
#prichmat <- matrix(0,mxn,prec)			# probability of n lineage being dead without issue by the end
#pfextmat <- matrix(0,mxn,prec)			# probability of that clade is extinct by then
postpmat <- matrix(0,mxn,prec)			# probability of having n lineages but final extinction
#while (pn>=(mxpn/100000))	{
for (i in 1:prec)	{
	dt <- (i-1)*chrons/prec									# start at beginning
	rt <- chrons-dt
	pext1 <- prob_paraclade_extinction(1,p,q,rt)			# probability of a lineage being dead without issue by the end
	n <- prich <- mxprich <- 0
	for (n in 1:mxn)	{
		pextn <- pext1^n								# prob that these species and descendants will be extinct at end
#		pfextmat[n,i] <- pextn
		pn <- prob_n_species_at_time_tm(n,1,p,q,dt)	# prob that we have n species
		# ((p[n taxa extant] x p[go extinct])/p[still alive then])
#		prichmat[n,i] <- pn
		postpmat[n,i] <- prich <- (pn*pextn) 			#prop cases of survivors with n species
		if (mxprich<prich)	mxprich <- prich
		}
	}
return(postpmat)
}

matrix_prob_n_species_over_time_t <- function(a,p,q,chrons,prec)	{
mxn <- 0
a <- 1
pn <- mxpn <- prob_paraclade_extinction(a,p,q,chrons)
while (pn>=(mxpn/1000000))	{
	mxn <- mxn+1
	pn <- prob_n_species_at_time_tm(mxn,a,p,q,chrons)
	if (mxpn<pn)	mxpn <- pn
	}
prich <- matrix(0,mxn+1,prec)
for (b in 1:prec)	{
	dt <- chrons*b/prec
	prich[1,b] <- pn <- mxpn <- prob_paraclade_extinction(a,p,q,dt)	# probability of zero taxa
	for (n in 1:mxn)	{
		nn <- n+1
		prich[nn,b] <- pn <- prob_n_species_at_time_tm(n,a,p,q,dt)	# probability of n>0 taxa
		}
	}
return(prich)
}

prob_sampling_descendants_over_time_t <- function(p,q,r,chrons)  {
pn <- mxps <- ps <- 0
n <- 1
ev_time <- vector(length=1)
while (ps>=(mxps/100000))	{
	expected <- vector(length=100)
	for (i in 1:100)	{
		dt <- (chrons*i)/100
		if (p!=q)	{
			alpha <- (q*(exp((p-q)*dt)-1))/(p*(exp((p-q)*dt)-q))	# Raup 1985 Eq. A13
			beta <- alpha*(p/q)
			expected[i] <- (1-alpha)*(1-beta)*beta^(n-1)			# Raup 1985 Eq. A17
			}	else {
			expected[i] <- ((p*dt)^(n-1))/((1+(p*dt))^(n+1))		# Raup 1985 Eq. A15
			}	# end case where turnover rates are unequal and we start with 1.  This should be the most commonly used
		}
	if (n==1)	{
		ev_time[n] <- chrons*mean(expected)
		} else {
		ev_time <- c(ev_time,chrons*n*mean(expected))
		}
	ps <- 1-exp(-r*ev_time[n])
	if (mxps<ps)	mxps <- ps
#	ps
#	n <- n+1
	}
t_time <- sum(ev_time)		# expected total time encompassed by lineage and possible descendants
pn <- 1-exp(-r*t_time)		# 1 minus probabiity of zero finds
return(pn)
}

integrate_prob_final_S_per_bin_continuous <- function(p,q,chrons)	{
# add division of stage duration to prob_n_species routines
# get P[miss]
n <- 1
x <- mx <- pn <- (integrate(prob_n_species_at_time_tm,0,chrons,n=n,a=1,p=p,q=q)$value)/chrons
while (x>(mx/1000000))	{
	n <- n+1
	x <- (integrate(prob_n_species_at_time_tm,0,chrons,n=n,a=1,p=p,q=q)$value)/chrons
	pn <- c(pn,x)
	if (mx<x)	mx <- x
	}
return(pn)
}

integrate_prob_find_per_bin_continuous_cond <- function(p,q,r,chrons)	{
# add division of stage duration to prob_n_species routines
# get P[miss]
n <- 1
pn <- (integrate(prob_n_species_at_time_tm,0,chrons,n=n,a=1,p=p,q=q)$value)/chrons
pf <- (1-exp(-r*chrons*n))		# poisson probability # probability of 1
mxpnf <- pnf <- pn*pf
#dbt <- pn
#dbf <- pf
while ((pn*pf)>(mxpnf/1000000))	{
	n <- n+1
	pn <- (integrate(prob_n_species_at_time_tm,0,chrons,n=n,a=1,p=p,q=q)$value)/chrons
	pf <- (1-exp(-r*chrons*n))
	pnf <- c(pnf,pn*pf)
#	dbt <- c(dbt,pn)	# for debugging
#	dbf <- c(dbf,pf)	# for debugging
	}
return(sum(pnf))
}

integrate_prob_find_per_bin_continuous_sum <- function(p,q,r,chrons)	{
# add division of stage duration to prob_n_species routines
# get P[miss]
n <- 1
ps <- (integrate(prob_n_species_at_time_tm,0,chrons,n=n,a=1,p=p,q=q)$value)/chrons
an <- chrons*ps	# total evolutionary time expected from richness = 1
paf <- vector(length=1)
mxpaf <- paf[1] <- (1-exp(-r*chrons*n*an))	# poisson expected sum of time from 1's
#dbt <- pn
#dbf <- pf
while (paf[n]>(mxpaf/1000000))	{
	n <- n+1
	pn <- (integrate(prob_n_species_at_time_tm,0,chrons,n=n,a=1,p=p,q=q)$value)/chrons
	ps <- c(ps,pn)
	an <- c(an,n*chrons*pn)
	paf <- c(paf,(1-exp(-r*chrons*n*pn)))
	if (mxpaf<paf[n])	mxpaf <- paf[n]
#	dbt <- c(dbt,pn)	# for debugging
#	dbf <- c(dbf,pf)	# for debugging
	}
return((1-exp(-r*sum(an))))
}

#integrate(prob_paraclade_extinction,lower=0,upper=chrons,a=1,p=p,q=q)$value/chrons

per_bin_probs_subsequent_sampling_continuous <- function(pp,qq,rr,dt)	{
#### NEW IDEA!  2016-03-15
# For b=bin:1
#   get P[0,1,2,3,4,5] surivving from b-1
#	get P[sample | 1] in bin b
#   get ∑P[n surviving] x (1-P[sample])^n
#   This gives the probability of one lineage starting the interval being sampled
#	Keep going downwards
bins <- length(pp)
# get probabilities of finding a lineage entering an interval or any of its descendants in that interval
prob_fd <- vector(length=bins)
for (b in 1:bins)	prob_fd[b] <- integrate_prob_find_per_bin_continuous_cond(p=pp[b],q=qq[b],r=rr[b],chrons=dt[b])
# get the probability of finding a lineage entering be during (fd) or after (fa) the bin
prob_f <- prob_fd
for (b in (bins-1):1)	{
	pS <- integrate_prob_final_S_per_bin_continuous(pp[b],qq[b],dt[b])
	psSurv <- length(pS)
	dbpf <- vector(length=psSurv)
	prob_fa <- 0
	for (s in 1:psSurv)	{
		dbpf[s] <- pfind <- 1-(1-prob_f[b+1])^s 	# 1 - prob missing s times
		prob_fa <- prob_fa+(pS[s]*pfind)
		}
	# P finding is 1- (p[miss then] x p[miss after])
	prob_f[b] <- 1-((1-prob_fd[b])*(1-prob_fa))
	}
return(prob_f)
}

per_bin_success_estimation <- function(a,pp,qq,rr,tt)	{
bins <- length(tt)-1

cll <- 1
for (b1 in 1:(bins-1))	{
	pouts <- prob_richnesses_at_time_tm(1,pp[b1],qq[b1],abs(tt[b1]-tt[b1+1]))
	mxsb1 <- length(pouts)
	if (b1==1)	{
		psurv <- matrix(0,mxsb1,1)
		psurv[,1] <- pouts
		}
	if (bins<10)	{
		header1 <- paste("bin_",b,sep="")
		}	else {
		header1 <- paste("bin_0",b,sep="")
		}
	for (b2 in (b1+1):bins)	{
		cll <- cll+1
		if (bins<10)	{
			header2 <- paste("bin_",b2,sep="")
			}	else {
			header2 <- paste("bin_0",b2,sep="")
			}
		if (b==1 && b2==2)	{
			header <- paste(header1,"_to_",header2,sep="")
			} else {
			header <- c(header,paste(header1,"_to_",header2,sep=""))
			}

		for (i in 2:mxsb1)	{	# start at 2 because i=1 gives P[extinct] & i=2 gives P[oS=1]
			a <- i-1
			pouts2 <- prob_richnesses_at_time_tm(a,pp[b2],qq[b2],abs(tt[b2]-tt[b2+1]))
			md <- length(pouts2)
			if (a==1)	{
				pouts2a <- pouts[i]*pouts2		# pouts[i] gives P[oS=a=(i-1)]
				mxsb2 <- md
				}	else {
				if (mxsb2<md)	{
					dummy <- vector(length=(md-mxsb2))
					pouts2a <- c(pouts2a,dummy)
					mxsb2 <- md
					} else if (mxsb2>md) {
					dummy <- vector(length=(mxsb2-md))
					pouts2 <- c(pouts2,dummy)
					}
				pouts2a <- pouts2a+(pouts[i]*pouts2)	# pouts[i] gives P[oS=a=(i-1)]
				}
			}
		# adjust matrix or vector size to accommodate new possibilities
		if (mxsb2>mxsb1)	{
			dummy <- matrix(0,abs(mxsb2-mxsb1),cll-1)
			psurv <- rbind(psurv,dummy)
			mxsb1 <- mxsb2
			}	else if (mxsb2<mxsb1)	{
			dummy <- vector(length=abs(mxsb2-mxsb1))
			pouts2a <- c(pouts2a,dummy)
			}
		psurv <- cbind(psurv,pouts2a)
		pouts <- pouts2a	# prob of 0…n survivors after b2
		}
	}

for (b in 1:bins)	{
	if (b==1)	{
		# note: psucc[1]=p[extinction]
		psucc <- prob_richnesses_at_time_tm(1,pp[b],qq[b],abs(tt[b]-tt[b+1]))
		mxs <- length(psucc)
		}	else {
		xx <- prob_richnesses_at_time_tm(1,pp[b],qq[b],abs(tt[b]-tt[b+1]))
		bns <- length(xx)
		if (bns<mxs)	{
			dummy <- vector(length=(mxs-bns))
			xx <- c(xx,dummy)
			psucc <- cbind(psucc,xx)
			} else if (bns>mxs)	{
			if (b>2)	{
				dummy <- matrix(0,(bns-mxs),(b-1))
				psucc <- rbind(psucc,dummy)
				}	else {
				dummy <- vector(length=(bns-mxs))
				psucc <- c(psucc,dummy)
				}
			psucc <- cbind(psucc,xx)
			mxs <- bns
			}	else {
			psucc <- cbind(psucc,xx)
			}
		}
	}
if (bins<10)	{
	headers <- "bin_1"
	for (i in 2:bins)	headers <- c(headers,(paste("bin_",i,sep="")))
	}	else {
	headers <- "bin_01"
	for (i in 2:9)	headers <- c(headers,(paste("bin_0",i,sep="")))
	for (i in 10:bins)	headers <- c(headers,(paste("bin_",i,sep="")))
	}
progeny <- vector(length=mxs)
for (i in 1:mxs)	progeny[i] <- i-1

colnames(psucc) <- headers
rownames(psucc) <- progeny

return(psucc)
}

prob_sampling_descendants_over_multiple_bins <- function(dt,ps,qs,rs,ts)  {
# dt[i] is the divergence time in question
# ps[i] is origination rate in bin i
# qs[i] is extinction rate in bin i
# rs[i] is sampling rate in bin i
# ts is timescale	
pn <- mxps <- ps <- 0
n <- 1
ev_time <- vector(length=1)
while (ps>=(mxps/100000))	{
	expected <- vector(length=100)
	for (i in 1:100)	{
		dt <- (chrons*i)/100
		if (p!=q)	{
			alpha <- (q*(exp((p-q)*dt)-1))/(p*(exp((p-q)*dt)-q))	# Raup 1985 Eq. A13
			beta <- alpha*(p/q)
			expected[i] <- (1-alpha)*(1-beta)*beta^(n-1)					# Raup 1985 Eq. A17
			}	else {
			expected[i] <- ((p*dt)^(n-1))/((1+(p*dt))^(n+1))					# Raup 1985 Eq. A15
			}	# end case where turnover rates are unequal and we start with 1.  This should be the most commonly used
		}
	if (n==1)	{
		ev_time[n] <- chrons*mean(expected)
		} else {
		ev_time <- c(ev_time,chrons*n*mean(expected))
		}
	ps <- 1-exp(-r*ev_time[n])
	if (mxps<ps)	mxps <- ps
	ps
	n <- n+1
	}
t_time <- sum(ev_time)		# expected total time encompassed by lineage and possible descendants
pn <- 1-exp(-r*t_time)		# 1 minus probabiity of zero finds
return(pn)
}

prob_sampling_descendants_in_interval <- function(p,q,r,chrons)	{
n <- 0
tpd <- pext <- prob_n_species_at_time_tm(n,1,p,q,chrons)
et <- 0.5*chrons
ps <- psm <- mxpsm <- (1-exp(-r*et))*pext

while (psm>=(mxpsm/10000))	{
	n <- n+1
	pds <- prob_n_species_at_time_tm(n,1,p,q,chrons)
	tpd <- tpd+pds
	et <- chrons*mean(c(1,n))
	psm <- (1-exp(-r*et))*pds
	ps <- ps+psm
	if (mxpsm<psm)	mxpsm <- psm
	}
#?do this in increments of 0.01 of chrons?
n <- n+1
pds <- (1-tpd)
et <- chrons*mean(c(1,n))
psm <- (1-exp(-r*et))*pds
ps <- ps+psm

return(ps)
}

expected_evolutionary_time_over_time_t <- function(p,q,chrons)  {
pn <- mxps <- ps <- 0
n <- 1
ev_time <- vector(length=1)
while (ps>=(mxps/100000))	{
	expected <- vector(length=100)
	for (i in 1:100)	{
		dt <- (chrons*i)/100
		if (p!=q)	{
			alpha <- (q*(exp((p-q)*dt)-1))/(p*(exp((p-q)*dt)-q))	# Raup 1985 Eq. A13
			beta <- alpha*(p/q)
			expected[i] <- (1-alpha)*(1-beta)*beta^(n-1)			# Raup 1985 Eq. A17
			}	else {
			expected[i] <- ((p*dt)^(n-1))/((1+(p*dt))^(n+1))		# Raup 1985 Eq. A15
			}	# end case where turnover rates are unequal and we start with 1.  This should be the most commonly used
		}
	if (n==1)	{
		ev_time[n] <- chrons*mean(expected)
		} else {
		ev_time <- c(ev_time,chrons*n*mean(expected))
		}
	ps <- sum(expected)
	if (mxps<ps)	mxps <- ps
	n <- n+1
	sum(ev_time)
	}
t_time <- sum(ev_time)		# expected total time encompassed by lineage and possible descendants
return(t_time)
}

prob_missing_descendants_in_interval <- function(p,q,r,chrons,prec)	{
pmissdt <- vector(length=prec)
debug <- vector(length=20)
for (i in 1:prec)	{
	dt <- i*chrons/prec
	n <- 0
	pds <- mxpds <- 0
	while (pds>=(mxpds/100000))	{
		debug[n+1] <- pds <- prob_n_species_at_time_tm(n,1,p,q,dt)
		if (n==0)	{
			pmissdt[i] <- pmissdt[i]+pds
			}	else {
			pmissdt[i] <- pmissdt[i]+(pds*exp(-r*n*dt))
			}
		if (mxpds<pds)	mxpds <- pds
		n <- n+1
		}
	}
xxx <- vector(length=prec)
for (i in 1:prec)	xxx[i] <- i/prec
plot(xxx,pmissdt,pch=21,cex=0.1)
return(mean(pmissdt))
}

expected_evolutionary_time_given_extinction <- function(p,q,chrons,prec)  {
pext1 <- vector(length=(prec-1))
pclsurv <- vector(length=(prec-1))
for (i in 1:(prec-1))	{
	dt <- i*chrons/prec
	rt <- chrons-dt
	pext1[i] <- prob_paraclade_extinction(1,p,q,rt)
	pclsurv[i] <- 1-prob_paraclade_extinction(1,p,q,dt)
	}
#plot(50*xxx[1:999],pclsurv,ylim=c(0,1))
#points(50*xxx[1:999],pext1)

n <- 1
pn <- mxpn <- 0
debug <- vector(length=(prec-1))
while (pn>=(mxpn/100000))	{
	if (n==1)	{
		pdoomedlines <- vector(length=1)
		}	else {
		pdoomedlines <- c(pdoomedlines,0)
		}
	for (i in 1:(prec-1))	{
		dt <- i*chrons/prec
		rt <- chrons-dt
		pext <- pext1[i]^n
		pds <- prob_n_species_at_time_tm(n,1,p,q,dt)
		debug[i] <- pext*pds
		pdoomedlines[n] <- pdoomedlines[n]+(pext*pds)/(prec-1)
		}
	if (n==1)	{
		plot(50*xxx[1:999],debug,pch=21,bg="red",col="red",cex=0.2)
		}	else if (n==2) {
		points(50*xxx[1:999],debug,pch=21,bg="green",col="green",cex=0.2)
		}	else if (n==3) {
		points(50*xxx[1:999],debug,pch=21,bg="purple2",col="purple2",cex=0.1)
		}	#else if (n==4) {
#		points(50*xxx[1:999],debug,pch=21,bg="purple2",col="purple2",cex=0.1)
#		}
	pn <- pdoomedlines[n]
	if (mxpn<pn)	mxpn <- pn
	n <- n+1
	}
ev_time <- 0
for (s in 1:(n-1))	ev_time <- ev_time+(chrons*s*pdoomedlines[s])
return(ev_time)
}

expected_evolutionary_time_given_extinction_fuzzy <- function(aa,p,q,chrons,prec)  {
# get prob n species given starting at 1, 2, 3, 4, etc. species
# weight by probability of n species given prior intervals
# aa: vector giving probability of a species entering the interval
# get probability  of n species
# REMEMBER: Proportion of extinct taxa with richness S in any time slice must 
# include proportion with zero already!!!! 
# simple solution: 
#	1) construct matrix of probability of 0…n lineages at time dt
#	2) construct matrix of probability of extinction given 0…n extinctions given
#		remaining time rt & richness
#	3) product of these is P[extinction | S] x P[S | turnover]
pn <- mxpn <- 0
mxn <- 1
mxa <- length(aa)	# number of possible introductory lineages that are reasonably probable
	#while (pn>=(mxpn/100000))	{
ev_time <- 0
mxnsf <- 0
for (a in 1:mxa)	{
	debug <- vector(length=prec)
	prichmat <- matrix(0,mxn,prec)
	pfextmat <- matrix(0,mxn,prec)
	postpmat <- matrix(0,mxn,prec)
	for (i in 1:prec)	{
		dt <- (i-1)*chrons/prec									# start at beginning
		rt <- chrons-dt
		pext1 <- prob_paraclade_extinction(1,p,q,rt)			# probability of a lineage being dead without issue by the end
	#	pcladeextn <- prob_paraclade_extinction(1,p,q,dt)		# probability clade is already extinct
	#	pstillalive <- 1-prob_paraclade_extinction(1,p,q,dt)	# probability clade is still extant
	#	pcladesurv <- 1-pext1									# probability that clade is going to survive at this point
		n <- prich <- mxprich <- 0
		while (prich>=(mxprich/100000))	{
			n <- n+1
			if (n>mxn)	{
				added <- vector(length=prec)
				prichmat <- rbind(prichmat,added)
				pfextmat <- rbind(pfextmat,added)
				postpmat <- rbind(postpmat,added)
				mxn <- n
				}
			pextn <- pext1^n								# prob that these species and descendants will be extinct at end
			pfextmat[n,i] <- pextn
			pn <- prob_n_species_at_time_tm(n,1,p,q,dt)	# prob that we have n species
			# ((p[n taxa extant] x p[go extinct])/p[still alive then])
			prichmat[n,i] <- pn
			postpmat[n,i] <- prich <- (pn*pextn) 			#prop cases of survivors with n species
			if (mxprich<prich)	mxprich <- prich
	#		sum(postpmat[,i])
			}
		}
	if (mxnsf<mxn)	{
		if (mxnsf==0)	{
			ss_probs <- vector(length=mxn)
			}	else {
			dummy <- vector(length=(mxn-mxnsf))
			ss_probs <- c(ss_probs,dummy)
			}
		mxnsf <- mxn
		}
	s_probs <- vector(length=mxn)		# keep these separate for debugging purposes
	for (s in 1:mxn)	{
		s_probs[s] <- s*mean(postpmat[s,])
		ss_probs[s] <- ss_probs[s]+(aa[a]*s_probs[s])
		}
	}
# DEBUGGING CRAP: CUT THIS OUT
#palrext <- vector(length=prec)
#xxx <- vector(length=prec)
#for (i in 1:prec)	{
#	palrext[i] <- prob_n_species_at_time_tm(0,1,p,q,(chrons*(i-1)/prec))
#	xxx[i] <- (i-1)/prec
#	}
#plot(xxx,(1-palrext)-min(1-palrext),xlab="Prop of Way through Interval",ylab="P[extinct by end]",pch=21,bg="black",col="black",cex=2)
#lines(xxx,postpmat[1,],bg="red",col="red",lwd=10)
#lines(xxx,postpmat[2,],bg="orange",col="orange",lwd=10)
#lines(xxx,postpmat[3,],bg="cyan",col="cyan",lwd=7.7)
return(chrons*sum(ss_probs))
}

probs_richnesses_in_intervals_varying_rates <- function(start,end,pp,qq,dt)	{
p_S <- prob_richnesses_at_time_tm(1,pp[start],qq[start],dt[start])
tx <- length(p_S)
sttx <- tx-1
p_final_S <- matrix(0,tx,1)
p_final_S[,1] <- p_S
for (b in (start+1):end)	{
	pxn <- mxn <- pnn <- mxpnn <- 0
	nn <- sttx	# 2016-02-29: sttx is the number of lineages; tx is one higher to allow for zero
	# get the maximum number of remotely realistic species
	while (pxn>=(mxpnn/1000000))	{
		pxn <- prob_n_species_at_time_tm(nn,sttx,pp[b],qq[b],dt[b])
		nn <- nn+1
		if (mxpnn<pxn)	mxpnn <- pxn
		}
	poss_outs <- matrix(0,sttx+1,nn)
	for (a in 1:sttx)	{
		xxxx <- prob_all_outcomes_given_starting_species(0,nn-1,a,pp[b],qq[b],dt[b])
		poss_outs[a,] <- xxxx
		}
	
	p_ext <- p_final_S[1,b-1]		# prob. that clade already is extinct
	for (a in 1:(sttx-1))	p_ext <- p_ext+p_final_S[a+1,b-1]*poss_outs[a,1]

	pnn <- vector(length=nn)
	pnn[1] <- p_ext
	for (m in 2:nn)	{
		n <- m-1		# n = St; m is the cell number for St=m in poss_outs to accomodate zero
		for (a in 1:sttx)	{
			aa <- a+1	# aa is the cell number for So=a in p_final_S to accomodate zero
			# prob of starting with So=a x prob of getting St=n given So=a
			pnn[m] <- pnn[m]+(p_final_S[aa,b-1]*poss_outs[a,m])
			}
		}

	if (nn==tx)	{
		p_final_S <- cbind(p_final_S,pnn)
		} else if (nn>tx)	{
		add <- abs(nn-tx)
		dummy <- matrix(0,add,b-1)
		p_final_S <- rbind(p_final_S,dummy)
		p_final_S <- cbind(p_final_S,pnn)
		tx <- nn
		}	else if (nn<tx)	{
		add <- abs(tx-nn)
		dummy <- vector(length=add)
		pnn <- c(pnn,dummy)
		p_final_S <- cbind(p_final_S,pnn)
		}	else {
		p_final_S <- cbind(p_final_S,pnn)
		}
	sttx <- nn
	}

return(p_final_S)
}

expected_sampled_descendants_integrated <- function(p,q,r,prec)	{
#prec: decimal for how finely to divide time units
if (prec>1)	prec <- 1/prec
f <- prob_sampling_clade_bapst(p,q,r)
mq <- q*prec
pnp <- np <- mxpnp <- 0
pnode <- (f*p)/((f*p)+q)	# probability that there are sampled descendants
while (pnp>=(mxpnp/100000))	{
	np <- np+1
	if (np==1)	{
		ppoly <- vector(length=np)
		} else {
		ppoly <- c(ppoly,0)
		}
	pp <- mxpp <- 0
	i <- chrons <- 0
	while (pp>=(mxpp/100000))	{
		chrons <- chrons+prec
		pd <- mq*((1-mq)^i)		# probablity of species living this long
		pp <- dpois(np,p*f*chrons)*pd	# probability of np sampled descendants given p, q, r & chrons
		ppoly[np] <- ppoly[np]+pp
		if (mxpp<pp)	mxpp <- pp
		i <- i+1
		}
	ppoly[np] <- ppoly[np]
	pnp <- ppoly[np]
	if (mxpnp<ppoly[np])	mxpnp <- ppoly[np]
	}
}

expected_sampled_descendants <- function(p,q,r,prec)	{
if (prec>1)	prec <- 1/prec
f <- prob_sampling_clade_bapst(p,q,r)
mxpnp <- 0
np <- 1
pnode <- (f*p)/((f*p)+q)	# probability of sampled descendants
pdesc <- vector(length=np)
pdesc[np] <- (1-pnode)	# probability of one sampled descendant given that there are descendants
while (pdesc[np]>=(mxpnp/100000))	{
	np <- np+1
	pdesc <- c(pdesc,(1-pnode)*pnode^(np-1))	# probability of np sampled descendants given that there are descendants
	if (mxpnp<pdesc[np])	mxpnp <- pdesc[np]
	}
return(pdesc)
}

probability_of_Sfurcation <- function(S,p,q,r)	{
# this is the probability of one sampled ancestor + S-1 sampled descendants or no sampled ancestor and S sampled descendants
f <- prob_sampling_clade_bapst(p,q,r)	# prob the descendant or any deriviative is sampled
#psamp_anc_span <- prob_sampling_taxon(q,r)*prob_sampling_taxon(p,r)	# prob. that ancestor is sampled before and after
psamp_anc_ever <- prob_sampling_taxon(q,r)	# prob of sampling any taxon
psamp_anc_after <- prob_sampling_taxon(q,r)*(1-prob_sampling_taxon(p,r)) # prob of sampling ancestor only after speciation
pnode <- (f*p)/((f*p)+q)		# probability of at least one sampled derivative of a taxon

pStomy <- (1-pnode)*(pnode^(S-1))
return(pStomy)
#if (S>1)	{
#	p_incl_anc <- psamp_anc_ever*(1-pnode)*(pnode^(S-2))
#	p_excl_anc <- (1-psamp_anc_ever)*(1-pnode)*(pnode^(S-1))	# prob of sampling S descendants 
#	}	else {
#	pf1 <- p/(p+q)	# probability of descendants
#	p_incl_anc <- psamp_anc_ever*(pf1-pnode)
#	p_excl_anc <- (1-psamp_anc_ever)
#	}
#return(p_excl_anc+p_incl_anc)
}

furcation_probabilities <- function(p,q,r)	{
psamp_anc <- r/(r+q)
f <- prob_sampling_clade_bapst(p,q,r)	# prob the descendant or any deriviative is sampled
pone <- (f*p)/((f*p)+q)
#p_invis <- ((1-psamp_anc)*pone)+(psamp_anc*(1-pone))
#p_invis <- ((1-psamp_anc)*pone)
p_invis <- psamp_anc*(1-pone)
p_visib <- 1-p_invis

S <- 2
ptomy <- mxp <- 0
pS <- vector(length=2)
pS[S] <- 1-pone
while (ptomy>=(mxp/1000000))	{
	S <- S+1
#	ptomy <- probability_of_Sfurcation(S,p,q,r)
	ptomy <- (1-pone)*(pone^(S-2))
	pS <- c(pS,ptomy)
	if(mxp<pS[S])	mxp <- pS[S]
	}
#xx <- vector(length=S)
#for (i in 1:S)	xx[i] <- i
#plot(xx,log(pS),pch=21,bg="red")
return(pS)
}

prob_n_species_and_no_sampling_at_time_tm <- function(n,a,p,q,r,chrons)  {
# 2016-03-14: uses integration!  Keep this.
if (n==0)	{
	pn <- prob_paraclade_extinction(a,p,q,chrons)
	pn <- pn^a
	} else {
	if (a==1)	{
		if (p!=q)	{
			alpha <- prob_paraclade_extinction(1,p,q,chrons)		# Raup 1985 Eq. A13
#				alpha <- (q*(exp((p-q)*chrons)-1))/(p*exp((p-q)*chrons)-q)		# Raup 1985 Eq. A13
			beta <- alpha*(p/q)
			pn <- (1-alpha)*(1-beta)*beta^(n-1)					# Raup 1985 Eq. A17
			}	else {
			pn <- ((p*chrons)^(n-1))/((1+(p*chrons))^(n+1))					# Raup 1985 Eq. A15
			}	# end case where turnover rates are unequal and we start with 1.  This should be the most commonly used
		}	else {
			pt <- 0
			if (p!=q)	{
				alpha <- prob_paraclade_extinction(1,p,q,chrons)			# Raup 1985 Eq. A13
#				alpha <- (q*(exp((p-q)*chrons)-1))/(p*exp((p-q)*chrons)-q)		# Raup 1985 Eq. A13
				beta <- alpha*(p/q)
				for (j in 0:min(a,n))	{
					pt <- pt+(choose(a,j)*choose((a+n-j-1),(a-1))*(alpha^(a-j))*(beta^(n-j))*((1-alpha-beta)^j))
					}												# Raup 1985 Eq. A16 pt 2
				pn <- pt
			}	else {
				po <- ((p*chrons)/(1+(p*chrons))^(a+n))							# Raup 1985 Eq. A16 pt 1
				for (j in 1:min(a,n))	{
					pt <- pt+(choose(a,j)*choose(n-1,j-1)*((p*chrons)^(-2*j)))
				}												# Raup 1985 Eq. A16 pt 2
				pn <- po*pt
			}	# end case where we start with 2+ taxa and turnover rates are equal
		}	# end case where we start with 2+ taxa
	}	# end case of non-extinction
pmiss <- pn*(exp(-r))^n
return(pmiss)
}

prob_n_species_and_sampling_at_time_tm <- function(n,a,p,q,r,chrons)  {
# 2016-03-14: uses integration!  Keep this.
if (n==0)	{
	pn <- prob_paraclade_extinction(a,p,q,chrons)
	pn <- pn^a
	} else {
	if (a==1)	{
		if (p!=q)	{
			alpha <- prob_paraclade_extinction(1,p,q,chrons)		# Raup 1985 Eq. A13
#				alpha <- (q*(exp((p-q)*chrons)-1))/(p*exp((p-q)*chrons)-q)		# Raup 1985 Eq. A13
			beta <- alpha*(p/q)
			pn <- (1-alpha)*(1-beta)*beta^(n-1)					# Raup 1985 Eq. A17
			}	else {
			pn <- ((p*chrons)^(n-1))/((1+(p*chrons))^(n+1))					# Raup 1985 Eq. A15
			}	# end case where turnover rates are unequal and we start with 1.  This should be the most commonly used
		}	else {
			pt <- 0
			if (p!=q)	{
				alpha <- prob_paraclade_extinction(1,p,q,chrons)			# Raup 1985 Eq. A13
#				alpha <- (q*(exp((p-q)*chrons)-1))/(p*exp((p-q)*chrons)-q)		# Raup 1985 Eq. A13
				beta <- alpha*(p/q)
				for (j in 0:min(a,n))	{
					pt <- pt+(choose(a,j)*choose((a+n-j-1),(a-1))*(alpha^(a-j))*(beta^(n-j))*((1-alpha-beta)^j))
					}												# Raup 1985 Eq. A16 pt 2
				pn <- pt
			}	else {
				po <- ((p*chrons)/(1+(p*chrons))^(a+n))							# Raup 1985 Eq. A16 pt 1
				for (j in 1:min(a,n))	{
					pt <- pt+(choose(a,j)*choose(n-1,j-1)*((p*chrons)^(-2*j)))
				}												# Raup 1985 Eq. A16 pt 2
				pn <- po*pt
			}	# end case where we start with 2+ taxa and turnover rates are equal
		}	# end case where we start with 2+ taxa
	}	# end case of non-extinction
pfind <- pn*(1-((exp(-r))^n))
return(pfind)
}

### WHAT THE FUCK?!?!?
prob_sampling_divergences_per_bin_given_varying_continuous_rates <- function(pp,qq,rr,dt)	{
# t0: divergence time
# pp: vector of origination rates
# qq: vector of extinction rates
# rr: vector of sampling rates
# ts: vector giving onset of different intervals
if (length(tt)==length(pp))	tt <- c(tt,0)
bins <- length(pp)

# get probabilities of n surviving descendants given that an interval starts with 1
pS_end <- vector(length=bins)
for (b in bins:1)	{
	dt <- abs(tt[b+1]-tt[b])
	pS_end[b] <- probs_richnesses_in_intervals_varying_rates(b,b+1,pp[b],qq[b],dt)
	}
}

prob_total_extinction <- function(p,q,chrons)	{
if (p==q)	{
	pext <- (p*chrons)/(1+(p*chrons))
	}	else {
	pext <- q*(exp((p-q)*chrons)-1)/(p*(exp((p-q)*chrons)-q))
	}
return(pext)
}

prob_finding_two_plus_descendants <- function(p,q,r,s)	{
pfind <- prob_sampling_clade_bapst(p,q,r)	# P clades is found
pmiss <- 1-pfind
p0 <- pmiss^s					# P all daughters missed
p1 <- s*pfind*(pmiss^(s-1))		# P all but one daughter missed
p2p <- 1-(p0+p1)				# P 2+ daughters sampled
return(p2p)
}

ancestor_effect <- function(p,q,r,chrons)	{
s <- 2
p2p <- prob_finding_two_plus_descendants(p,q,r,s)
nsm <- ism <- sum <- ((p*chrons)^s)*p2p/factorial(s)
while(nsm>(ism/100000))	{
	s <- s+1
	p2p <- prob_finding_two_plus_descendants(p,q,r,s)
	nsm <- ((p*chrons)^s)*p2p/factorial(s)
	sum <- sum+nsm
	}
return (sum)
}

prob_amalgamate_at_node <- function(p,q,r,chrons)	{
pfind <- prob_sampling_clade_bapst(p,q,r)	# P clades is found
pgam <- 1-exp(-((pfind*p))*chrons)		# P no finds and no sampled sister taxa
#pgl <- dpois(1,(pfind*p)*2*chrons)*dpois(0,r*chrons)
return(pgam)
}

prob_missing_paraclade_continuous <- function(p,q,r,chrons,prec)	{
#written 2016-02-29
pmiss <- vector(length=prec)
pm_dt <- exp(-r*(chrons/prec))
#cps <- prob_richnesses_at_time_tm(1,p,q,chrons)	# conditional probabilities of richnesses
#mxtx <- length(cps)
for (it in 1:prec)	{
	dt <- chrons*it/prec
	cps <- prob_richnesses_at_time_tm(1,p,q,dt)	# conditional probabilities of richnesses
	mxtx <- length(cps)
	for (cs in 1:mxtx)	{
		s <- cs-1
		pmiss[it] <- pmiss[it]+cps[cs]*((pm_dt)^s)
#		cps <- prob_n_species_at_time_tm(s,1,p,q,dt)	# conditional probabilities of richnesses
#		pmiss[it] <- pmiss[it]+cps*((pm_dt)^s)
		}
	}
p_unsamp <- exp(sum(log(pmiss)))
}

prob_missing_paraclade_continuous_multibin <- function(pp,qq,rr,tt,prec)	{
# written 2016-02-29.  pp, qq, rr and tt are vectors giving origination, extinction & sampling rates in different bins
#	prec gives how precisely we divide up the time unit for "integrating" over interval
# output: vector giving the probability of a single lineage present at the outset of a bin failing to be sampled OR to hav
#	any descendants sampled given the origination, extinction & prservation rates in that bin and subsequent bins.  
bins <- length(pp)
p_miss_all <- vector(length=bins)
for (b in bins:1)	{
	p_miss_one_now <- prob_missing_paraclade_continuous(pp[b],qq[b],rr[b],tt[b],prec)
	if (b<bins)	{
		p_miss_one_next <- prob_missing_paraclade_continuous(pp[b+1],qq[b+1],rr[b+1],tt[b+1],prec)
		p_miss_surv <- prob_paraclade_extinction(1,pp[b],qq[b],tt[b])	# definitely miss survivors if there are none
		s <- ps <- mxps <- 0
		pfindafter <- 0
		while (ps>=(mxps/100000))	{
			s <- s+1
			ps <- prob_n_species_at_time_tm(s,1,pp[b],qq[b],tt[b])	# from Raup (1985)
			pm <- p_miss_one_next^s				# prob. of missing all s is P[miss one]^s
			p_miss_surv <- p_miss_surv+(ps*pm)	# prob. of missing all is ∑P[s | prior p,q,chrons] x P[0 finds | r,chrons]^s
			pfindafter <- pfindafter+ps*(1-pm)
			if (mxps<ps)	mxps <- ps
			}
		}	else p_miss_one_next <- 1.0

	p_miss_all[b] <- p_miss_one_now*p_miss_one_next
	}
return(p_miss_all)
}

gap_likelihoods_given_variable_sampling_over_time <- function(timescale, rr, gap_top, gap_bottom, prec=0.1)	{
stages <- length(rr)
st <- stages
while (timescale[st]<=gap_top)		st <- st-1
sb <- 1
while (timescale[sb]>gap_bottom)	sb <- sb+1
steps <- ceiling(abs(gap_bottom-gap_top)/prec)
time_steps <- create_rank_vector(steps)*prec

for (ss in st:sb)	{
	onset <- min(timescale[ss],gap_bottom)
	end <- max(timescale[ss+1],gap_top)
	gaps <- gap_likelihoods_per_sampling_interval(onset,end,rr[ss],prec)
	if (ss==st)	{
		gap_likelihoods <- gaps
		}	else	{
		gap_likelihoods <- c(gap_likelihoods,gaps*min(gap_likelihoods))	
		}
	}
names(gap_likelihoods) <- time_steps
return(gap_likelihoods)
}

gap_likelihoods_per_sampling_interval <- function(onset,end,r,prec)	{
substeps <- ceiling(abs(onset-end)/prec)
gapls <- vector(length=substeps)
for (chrons in 1:substeps)	gapls[chrons] <- exp(-(r*prec)*chrons)
return(gapls)
}

gap_loglikelihoods_per_sampling_interval <- function(onset,end,r,prec)	{
substeps <- ceiling(abs(onset-end)/prec)
gaplnls <- vector(length=substeps)
for (chrons in 1:steps)	gaplnls <- -(r*prec)*chrons
return(gaplnls)
}

# count the survivors
#numerare_concludes_reliquos  <- function(taxon_ranges,interval_names="")	{
count_survivors  <- function(taxon_ranges,interval_names="")	{
if (!is.data.frame(taxon_ranges))	{
	taxon_ranges <- data.frame(fa=taxon_ranges[,1],la=taxon_ranges[,2],row.names = rownames(taxon_ranges));
	}
ntaxa <- nrow(taxon_ranges);
st_a <- min(taxon_ranges[,1]);
st_z <- max(taxon_ranges[,2]);
survivors <- c();
survive_2 <- c();
for (st in st_a:(st_z-1))	{
#	xxxx <- subset(taxon_ranges,taxon_ranges$fa<=st);
#	yyyy <- subset(xxxx,xxxx$la>st);
#
#	survive_2 <- c(survive_2,nrow(yyyy));
	survivors <- c(survivors,sum((1:ntaxa)[taxon_ranges[,1]<=st] %in% (1:ntaxa)[taxon_ranges[,2]>st]));
	}

if (length(interval_names)>0)
	names(survivors) <- interval_names[st_a:(st_z-1)];
return(survivors);
}

# count the preceders
count_preceders  <- function(taxon_ranges,interval_names="")	{
if (!is.data.frame(taxon_ranges))	{
	taxon_ranges <- data.frame(fa=taxon_ranges[,1],la=taxon_ranges[,2],row.names = rownames(taxon_ranges));
	}
ntaxa <- nrow(taxon_ranges);
st_a <- min(taxon_ranges$fa);
st_z <- max(taxon_ranges$la);
preceders <- c();
#precies <- c();
for (st in (st_a+1):st_z)	{
	preceders <- c(preceders,
sum((1:ntaxa)[taxon_ranges$fa<st] %in% (1:ntaxa)[taxon_ranges$la>=st]));
	
#	xxxx <- subset(taxon_ranges,taxon_ranges$la>=st);
#	yyyy <- subset(xxxx,xxxx$fa<st);
#	precies <- c(precies,nrow(yyyy))
	}

if (length(interval_names)>0)
	names(preceders) <- interval_names[(st_a+1):st_z];
return(preceders);
}
#cbind(survivors,preceders)

count_shared <- function(samples_A,samples_B)	{
present_A <- present_B <- array(0,length(samples_A))
present_A[(1:length(samples_A))[samples_A>0]] <- 1
present_B[(1:length(samples_B))[samples_B>0]] <- 1
return(sum((present_A>0)*(present_B>0)))
}

### gives taxa shared between bin X and X+1
# use for extinction for X and origination for X+1
count_two_timers_after <- function(finds_per_bin)	{
bins <- ncol(finds_per_bin);
two_timers <- c();
for (b in 1:(bins-1))
	two_timers <- c(two_timers,count_shared(finds_per_bin[,b],finds_per_bin[,b+1]));
names(two_timers) <- colnames(finds_per_bin)[1:(bins-1)];
return(two_timers)
}

count_two_timers_before <- function(finds_per_bin)	{
bins <- ncol(finds_per_bin);
two_timers <- c();
for (b in 2:bins)
	two_timers <- c(two_timers,count_shared(finds_per_bin[,b-1],finds_per_bin[,b]));
names(two_timers) <- colnames(finds_per_bin)[2:bins];
return(two_timers)
}

### gives taxa shared between bin X and X+2
# use for extinction for X and origination for X+2
count_three_timers_after <- function(finds_per_bin)	{
bins <- ncol(finds_per_bin);
three_timers <- c();
for (b in 1:(bins-2))
	three_timers <- c(three_timers,count_shared(finds_per_bin[,b],finds_per_bin[,b+2]));
names(three_timers) <- colnames(finds_per_bin)[1:(bins-2)];
return(three_timers)
}

count_three_timers_before <- function(finds_per_bin)	{
bins <- ncol(finds_per_bin);
three_timers <- c();
for (b in 3:bins)
	three_timers <- c(three_timers,count_shared(finds_per_bin[,b-2],finds_per_bin[,b]));
names(three_timers) <- colnames(finds_per_bin)[3:bins];
return(three_timers)
}

### gives taxa shared between bin X and X+2 but absent in X+1
# use for extinction for X and origination for X+2 and sampling in X+1
count_three_timers_with_gap_after <- function(finds_per_bin)	{
bins <- ncol(finds_per_bin)
three_timers_with_gap <- vector(length=bins)
for (b in 1:(bins-2))
	three_timers_with_gap[b] <- sum((finds_per_bin[,b]>0)*(finds_per_bin[,b+1]==0)*(finds_per_bin[,b+2]>0))
return(three_timers_with_gap)
}

count_gaps_in_ranges <- function(taxon_ranges,finds_per_bin)	{
if (!is.data.frame(taxon_ranges))	{
	taxon_ranges <- data.frame(fa=as.numeric(taxon_ranges[,1]),fa=as.numeric(taxon_ranges[,2]));
	rownames(taxon_ranges) <- rownames(finds_per_bin);
	}
combined <- cbind(finds_per_bin,taxon_ranges);
bins <- ncol(finds_per_bin);
fa <- match("fa",colnames(combined));
la <- match("la",colnames(combined));
gappers <- vector(length=bins)
for (b in 2:(bins-1))	{
	relevant <- subset(combined,combined$fa<b);
	relevant <- subset(relevant,relevant$la>b);
	gappers[b] <- length(relevant[,b])-sum(relevant[,b]>0);
	}
names(gappers) <- colnames(finds_per_bin);
return(gappers)
}

#finds_per_bin <- finds_by_bin_test
count_gaps_per_taxon_in_ranges <- function(finds_per_bin,cutoff=1000)	{
# cutoff is the maximum gap size before we assume polyphyly
notu <- nrow(finds_per_bin)
nbins <- ncol(finds_per_bin)
first_last <- accio_classic_first_last_appearances(finds_per_bin)
taxon_gaps_per_bin <- array(0,dim=c(notu,nbins))
colnames(taxon_gaps_per_bin) <- colnames(finds_per_bin)
gap_length <- array(0,notu)
alterations <- 0
for (n in 1:notu)	{
	if (first_last[n,1]>0)	{
		gaps <- (first_last[n,1]:first_last[n,2])[finds_per_bin[n,first_last[n,1]:first_last[n,2]]==0]
		gap_length[n] <- length(gaps)
		if (gap_length[n]>0)	{
			taxon_gaps_per_bin[n,gaps] <- 1
			# now, remove overlong gaps
			if (gap_length[n]>=cutoff)	{
#				run <- 1
				gapsizes <- c()
				brk <- 1
				for (g in 2:gap_length[n]) if (gaps[g]>(gaps[g-1]+1))	brk <- brk+1
				gapborders <- array(0,dim=c(brk,2))
				gapborders[brk,2] <- gaps[gap_length[n]]
				brk <- 1
				gapborders[brk,1] <- gaps[1]
				for (g in 1:(gap_length[n]-1)) {
					if (gaps[g+1]>(gaps[g]+1))	{
						gapborders[brk,2] <- gaps[g]
						brk <- brk+1
						gapborders[brk,1] <- gaps[g+1]
						}
					}
				for (br in 1:brk)	{
					if (cutoff <= (1+gapborders[br,2]-gapborders[br,1]))	{
						taxon_gaps_per_bin[n,gapborders[br,1]:gapborders[br,2]] <- 0
						alterations <- alterations+1
						}
					}
#				g <- 2
#				brk <- 1
#				while (g <= gap_length[n])	{
#					gb1 <- g-1
#					while (gaps[g]==gaps[g-1]+1 && g<=gap_length[n])	{
#						run <- run+1
#						gb2 <- g
#						g <- g+1
#						}
#					gapborders[brk,] <- c(gb1,gb2)
#					gapsizes <- c(gapsizes,run)
#					}
#				for (gg in 1:length(gapsizes))	{
#					if (gapsizes[gg]>= cutoff)	{
#						taxon_gaps_per_bin[n,gaps[gapborders[gg,1]:gapborders[gg,2]]] <- 0
#						gap_length[n] <- gap_length[n]-(1+abs(gapborders[gg,2]-gapborders[gg,1]))
#						}
#					}
				}
			}
		}
	}
return(taxon_gaps_per_bin)
}

accio_classic_first_last_appearances <- function(finds_per_bin)	{
notu <- nrow(finds_per_bin)
nbins <- ncol(finds_per_bin)
first_last <- c()
for (n in 1:notu)	{
	bins <- (1:nbins)[finds_per_bin[n,]>0]
	if (length(bins)>0)	{
		first_last <- rbind(first_last,c(min(bins),max(bins)))
		} else {
		first_last <- rbind(first_last,c(0,0))
		}
	}
return(first_last)
}

#taxon_finds_per_bin <- finds_per_bin[257,]
gaps_in_taxon_range <- function(taxon_finds_per_bin)	{
bin_finds <- (1:length(taxon_finds_per_bin))[taxon_finds_per_bin>0];
gaps <- array(0,dim=length(taxon_finds_per_bin));
if (length(bin_finds)>0)	{
	bin_range <- min(bin_finds):max(bin_finds);
	bin_gaps <- bin_range[!bin_range %in% bin_finds]
	gaps[bin_gaps] <- 1;
	}
return(gaps);
}

#taxon_finds_per_bin <- finds_per_bin[2,]
finds_within_taxon_range <- function(taxon_finds_per_bin)	{
#print(taxon_finds_per_bin);
bin_finds <- (1:length(taxon_finds_per_bin))[taxon_finds_per_bin>0];
filler <- array(0,dim=length(taxon_finds_per_bin));
if (length(bin_finds)>2)	{
	bin_range <- (1+min(bin_finds)):(max(bin_finds)-1);
	filler[bin_range[taxon_finds_per_bin[bin_range]>0]] <- 1;
	}
return(filler);
}

numerare_gap_fillers_in_ranges <- function(taxon_ranges,finds_per_bin)	{
# taxon_ranges: a taxon x 2 matrix where column 1 gives 1st appearance bin & column 2
#	gives 2nd appearance bin. Bin 10 must precede Bin 11.
# finds_per_bin: taxon by bins matrix where finds_per_bin[n,b]≥1 means that the taxon n
#	is present in bin b. Occurrence numberes are fine: there is no need to convert them to 1 
#
#max(taxon_ranges[,2])
#colSums(finds_per_bin)
#min(taxon_ranges[,1])
#taxon_finds_per_bin <- finds_per_bin[1:25,];
#filler_matrix <- sapply(taxon_finds_per_bin,finds_within_taxon_range);
ntaxa <- nrow(taxon_ranges);
colnames(taxon_ranges) <- c("fa","la");
combined <- cbind(finds_per_bin,taxon_ranges);
mx_bins <- ncol(finds_per_bin);
bin_finds <- colSums(finds_per_bin);
mn_bins <- 1;
while(bin_finds[mn_bins]==0)	mn_bins <- mn_bins+1;
fa <- match("fa",colnames(combined));
la <- match("la",colnames(combined));
fillers <- vector(length=mx_bins);
for (b in (mn_bins+1):(mx_bins-1))	{
#	b <- b+1;
	appeared <- (1:ntaxa)[taxon_ranges[,1]<b];
	aragorns <- appeared[appeared %in% (1:ntaxa)[taxon_ranges[,2]>b]];
	fillers[b] <- sum(finds_per_bin[aragorns,b]>0);
#	if (ncol(relevant)>0)	fillers[b] <- sum(relevant[,b]>0);
	}
names(fillers) <- colnames(finds_per_bin)[1:mx_bins]
return(fillers)
}

numerare_gap_fillers_in_ranges_n_timers <- function(taxon_ranges,finds_per_bin,bins_pre=0,bins_aft=0)	{
# taxon_ranges: a taxon x 2 matrix where column 1 gives 1st appearance bin & column 2
#	gives 2nd appearance bin. Bin 10 must precede Bin 11.
# finds_per_bin: taxon by bins matrix where finds_per_bin[n,b]≥1 means that the taxon n
#	is present in bin b. Occurrence numberes are fine: there is no need to convert them to 1 
# bins_pre: look for range-throughs n bins earlier that are unsampled
#	this is used for origination
# bins_aft: look for range-throughs n bins late that are unsampled
#	this is used for extinction
taxon_ranges <- data.frame(fa=taxon_ranges[,1],la=taxon_ranges[,2],row.names = rownames(taxon_ranges));
ntaxa <- nrow(taxon_ranges);
combined <- cbind(finds_per_bin,taxon_ranges);
mx_bins <- ncol(finds_per_bin);
bin_finds <- colSums(finds_per_bin);
mn_bin_no <- min(taxon_ranges$fa);
mx_bin_no <- mn_bin_no+ncol(finds_per_bin)-1;

if (mn_bin_no > 1)	{
	match_test <- colSums(finds_per_bin);
	if (sum(match_test[1:mn_bin_no]==0)==(mn_bin_no-1))	{
		ad_hoc_strat_scale <- colnames(finds_per_bin);
		} else	{
		add_st <- (mn_bin_no-1)-sum(match_test[1:mn_bin_no]==0);
		ad_hoc_strat_scale <- c(rep("",add_st),colnames(finds_per_bin));
		}
	}
#b <- match("Hi",ad_hoc_strat_scale)
#fillers <- vector(length=mx_bins);
fillers <- vector(length=ncol(finds_per_bin));
names(fillers) <- colnames(finds_per_bin);
#if (bins_pre>0)	
for (b in (mn_bin_no+bins_pre):(mx_bin_no-bins_aft))	{
#	b <- b+1;
	bb <- match(ad_hoc_strat_scale[b],colnames(finds_per_bin));
	if (bins_pre > bins_aft)	{
#		bbb <- bb-bins_pre;
		bbb <- bb-1;
		} else if (bins_pre < bins_aft)	{
#		bbb <- bb+bins_aft;
		bbb <- bb+1;
		} else	{
		bbb <- bb;
		}
	already_born <- subset(taxon_ranges,taxon_ranges$fa<=(b-bins_pre));
	spanners <- subset(already_born,already_born$la>(b+bins_aft));
	# get species that span the bins
	three_timers <- match(rownames(spanners),rownames(finds_per_bin));
#	finds_per_bin[three_timers,(b-1):(b+2)]
#	ntx <- (1:ntaxa)[rownames(finds_per_bin) %in% rownames(already_born[three_timers,])];
	fillers[bb] <- sum(finds_per_bin[three_timers,bbb]>0);
	}
return(fillers)
}

numerare_gaps_in_ranges_n_timers <- function(taxon_ranges,finds_per_bin,bins_pre=1,bins_aft=1)	{
# taxon_ranges: a taxon x 2 matrix where column 1 gives 1st appearance bin & column 2
#	gives 2nd appearance bin. Bin 10 must precede Bin 11.
# finds_per_bin: taxon by bins matrix where finds_per_bin[n,b]≥1 means that the taxon n
#	is present in bin b. Occurrence numberes are fine: there is no need to convert them to 1 
#
taxon_ranges <- data.frame(fa=taxon_ranges[,1],la=taxon_ranges[,2],row.names = rownames(taxon_ranges));
ntaxa <- nrow(taxon_ranges);
combined <- cbind(finds_per_bin,taxon_ranges);
mx_bins <- ncol(finds_per_bin);
bin_finds <- colSums(finds_per_bin);
mn_bin_no <- min(taxon_ranges$fa);
mx_bin_no <- mn_bin_no+ncol(finds_per_bin)-1;

if (mn_bin_no > 1)	{
	match_test <- colSums(finds_per_bin);
	if (sum(match_test[1:mn_bin_no]==0)==(mn_bin_no-1))	{
		ad_hoc_strat_scale <- colnames(finds_per_bin);
		} else	{
		add_st <- (mn_bin_no-1)-sum(match_test[1:mn_bin_no]==0);
		ad_hoc_strat_scale <- c(rep("",add_st),colnames(finds_per_bin));
		}
	}

gappers <- vector(length=ncol(finds_per_bin));
names(gappers) <- colnames(finds_per_bin);
for (b in (mn_bin_no+bins_pre):(mx_bin_no-bins_aft))	{
#	b <- b+1;
	bb <- match(ad_hoc_strat_scale[b],colnames(finds_per_bin));
	if (bins_pre > bins_aft)	{
		bbb <- bb-bins_pre;
		} else if (bins_pre < bins_aft)	{
		bbb <- bb+bins_aft;
		} else	{
		bbb <- bb;
		}
	already_born <- subset(taxon_ranges,taxon_ranges$fa<=(b-bins_pre));
	not_yet_dead <- subset(already_born,already_born$la>(b+bins_aft));
	three_timers <- match(rownames(not_yet_dead),rownames(finds_per_bin));
#	ntx <- (1:ntaxa)[rownames(finds_per_bin) %in% rownames(already_born[three_timers,])];
	gappers[bb] <- sum(finds_per_bin[three_timers,bbb]==0);
	}

return(gappers)
}

numerare_gaps_in_ranges <- function(taxon_ranges,finds_per_bin)	{
# taxon_ranges: a taxon x 2 matrix where column 1 gives 1st appearance bin & column 2
#	gives 2nd appearance bin. Bin 10 must precede Bin 11.
# finds_per_bin: taxon by bins matrix where finds_per_bin[n,b]≥1 means that the taxon n
#	is present in bin b. Occurrence numberes are fine: there is no need to convert them to 1 
#
#max(taxon_ranges[,2])
#colSums(finds_per_bin)
#min(taxon_ranges[,1])
#taxon_finds_per_bin <- finds_per_bin[1:25,];
#filler_matrix <- sapply(taxon_finds_per_bin,finds_within_taxon_range);
ntaxa <- nrow(taxon_ranges);
colnames(taxon_ranges) <- c("fa","la");
combined <- cbind(finds_per_bin,taxon_ranges);
mx_bins <- ncol(finds_per_bin);
bin_finds <- colSums(finds_per_bin);
mn_bins <- 1;
while(bin_finds[mn_bins]==0)	mn_bins <- mn_bins+1;
fa <- match("fa",colnames(combined));
la <- match("la",colnames(combined));
gappers <- vector(length=mx_bins);
for (b in (mn_bins+1):(mx_bins-1))	{
#	b <- b+1;
	appeared <- (1:ntaxa)[taxon_ranges[,1]<b];
	aragorns <- appeared[appeared %in% (1:ntaxa)[taxon_ranges[,2]>b]];
	gappers[b] <- sum(finds_per_bin[aragorns,b]==0);
#	if (ncol(relevant)>0)	fillers[b] <- sum(relevant[,b]>0);
	}
names(gappers) <- colnames(finds_per_bin)[1:mx_bins]
return(gappers)
}

accio_observed_richness <- function(finds_per_bin)	{
# finds_per_bin: taxon by bins matrix where finds_per_bin[n,b]≥1 means that the taxon n
#	is present in bin b. Occurrence numberes are fine: there is no need to convert them to 1 
return(colSums(finds_per_bin>0))
}

accio_synoptic_richness <- function(taxon_ranges,interval_names="")	{
# taxon_ranges: a taxon x 2 matrix where column 1 gives 1st appearance bin & column 2
#	gives 2nd appearance bin. Bin 10 must precede Bin 11.
#ntaxa <- nrow(taxon_ranges);
bins <- max(taxon_ranges[,2]);
mnbin <- min(taxon_ranges[taxon_ranges[,1]!=0,1]);
jacks <- c();

for (b in mnbin:bins)	jacks <- c(jacks,sum((taxon_ranges[,1]<=b)*(taxon_ranges[,2]>=b)));

if (length(interval_names)>0)	names(jacks) <- interval_names[mnbin:bins];
return(jacks)
}

# written 2019-08-15 to simplify life
# modified 2019-08-31 to simplify life
accio_stratigraphic_ranges_from_sampled_in_bin <- function(finds_per_bin)	{
#present_in_bin <- 1*(finds_per_bin>=0.5);
ttl_bins <- 1:ncol(finds_per_bin);
ttl_taxa <- nrow(finds_per_bin);
taxon_ranges <- c();
for (tx in 1:ttl_taxa)	{
	p_i_b <- 1*(finds_per_bin[tx,]>=0.5);
	if (sum(p_i_b)==0)
		p_i_b[match(max(finds_per_bin[tx,]),finds_per_bin[tx,])] <- 1;
	p_i_b <- ttl_bins*p_i_b;
	p_i_b <- p_i_b[p_i_b>0];
	taxon_ranges <- rbind(taxon_ranges,c(min(p_i_b),max(p_i_b)));
	}
rownames(taxon_ranges) <- rownames(finds_per_bin);
return(taxon_ranges);
}

# written 2019-08-15 to simplify life
accio_synoptic_richness_from_sampled_in_bin <- function(finds_per_bin)	{
taxon_ranges <- accio_stratigraphic_ranges_from_sampled_in_bin(finds_per_bin);
synoptic_richness <- vector(length=max(max(taxon_ranges),ncol(finds_per_bin)));
names(synoptic_richness) <- colnames(finds_per_bin);
for (tx in 1:nrow(taxon_ranges))	{
	synoptic_richness[taxon_ranges[tx,1]:taxon_ranges[tx,2]] <- synoptic_richness[taxon_ranges[tx,1]:taxon_ranges[tx,2]]+1;
	}
return(synoptic_richness);
}

count_singletons_per_bin <- function(taxon_ranges,interval_names="")	{
# taxon_ranges: a taxon x 2 matrix where column 1 gives 1st appearance bin & column 2
#	gives 2nd appearance bin. Bin 10 must precede Bin 11.
bins <- max(taxon_ranges[,2])
mnbin <- min(taxon_ranges[,1])
sgtpepper <- vector(length=bins)
for (b in mnbin:bins)
	sgtpepper[b] <- sum((taxon_ranges[,1]==b)*(taxon_ranges[,2]==b))
if (length(interval_names)>=bins)
	names(sgtpepper) <- interval_names[mnbin:bins];
return(sgtpepper)
}

count_very_singletons_per_bin <- function(finds_per_bin)	{
return(colSums((rowSums(finds_per_bin)==1)*(finds_per_bin>0)))
}

# taxa that cross lower boundary and upper boundary
count_range_throughs <- function(taxon_ranges,interval_names="")	{
bins <- max(taxon_ranges[,2]);
mnbin <- min(taxon_ranges[,1]);
foote_bt <- vector(length=bins);
for (b in (mnbin+1):(bins-1))	{
	foote_bt[b] <- sum((taxon_ranges[,1]<b)*(taxon_ranges[,2]>b));
	}
if (length(interval_names)>=bins)
	names(foote_bt) <- interval_names[mnbin:bins];
return(foote_bt)
}

# taxa that cross lower boundary and upper boundary
accio_foote_bt <- function(taxon_ranges,interval_names="",TRUNC=FALSE)	{
if (!is.data.frame(taxon_ranges))
	taxon_ranges <- data.frame(taxon_ranges);
colnames(taxon_ranges) <- c("fa","la");
bins <- max(taxon_ranges[,2]);
mnbin <- min(taxon_ranges[,1]);
foote_bt <- vector(length=bins);
names(foote_bt) <- interval_names[1:bins];
#for (b in (mnbin+1):(bins-1))	{
#	preceded <- subset(taxon_ranges,taxon_ranges$fa<b);
#	spanned <- subset(preceded,preceded$la>b);
#	foote_bt[b] <- nrow(spanned);
#	}
#foote_bta <- vector(length=bins);
for (b in (mnbin+1):(bins-1))
	foote_bt[b] <- sum((taxon_ranges[,1]<b)*(taxon_ranges[,2]>b));
#if (length(interval_names)>=bins)
#	names(foote_bt) <- interval_names[mnbin:bins];
if (TRUNC)
	foote_bt <- foote_bt[mnbin:length(foote_bt)];
return(foote_bt)
}

# taxa that cross lower boundary but not upper boundary
accio_foote_bL <- function(taxon_ranges,interval_names="")	{
taxon_ranges <- data.frame(taxon_ranges);
colnames(taxon_ranges) <- c("FA","LA")
bins <- max(taxon_ranges[,2]);
mnbin <- min(taxon_ranges[,1])-1;
relevant_ranges <- subset(taxon_ranges,taxon_ranges$FA<taxon_ranges$LA);
foote_bL <- hist(relevant_ranges$LA,breaks=((mnbin-1):bins),plot=FALSE)$counts;
if (length(interval_names)>=bins)
	names(foote_bL) <- interval_names[mnbin:bins];
#for (tx in 1:notu)	{
#	if (taxon_ranges$FA[tx] < taxon_ranges$LA[tx])
#		foote_bL[taxon_ranges$LA[tx]] <- foote_bL[taxon_ranges$LA[tx]];
#	}
return(foote_bL);
}

# taxa that cross upper boundary but not lower boundary
accio_foote_Ft <- function(taxon_ranges,interval_names="")	{
taxon_ranges <- data.frame(taxon_ranges);
colnames(taxon_ranges) <- c("FA","LA")
bins <- max(taxon_ranges[,2]);
mnbin <- min(taxon_ranges[,1])-1;
relevant_ranges <- subset(taxon_ranges,taxon_ranges$FA<taxon_ranges$LA);
foote_Ft <- hist(relevant_ranges$FA,breaks=((mnbin-1):bins),plot=FALSE)$counts;
if (length(interval_names)>=bins)
	names(foote_Ft) <- interval_names[mnbin:bins];
#bins <- max(taxon_ranges[,2]);
#mnbin <- min(taxon_ranges[,1]);
#foote_Ft <- vector(length=bins);
#for (b in mnbin:(bins-1))	{
#	foote_Ft[b] <- sum((taxon_ranges[,1]==b)*(taxon_ranges[,2]>b));
#	}
return(foote_Ft);
}

# function to calculate probability of observing n of N shared taxa give sampling
#	x probability of N survivors given extinction/origination/beta & S1 original taxa
prob_observing_n_shared_given_turnover_sampling_and_hypothesized_shared <- function(n,S1,S2,freq_turn,pmiss)	{
# S1 = observed in bin
# S2 = hypothesized in other bin
# n = observed shared
# freq_turn=expected n=(S1-S2)S1
# pfind = prob finding taxon
#return(dbinom(S1-S2,S1,freq_turn)*dbinom(n,S2,pfind))
return(dbinom(S1-S2,S1,freq_turn)*pmiss^(S2-n))
}

# 
likelihood_turnover_rate_given_sampling <- function(rate,pmiss,S1,two_timer,gap_filler,continuous)	{
# rate: per-lineage rate
#	poisson if continuous==TRUE; binomial if continuous==FALSE
# pmiss: probability that an unobserved taxon is present-but-unsampled rather than non-existent
# S1: observed richness in "test" interval
# two_timer: shared richness
# gap_filler: number of taxa inferred to be present because they are there before and after
# continuous: TRUE for continuous turnover, FALSE for discrete turnover
if (continuous==TRUE)	{
	freq <- Poisson_rate_to_probability(rate)
	}	else {
	freq <- rate
	}	# get expected frequency of shared taxa
hS <- seq(two_timer+gap_filler,S1,by=1)	# range of possible shared taxa; max=S1, min=observed+directly inferred
pt <- prob_observing_n_shared_given_turnover_sampling_and_hypothesized_shared(n=two_timer,S1=S1,S2=hS,freq_turn=freq,pmiss=pmiss)

return(max(10^-320,sum(pt)))
}

accio_best_turnover_given_sampling <- function(pmiss,S1,two_timer,gap_filler,continuous)	{
# rate gives the proportion NOT shared between bins
#		originations or extinctions for historical data
# pmiss: typical probability of failing to sample a taxon;
# S1: observed standing richness (no unsampled range-throughs);
# two_timer: number shared with the prior interval;
# gap_filler: number of unsampled taxa spanning a gap (synoptic - observed richness);
# continuous: if T, then continuous diversification assumed; otherwise, pulsed.
if (continuous==TRUE)	{
	if ((two_timer+gap_filler)>0)	{
		rate <- -log((two_timer+gap_filler)/S1);
		}	else {
		rate <- -log(pmiss/S1);
		}
	}	else {
	rate <- 1-(two_timer/S1);
	}
# we expect to sample (1-pmiss) of the taxa.
minrate <- 0;
maxrate <- rate;	# no point in considering a rate higher than face value
if (maxrate==0)
	maxrate <- 0.1;	# added 2019-09-08

cl <- list(fnscale=-1);
if ((two_timer+gap_filler)<S1)	{
	accio <- optim(rate,fn=likelihood_turnover_rate_given_sampling,method="L-BFGS-B",pmiss=pmiss,S1=S1,two_timer=two_timer,gap_filler=gap_filler,continuous=continuous,lower=0,upper=maxrate,control=cl)
	best_diversification <- max(0,accio$par)
	L_best_diversification <- accio$value
	}	else {
	best_diversification <- 0.0
	L_best_diversification <- 1.0
	}
output <- data.frame(ML_Rate=as.numeric(best_diversification),MlnL_Rate=as.numeric(log(L_best_diversification)));
return(output)
}

# added 2019-08-15
accio_diversification_likelihoods <- function(occr_per_bin,prob_miss_interval,first_bin=1)	{
# added 2019-08-15
# occr_per_bin: finds per bin.  Numbers can be fractions reflecting uncertaing
# prob_miss_interval: probability of failing to find a taxon in an interval;
# first_bin: which bin to begin tallying diversity; defaults to 1
n_intervals <- ncol(occr_per_bin);
interval_names <- colnames(occr_per_bin);

sampled_in_bin <- 1*(occr_per_bin>=0.5);
tnotu <- nrow(occr_per_bin);
missed_taxa <- (1:tnotu)[rowSums(sampled_in_bin)==0];
for (mt in missed_taxa)	{
	if (max(occr_per_bin[mt,])<0)
		sampled_in_bin[mt,match(max(occr_per_bin[mt,]),occr_per_bin[mt,])] <- 1;
	}

#cbind(as.character(stg_slices[1:lst_sl,1]),sampled_in_bin_spc)
taxon_gaps_per_slice <- count_gaps_per_taxon_in_ranges(finds_per_bin=sampled_in_bin,cutoff=5);
rownames(taxon_gaps_per_slice) <- rownames(sampled_in_bin);
two_timer_after <- count_two_timers_after(sampled_in_bin);
three_timer_after <- count_three_timers_after(sampled_in_bin);
three_not_two_timer_after <- count_three_timers_with_gap_after(sampled_in_bin);
two_timer_pre <- c(0,two_timer_after[1:(n_intervals-1)]);
three_timer_pre <- c(0,0,three_timer_after[1:(n_intervals-2)]);
three_not_two_timer_pre <- c(0,0,three_not_two_timer_after[1:(n_intervals-2)]);

diversification_rates <- data.frame(ML_Orig=as.numeric(rep(0,n_intervals)),lnL_ML_Orig=as.numeric(rep(0,n_intervals)),LB_Sup_Orig=as.numeric(rep(0,n_intervals)),UB_Sup_Orig=as.numeric(rep(0,n_intervals)),
									ML_Extn=as.numeric(rep(0,n_intervals)),lnL_ML_Extn=as.numeric(rep(0,n_intervals)),LB_Sup_Extn=as.numeric(rep(0,n_intervals)),UB_Sup_Extn=as.numeric(rep(0,n_intervals)));
rownames(diversification_rates) <- interval_names;

observed_in_bin <- colSums(sampled_in_bin);
### start here!!!!
for (b in first_bin:length(interval_names))	{
	# get the number of taxa observed in bin X
	# get number of those observed in prior bin
	# ask the probability of getting that given the sampling for that bin & extinction rates
	update <- paste("Doing Best Overall Turnover for ",interval_names[b],sep="");
	print(update)
	a <- b-1
	c <- b+1
	S1 <- sampled_in_bin[b]; 	#fixed 2019-09-07
	observed <- (1:tnotu)[sampled_in_bin[,b]>0];
	if (b>(first_bin+1))	{
		pmiss <- prob_miss_interval[a]
		two_timer <- two_timer_pre[b]
		gap_filler <- sum(taxon_gaps_per_slice[observed,b-1])
		dynamic <- accio_best_turnover_given_sampling(pmiss,S1,two_timer,gap_filler,continuous=TRUE);
		diversification_rates$ML_Orig[b]=dynamic$ML_Rate;
		diversification_rates$lnL_ML_Orig[b]=dynamic$MlnL_Rate;
		diversification_rates[b,1:2] <- accio_best_turnover_given_sampling(pmiss,S1,two_timer,gap_filler,continuous=TRUE)
		if (diversification_rates$ML_Orig[b]>0)	{
			poss_rates <- seq(diversification_rates[b,1]/10,10*diversification_rates[b,1],by=diversification_rates[b,1]/100);
			} else	{
			poss_rates <- seq(0,1,by=0.05);
			}
		error_bars <- log(mapply(likelihood_turnover_rate_given_sampling,poss_rates,pmiss,S1,two_timer=two_timer_pre[b],gap_filler=three_not_two_timer_pre[b],continuous=TRUE))-diversification_rates[b,2]
		diversification_rates$LB_Sup_Orig[b] <- min(poss_rates[error_bars>-1]);
		diversification_rates$UB_Sup_Orig[b] <- max(poss_rates[error_bars>-1]);
		}
	if (b<(last_bin_backdrop-1))	{
		pmiss <- prob_miss_interval[c];
		two_timer <- two_timer_after[b];
		gap_filler <- sum(taxon_gaps_per_slice[observed,b+1]);
		dynamic <- accio_best_turnover_given_sampling(pmiss,S1,two_timer,gap_filler,continuous=TRUE);
		diversification_rates$ML_Extn[b]=dynamic$ML_Rate;
		diversification_rates$lnL_ML_Extn[b]=dynamic$MlnL_Rate;
		if (diversification_rates$ML_Extn[b]>0)	{
			poss_rates <- seq(diversification_rates$ML_Extn[b]/5,5*diversification_rates$ML_Extn[b],by=diversification_rates$ML_Extn[b]/100);
			} else	{
			poss_rates <- seq(0,1,by=0.05);
			}
		error_bars <- log(mapply(likelihood_turnover_rate_given_sampling,poss_rates,pmiss,S1,two_timer=two_timer_after[b],gap_filler=three_not_two_timer_after[b],continuous=TRUE))-diversification_rates$lnL_ML_Extn[b];
		diversification_rates$LB_Sup_Extn[b] <- min(poss_rates[error_bars>-1]);
		diversification_rates$UB_Sup_Extn[b] <- max(poss_rates[error_bars>-1]);
		}
	}
return(diversification_rates);
}

# added 2019-08-15
# modified 2019-08-31
accio_diversification_per_ma_likelihoods <- function(occr_per_bin,prob_miss_interval,time_scale,first_bin=1,last_bin=MAXNO)	{
# added 2019-08-15
# occr_per_bin: finds per bin.  Numbers can be fractions reflecting uncertaing
# prob_miss_interval: probability of failing to find a taxon in an interval;
# time_scale: time scale giving interval names matching those in column headings for occr_per_bin & onsets / ends of those intervals
# first_bin: which bin to begin tallying diversity; defaults to 1
n_intervals <- ncol(occr_per_bin);
interval_names <- colnames(occr_per_bin);

#accio_synoptic_richness_from_sampled_in_bin(finds_per_bin = occr_per_bin);
sampled_in_bin <- 1*(occr_per_bin>=0.5);
tnotu <- nrow(occr_per_bin);
missed_taxa <- (1:tnotu)[rowSums(sampled_in_bin)==0];
for (mt in missed_taxa)	{
	if (max(occr_per_bin[mt,])<0)
		sampled_in_bin[mt,match(max(occr_per_bin[mt,]),occr_per_bin[mt,])] <- 1;
	}

#cbind(as.character(stg_slices[1:lst_sl,1]),sampled_in_bin_spc)
taxon_gaps_per_slice <- count_gaps_per_taxon_in_ranges(finds_per_bin=sampled_in_bin,cutoff=5);
rownames(taxon_gaps_per_slice) <- rownames(sampled_in_bin);
two_timer_after <- count_two_timers_after(sampled_in_bin);
three_timer_after <- count_three_timers_after(sampled_in_bin);
three_not_two_timer_after <- count_three_timers_with_gap_after(sampled_in_bin);
two_timer_pre <- c(0,two_timer_after[1:(n_intervals-1)]);
three_timer_pre <- c(0,0,three_timer_after[1:(n_intervals-2)]);
three_not_two_timer_pre <- c(0,0,three_not_two_timer_after[1:(n_intervals-2)]);

diversification_rates <- data.frame(ML_Orig=as.numeric(rep(0,n_intervals)),lnL_ML_Orig=as.numeric(rep(0,n_intervals)),LB_Sup_Orig=as.numeric(rep(0,n_intervals)),UB_Sup_Orig=as.numeric(rep(0,n_intervals)),
									ML_Extn=as.numeric(rep(0,n_intervals)),lnL_ML_Extn=as.numeric(rep(0,n_intervals)),LB_Sup_Extn=as.numeric(rep(0,n_intervals)),UB_Sup_Extn=as.numeric(rep(0,n_intervals)));
rownames(diversification_rates) <- interval_names;

observed_in_bin <- colSums(sampled_in_bin);
if (last_bin>length(interval_names))
	last_bin <- length(interval_names);

### start here!!!!
for (b in first_bin:last_bin)	{
	# get the number of taxa observed in bin X
	# get number of those observed in prior bin
	# ask the probability of getting that given the sampling for that bin & extinction rates
	update <- paste("Doing Best Overall Turnover for ",interval_names[b],sep="");
	interval_duration <- abs(time_scale$ma_lb[match(interval_names[b],time_scale$interval)]-time_scale$ma_ub[match(interval_names[b],time_scale$interval)]);
	print(update);
	a <- b-1;
	c <- b+1;
	S1 <- observed_in_bin[b];	# fixed 2019-09-07
	observed <- (1:tnotu)[sampled_in_bin[,b]>0];
#	if (b>fst_sl && med_sampling[a]>0)	{
	if (b>(first_bin+1))	{
		pmiss <- prob_miss_interval[a];
		two_timer <- two_timer_pre[b]
		gap_filler <- sum(taxon_gaps_per_slice[observed,b-1]);
		dynamic <- accio_best_turnover_given_sampling(pmiss,S1,two_timer,gap_filler,continuous=TRUE);
		diversification_rates$ML_Orig[b]=dynamic$ML_Rate/interval_duration;
		diversification_rates$lnL_ML_Orig[b]=dynamic$MlnL_Rate;
		if (diversification_rates$ML_Orig[b]>0)	{
			poss_rates <- seq(dynamic$ML_Rate/5,5*dynamic$ML_Rate,by=dynamic$ML_Rate/100);
			} else	{
			poss_rates <- seq(0,1,by=0.05);
			}
		error_bars <- log(mapply(likelihood_turnover_rate_given_sampling,poss_rates,pmiss,S1,two_timer=two_timer_pre[b],gap_filler=three_not_two_timer_pre[b],continuous=TRUE))-diversification_rates$lnL_ML_Orig[b];
		diversification_rates$LB_Sup_Orig[b] <- min(poss_rates[error_bars>-1])/interval_duration;
		diversification_rates$UB_Sup_Orig[b] <- max(poss_rates[error_bars>-1])/interval_duration;
		}
	if (b<(last_bin-1))	{
		pmiss <- prob_miss_interval[c];
		two_timer <- two_timer_after[b];
		gap_filler <- sum(taxon_gaps_per_slice[observed,b+1]);
		dynamic <- accio_best_turnover_given_sampling(pmiss,S1,two_timer,gap_filler,continuous=TRUE);
		diversification_rates$ML_Extn[b]=dynamic$ML_Rate/interval_duration;
		diversification_rates$lnL_ML_Extn[b]=dynamic$MlnL_Rate;
		if (diversification_rates$ML_Extn[b]>0)	{
			poss_rates <- seq(dynamic$ML_Rate/5,5*dynamic$ML_Rate,by=dynamic$ML_Rate/100);
			} else	{
			poss_rates <- seq(0,1,by=0.05);
			}
		error_bars <- log(mapply(likelihood_turnover_rate_given_sampling,poss_rates,pmiss,S1,two_timer=two_timer_after[b],gap_filler=three_not_two_timer_after[b],continuous=TRUE))-diversification_rates$lnL_ML_Extn[b];
		diversification_rates$LB_Sup_Extn[b] <- min(poss_rates[error_bars>-1])/interval_duration;
		diversification_rates$UB_Sup_Extn[b] <- max(poss_rates[error_bars>-1])/interval_duration;
		}
	} #help
return(diversification_rates);
}

likelihood_turnover_rate_given_sampling_old <- function(rate,pmiss,S1,two_timer,gap_filler,continuous)	{
# rate: per-lineage rate
#	poisson if continuous==TRUE; binomial if continuous==FALSE
# pmiss: probability that an unobserved taxon is present-but-unsampled rather than non-existent
# S1: observed richness in "test" interval
# hS2: hypothesized richness, which must be at least known shared but no greater than S1
#	this usually will be a vector, with hS2[1] the minimum possible shared richness
# continuous: TRUE for continuous turnover, FALSE for discrete turnover
if (continuous==TRUE)	{
	prob <- Poisson_rate_to_probability(rate)
	}	else {
	prob <- rate
	}

if (prob>0)	{
#	debug <- paste("S1=",S1,"; ",length(probs),sep="")
#	print(debug)
	# prob of N taxa being there given rate x prob of observing X
	pturn <- sum(dbinom(S1-hS2,S1,prob)*dbinom(two_timer,S1-hS2,(1-pmiss)))
	}	else {
	# probability of observing X taxa if all of them were there
	pturn <- dbinom(two_timer,S1,(1-pmiss))
	}
return(sum(pturn))
}

accio_best_turnover <- function(median_sampling,stdev_mag,S1,two_timer,gap_filler,ncoll,continuous)	{
mxx <- 0.99*pnorm(-log(median_sampling)/log(stdev_mag))
pfind <- 1-integrate(prob_missing_given_lognormal,lower=0,upper=mxx,median_sampling=median_sampling,stdev_mag=stdev_mag,ncoll=ncoll)$value
# rate gives the proportion NOT shared between bins
#		originations or extinctions for historical data
if (continuous==TRUE)	{
#	rate <- init_rate <- 1-exp(log(two_timer)-log(S1))
	rate <- -log((two_timer+gap_filler)/S1)
	minrate <- log(min(two_timer/pfind,S1)/S1)/(-4/3)
#	minrate <- rate/2
	maxrate <- rate
	}	else {
	rate <- init_rate <- 1-(two_timer/S1)
	minrate <- 1-(min(two_timer/pfind,S1)/S1)
	maxrate <- rate
	}
hS2 <- seq((two_timer+gap_filler),S1,by=1)	# possible number of shared taxa

cl <- list(fnscale=-1)
#optim(c(br,mag),fn=loglikelihood_lognormal_occupancy_rates,method="L-BFGS-B",oS=oS,hS=hS,observed=observed,ncoll=ncoll,lower=c(min(sc),min(mm)),upper=c(max(sc),max(mm)),control=cl)
bt <- optim(rate,fn=likelihood_turnover_rate_given_sampling,method="L-BFGS-B",pmiss=1-pfind,S1=S1,two_timer=two_timer,gap_filler=gap_filler,continuous=TRUE,lower=minrate,upper=maxrate,control=cl)
return(c(bt$par,log(bt$value)))
}

#Xbt <- 10
#XFtbL <- 5
#Poisson_rate_to_probability(.4054651)
boundary_crosser_turnover <- function(Xbt,XFtbL)	{
# Xbt is the number of taxa ranging through an interval
# XFtbL is either XFt or XbL (i.e., those entering/exiting interval, but not leaving/entering interval)
return(-log(Xbt/(Xbt+XFtbL)))
}

assign_interval_rates_to_finer_time_scale_old <- function(bin_rates,bin_onsets,prec)	{
slices <- 1+(max(bin_onsets)-min(bin_onsets))/prec
ts <- bin_onsets
for (b in 2:length(bin_onsets))	ts[b] <- ts[b]+ts[b-1]
fine_rates <- vector(length=slices)
b <- 1
chrons <- bin_onsets[b]
for (i in 1:slices)	{
	if (chrons>=ts[b+1])	b <- b+1
	fine_rates[i] <- bin_rates[b]
	chrons <- chrons+prec
	}
return(fine_rates)
}

accio_first_differences <- function(time_series)	{
tt <- length(time_series)-1
first_diffs <- vector(length = tt)
for (i in 1:tt)
	first_diffs[i] <- time_series[i+1]-time_series[i]
return(first_diffs)
}

### richness estimators
chao2 <- function(abundance)	{
ntaxa <- sum(abundance>0)
q1 <- sum(abundance==1);
q2 <- sum(abundance==2);
x <- ntaxa+((q1*q1)/(2*(q2+1))) - (((q1*q2))/((2*(q2+1)*(q2+1))))
S <- round(x)
return (S)
}

### taxon sampling estimates crude
chao2_prop_taxon_sampling <- function(abundance)	{
ntaxa <- sum(abundance>0)
q1 <- sum(abundance==1);
q2 <- sum(abundance==2);
fr <- 1-((((q1*q1)/(2*(q2+1))) - (((q1*q2))/((2*(q2+1)*(q2+1)))))/(ntaxa+((q1*q1)/(2*(q2+1))) - (((q1*q2))/((2*(q2+1)*(q2+1))))))
return (fr)
}

jack5_prop_richness <- function(abundance)	{
ntaxa <- sum(abundance>0);
observed <- c();
for (nn in 1:5)
	observed <- c(observed,sum(abundance==nn));
ss <- 0
for (i in 1:length(observed))	ss <- ss+(i*observed[i])
fr <- ntaxa/(ntaxa + (observed[1]*(((5*ss)-15)/ss)) - observed[2]*((10*(ss*ss)-(70*ss)+125)/(ss*(ss-1))) + observed[3]*(((10*(ss^3))-(120*(ss^2))+(485*ss)-660))/(ss*(ss-1)*(ss-2)) - observed[4]*((ss-4)^4)/(ss*(ss-1)*(ss-2)*(ss-3)) + observed[5]*((ss-5)^5)/(ss*(ss-1)*(ss-2)*(ss-3)*(ss-4)));
return(fr);
}

coverage_good_1965 <- function(Xi)	{
# Good's U
# Xi: number of specimens/occurrences for each taxon i
O <- sum(Xi);		# ∑ occurrences
s1 <- sum(Xi==1);	# No. of singletons
return(1-(s1/O));
}

coverage_alroy_2010 <- function(Xi,single_pub_occ)	{
# Alroy 2010 Supplement
# Xi: number of specimens/occurrences for each taxon i
# single_pub: taxa known from only 1 publication
# richest_coll: occurrences from biggest (richest) collection
return(1-((sum(Xi) - (single_pub_occ + Xi[1]))/sum(Xi)));
}

coverage_chao_jost_2012 <- function(Xi,m)	{
# Chao & Jost 2012 eq. 4b
# m: smallest sample size
# Xi: number of specimens/occurrences for each taxon i
n <- sum(Xi)
S <- length(Xi)
summation <- 0
for (i in 1:S)	summation <- summation+((Xi[i]/n)*(choose(n-Xi[i],m)/choose(n,m)))
return(1-summation)
}

estimated_species_richness_at_sample_size_m <- function(So,m,n)	{

}


#find_data_set <- subset(find_data,find_data$fuzz==0);
#find_data_set_20 <- subset(find_data_set,find_data_set$bin_lb==20);
#taxa_in_20 <- unique(find_data_set_20$accepted_name);
#observations <- c();
#for (tt in 1:length(taxa_in_20))	{
#	species_finds_20 <- subset(find_data_set_20,find_data_set_20$accepted_name==taxa_in_20[tt]);
#	observations <- c(observations,length(unique(species_finds_20$collection_no)));
#	}
#observations <- sort(observations,decreasing=TRUE);
#ttl_collections_20 <- length(unique(find_data_set_20$collection_no));
#observations_a <- Kropotkin_Distribution(S=length(observations),N=sum(observations))
#observations_z <- Koch_Distribution(S=length(observations),N=sum(observations),C=ttl_collections_20)
#squares(observations)
#squares(observations=observations_a)
#squares(observations=observations_z)

squares <- function(observations)	{
So <- length(observations);
s1 <- sum(observations==1);
N2 <- sum(observations)^2
n2 <- observations^2
return(So+((s1^2)*sum(n2))/(N2-s1*So));
}

multiton_rescale <- function(observations)	{
N <- sum(observations)
So <- length(observations)
s1 <- sum(observations==1)
m <- So-s1
return((N/(s1+1))*(m/So))
}

capture_mark_recapture_with_distributed_sampling <- function(presences,sampling,sampling_dist="LogNormal",timers=3,metric)	{
# presences: matrix giving number of presences in each bin
#	will be converted to presence/absence
# sampling <- array of numbers describing sampling distribution
#		if just a vector, then it must be uniform!
# sampling_dist: the distribution model for sampling (default is lognormal)
# timers: number of consecutive intervals over which we look
# metric: origination or extinction.
presabs <- presences
presabs[presabs>0] <- 1
bins <- ncol(presabs)
notu <- nrow(presabs)
}

#### Simple Cal3

# 1-Phi given Didier et al.'s 2017 Quadratic
Prob_Missing_Clade_Didier <- function(p,q,r)	{
ab <- c((p+q+r+sqrt(((p+q+r)^2)-(4*p*q)))/(2*p),(p+q+r-sqrt(((p+q+r)^2)-(4*p*q)))/(2*p))
phi <- min(ab)
return(phi)
}

# Phi given Didier et al.'s 2017 Quadratic
Prob_Sampling_Clade_Didier <- function(p,q,r)	{
ab <- c((p+q+r+sqrt(((p+q+r)^2)-(4*p*q)))/(2*p),(p+q+r-sqrt(((p+q+r)^2)-(4*p*q)))/(2*p))
phi <- 1-min(ab)
return(phi)
}

# conditional probability of a paraclade ancestor being extant given paraclade is extant
prob_extant_paraclade_ancestor_extant <- function(lambda,mu,ts)	{
pexx <- prob_paraclade_survival(a=1,p=lambda,q=mu,chrons=ts)
pex <- exp(-mu*ts)
return(pex/pexx)
}

prob_lineage_unsampled_from_t1_to_t2_variable_r <- function(r,bs,prec=0.10)	{
rr <- assign_interval_rates_to_finer_time_scale(bin_rates = r,bin_spans = bs,prec)
bins <- length(rr)
p_ghost <- array(0,dim=c(bins,bins))
for (b1 in 1:(bins-1))	{
	p_ghost[b1,b1] <- 1.0
	for (b2 in (b1+1):bins)	{
		p_ghost[b1,b2] <- p_ghost[b2,b1] <- exp(sum(-1*prec*rr[b1:(b2-1)]))			# prob. 0 finds
		}
	}
p_ghost[bins,bins] <- 1.0
return(p_ghost)
}

# probability of paraclade surviving from divergence to sampling/splitting AND having zero samples
prob_range_extension_given_varying_diversification <- function(p,q,r,bs,prec=0.10)	{
#pp <- assign_interval_rates_to_finer_time_scale(bin_rates = p,bin_spans = bs,prec)
#qq <- assign_interval_rates_to_finer_time_scale(bin_rates = q,bin_spans = bs,prec)
tt <- rep(prec,bins)

prob_branch_surviving <- accio_prob_paraclade_survivorship_in_time_slices_given_shifting_rates(p,q,bs,prec=0.1,max_S=15)
rr <- assign_interval_rates_to_finer_time_scale(bin_rates = r,bin_spans = bs,prec)
bins <- length(rr)
p_ghost <- array(0,dim=c(bins,bins))
for (b1 in 1:(bins-1))	{
	for (b2 in (b1+1):bins)	{
		p_ghost[b1,b2] <- p_ghost[b2,b1] <- exp(sum(-1*prec*rr[b1:b2]))			# prob. 0 finds
		}
	}
prob_branch_from_divergence <- p_ghost*prob_branch_surviving
return(p_ghost*pexx)
}

# conditional probability of a paraclade ancestor being extant given paraclade is extant with shifting rates
prob_extant_paraclade_ancestor_extant_varying_diversification <- function(p,q,bs,prec=0.1,max_S=15)	{
# p: origination rates from different intervals
# q: extinction rates from same intervals
# bs: time-spans of intervals
pp <- assign_interval_rates_to_finer_time_scale(bin_rates = p,bin_spans = bs,prec)
qq <- assign_interval_rates_to_finer_time_scale(bin_rates = q,bin_spans = bs,prec)
tt <- rep(prec,bins)
# get probablity of paraclade survival over time interval with rate shifts
xxx <- accio_prob_richness_shifting_rates(pp,qq,tt,max_S)
pexx <- 1-xxx[1,length(xxx[1,])]	# prob. paraclade survival
pex <- exp(sum(-1*prec*qq))			# prob. progenitor survival
return(pex/pexx)
}

# get the probability of standing richness given shifts in diversification after clade onset
accio_prob_richness_shifting_rates_by_bin_partial <- function(bb,pp,qq,tt,max_S=50)	{
tb <- length(pp)
dd <- cbind(pp[bb:tb],qq[bb:tb])
dvrs <- c(pp[bb],qq[bb])
ttt <- array(tt[bb],dim=1)
j <- 1
for (i in 2:(tb+1-bb))	{
	if (sum(dd[i,]==dd[i-1,])==2)	{
		ttt[j] <- ttt[j]+tt[i]
		} else	{
		dvrs <- rbind(dvrs,c(pp[i],qq[i]))
		ttt <- c(ttt,tt[i])
		j <- j+1
		}
	}

n <- 0:max_S
p_s_outset <- array(0,dim=1+max_S)
p_s_outset[2] <- 1
ts <- seq(prec,ttt[1],by=prec)
p_rich_at_time_slice <- sapply(ts,accio_prob_richness_time_slices_within_interval,max_S,p=dvrs[1,1],q=dvrs[1,2],p_s_outset,prec=0.1)

for (b in 2:nrow(dvrs))	{
	p_s_outset <- p_rich_at_time_slice[,ncol(p_rich_at_time_slice)]
	ts <- seq(prec,ttt[b],by=prec)
	prats <- sapply(ts,accio_prob_richness_time_slices_within_interval,max_S,p=dvrs[b,1],q=dvrs[b,2],p_s_outset,prec=0.1)
	p_rich_at_time_slice <- cbind(p_rich_at_time_slice,prats)
	}
return(p_rich_at_time_slice)
#p_rich_at_time_slice[,ncol(p_rich_at_time_slice)]
#sapply(n,prob_n_species_at_time_tm,a=1,p,q,chrons=32.5)
}

accio_prob_richness_time_slices_within_interval <- function(ts,max_S=50,p,q,p_s_outset,prec=0.1)	{
a <- 1:max_S
ns <- 0:max_S
jjjj <- sapply(a,accio_paraclade_richness_probs,ns,p,q,chrons=ts)
hhhh <- rep(0,max_S+1)
hhhh[1] <- 1		# if it starts at 0, then it is guaranteed to stay at 0
jjjj <- cbind(hhhh,jjjj)

for (i in 1:(max_S+1)) jjjj[,i] <- jjjj[,i]*p_s_outset[i]
prob_rich <- rowSums(jjjj)
return(prob_rich)
}

accio_paraclade_richness_probs <- function(a,ns,p,q,chrons)	{
n <- ns
pr_s_n <- sapply(n,prob_n_species_at_time_tm,a,p,q,chrons)
return(pr_s_n)
}

# get the probability of standing richness given shifts in diversification after clade onset
accio_prob_richness_shifting_rates <- function(pp,qq,tt,max_S=50)	{
bins <- length(tt)
n <- 0:max_S
prob_S <- array(0,dim=c(max_S+1,bins))
# prob of 0…n species surviving to the next bin
#pr_s <- sapply(n,prob_n_species_at_time_tm,a=1,p=pp[b],q=qq[b],chrons=tt[b])
prob_S[2,1] <- 1
#for (b in 2:bins)	{
b <- 2
while (b<=bins)	{
	pr_s <- rep(0,1+max_S)
	prob_S[,b-1]
	for (a in 0:max_S)	{
		pr_s_a <- sapply(n,prob_n_species_at_time_tm,a,p=pp[b-1],q=qq[b-1],chrons=tt[b-1])
		pra <- prob_S[a+1,b-1]
		pr_s <- pr_s+pra*pr_s_a
#		print(c(sum(pr_s),sum(prob_S[1:(a+1),b-1])))
		}
	prob_S[,b] <- pr_s
	b <- b+1
	}
return(prob_S)	# probability of 0…max_S species in each time slice given rates
}

# get the probability of a paraclade going extinct over intervals with different origination & sampling
accio_prob_paraclade_extinction_in_time_slices_given_shifting_rates <- function(p,q,bs,prec=0.1)	{
pp <- assign_interval_rates_to_finer_time_scale(bin_rates = p,bin_spans = bs,prec)
qq <- assign_interval_rates_to_finer_time_scale(bin_rates = q,bin_spans = bs,prec)
bins <- length(pp)
tt <- rep(prec,bins)
prob_no_survivors <- array(0,dim=c(bins,bins-1))

for (b in 1:(bins-1))	{
	ps <- pp[b:bins]
	qs <- qq[b:bins]
	tms <- tt[b:bins]
	xxx <- accio_prob_richness_shifting_rates(ps,qs,tms,max_S=15)
	prob_no_survivors[b:bins,b] <- xxx[1,]
	}
return(prob_no_survivors)
}

# get the probability of a paraclade surviving over intervals with different origination & sampling
accio_prob_paraclade_survivorship_in_time_slices_given_shifting_rates <- function(p,q,bs,prec=0.1,max_S=15)	{
pp <- assign_interval_rates_to_finer_time_scale(bin_rates = p,bin_spans = bs,prec)
qq <- assign_interval_rates_to_finer_time_scale(bin_rates = q,bin_spans = bs,prec)
bins <- length(pp)
tt <- rep(prec,bins)
prob_survivors <- array(0,dim=c(bins,bins))

b <- 1:(bins-1)
#prob_extinct <- sapply(b,accio_prob_paraclade_extinct_shifting_rates_from_b_to_B,pp,qq,tt,max_S)
prob_survivors <- sapply(b,accio_prob_paraclade_extant_shifting_rates_from_b_to_B,pp,qq,tt,max_S)
for (b in 1:(bins-1))	{
	if (b%%25==0)	print(b)
	xxx <- accio_prob_paraclade_extinct_shifting_rates_from_b_to_B(b,pp,qq,tt,max_S)
	1-xxx[1,]
	ps <- pp[b:bins]
	qs <- qq[b:bins]
	tms <- tt[b:bins]
	xxx <- accio_prob_richness_shifting_rates(ps,qs,tms,max_S)
	prob_survivors[b:bins,b] <- 1-xxx[1,]
	prob_survivors[b,b:bins] <- 1-xxx[1,]
	}
return(prob_survivors)
}

accio_prob_paraclade_extinct_shifting_rates_from_b_to_B <- function(b,pp,qq,tt,max_S=50)	{
if (b%%25==0)	print(b)
bins <- length(pp)
xxx <- accio_prob_richness_shifting_rates(pp[b:bins],qq[b:bins],tt[b:bins],max_S)
prob_ex <- xxx[1,]
if (b>1)	prob_ex <- c(rep(0,b-1),prob_ex)
return(prob_ex)
}

accio_prob_paraclade_extant_shifting_rates_from_b_to_B <- function(b,pp,qq,tt,max_S=50)	{
if (b%%25==0)	print(b)
bins <- length(pp)
xxx <- accio_prob_richness_shifting_rates(pp[b:bins],qq[b:bins],tt[b:bins],max_S)
prob_srv <- 1-xxx[1,]
if (b>1)	prob_srv <- c(rep(1,b-1),prob_srv)
return(prob_srv)
}

### let's work it backwards....
# get p_surv or p_extn from bin[N,M…1];
#	at bin N-1, use P[0…max_S] at final slice to put conditional probs on outcomes afterwards
prob_paraclade_survives_from_t1_orig_to_t2 <- function(p,q,bs,prec=0.1,max_S=20)	{
tb <- length(p)
dd <- cbind(p,q)
dvrs <- c(p[1],q[1])
ttt <- array(bs[1],dim=1)
j <- 1
for (i in 2:tb)	{
	if (sum(dd[i,]==dd[i-1,])==2)	{
		ttt[j] <- ttt[j]+bs[i]
		} else	{
		dvrs <- rbind(dvrs,c(p[i],q[i]))
		ttt <- c(ttt,bs[i])
		j <- j+1
		}
	}
lmp_b <- nrow(dvrs)		# lumped bins

tts <- sum(bs)/prec		# total time slices
ttss <- round(sum(ttt[1:(lmp_b-1)])/prec,0)

# get the probability that a lineage present at any point in the first B-1 stages has 0…Max_S surivovors at the end of that stage
prob_successors <- array(0,dim=c(max_S+1,ttss))
n <- 0:max_S
for (b in 1:(lmp_b-1))	{
	k <- round(sum(ttt[1:b])/prec,0)
	if (b>1)	{
		j <- round(1+sum(ttt[1:(b-1)])/prec,0)
		} else	{
		j <- 1
		}
	for (h in j:k)	{
		prob_successors[,h] <- sapply(n,prob_n_species_at_time_tm,a=1,p=dvrs[b,1],q=dvrs[b,2],chrons=prec*(1+k-h))
		}
	}

# work backwards to get the probability that a species appearing in any time slice is extant in each future time slice
prob_survival_t1_to_t2 <- array(0,dim=c(tts,tts))
for (b in lmp_b:1)	{
	k <- round(sum(ttt[1:b])/prec,0)
	if (b>1)	{
		j <- round(1+sum(ttt[1:(b-1)])/prec,0)
		} else	{
		j <- 1
		}
	ttl_slc <- round(ttt[b]/prec,0)
	chrons <- seq(ttt[b]-prec,prec,by=-1*prec)
	pr_extn <- sapply(chrons,prob_paraclade_extinction,a=1,p=dvrs[b,1],q=dvrs[b,2])
	pr_surv <- 1-pr_extn
	pr_surv <- c(pr_surv,1)
	tbs <- length(pr_surv)
	
	for (h in k:j)	{
		# get prob of surviving from onset of cell a to the onset of cell b
		prob_survival_t1_to_t2[h:j,h] <- prob_survival_t1_to_t2[h,h:j] <- pr_surv[tbs:(1+(k-h))]
		if (b>1)	{
			for (bb in (b-1):1)	{
				kk <- round(sum(ttt[1:bb])/prec,0)
				if (bb>1)	{
					jj <- round(1+sum(ttt[1:(bb-1)])/prec,0)
					} else	{
					jj <- 1
					}
				for (hh in kk:jj)	{
					p_extn <- 0
					for (a in 0:max_S)	{
						p_extn <- p_extn+(prob_successors[a+1,hh]*(1-prob_survival_t1_to_t2[(kk+1),h])^a)
#						print(c(a,prob_successors[a+1,hh],(1-prob_survival_t1_to_t2[(kk+1),h])^a,prob_successors[a+1,kk]*(1-prob_survival_t1_to_t2[(kk+1),h])^a))
						}
					prob_survival_t1_to_t2[h,hh] <- prob_survival_t1_to_t2[hh,h] <- 1-p_extn
					}
				}
			}
		}
	}
### now I have prob successors! This should let me fill in backwards....
### remember, I only need P[extn | n species, chrons, q, p]
return(prob_survival_t1_to_t2)
}

accio_extinction_rate_given_median_duration <- function(median_duration)	{
mu <- -log(0.5)/median_duration
return(mu)
}

probability_of_lineage_extinction_given_extinction_rate_and_time <- function(mu,chrons)	{
return(1-exp(-mu*chrons))
}

probability_of_n_of_N_lineage_extinctions_given_extinction_rate_and_time <- function(n,N,mu,chrons,exact=FALSE)	{
q <- probability_of_lineage_extinction_given_extinction_rate_and_time(mu,chrons)
if (exact)	{
	# get exact probability of n of N extinctions
	return(dbinom(n,N,q))
	} else	{
	# get probability of such extreme (or more extreme) results
	m <- 0:n
	pnq <- sum(dbinom(m,N,q))
	if (pnq>0.5)	{
		return(1-pnq)
		} else	{
		return(pnq)
		}
	}
}

rate_of_extinction_given_n_of_N_lineage_extinctions_over_time <- function(n,N,chrons)	{
mu <- log(1+n/N)/chrons
return (mu)
}

# get the overlap of two temporal duratins, with older durations having lower numbers
accio_range_overlap <- function(fa1,la1,fa2,la2)	{
if ((fa1 <= fa2 && fa1 > la2) || (fa2 <= fa1 && fa2 > la1)) {
#if (fa1 > la2 || fa2 > la1)	{
	overlap <- c(min(fa1,fa2),max(la1,la2));
	}	else	{
	overlap <- c(0,0);
	}
return (overlap);
}

accio_range_overlaps_for_multiple_cases <- function(strat_ranges)	{
if (strat_ranges$fa1 > strat_ranges$la2 || strat_ranges$fa2 > strat_ranges$la1)	{
	overlap <- c(min(strat_ranges$fa1,strat_ranges$fa2),max(strat_ranges$la1,strat_ranges$la2));
	}	else	{
	overlap <- c(0,0);
	}
return (overlap);
}

#### routines to get the likelihood of an origination or extinction rate
## 		written initially to review a paper about treating diversification rates
##		as models of increasing complexity

loglike_turnover_given_observed_n_and_sampling <- function(diverse_rate,n,S1,S2,pmiss,interval_duration=1)	{
# S1 = observed in bin
# S2 = minimum shared in other bin
# n = observed shared in other bin
# freq_turn=expected n=(S1-S2)S1
# pfind = prob finding taxon
# interval_duration: if 1, then this is just a per-stage rate
freq_turn <- Poisson_rate_to_probability(diverse_rate*interval_duration);
if(sum(dbinom(S1-(S2:S1),S1,freq_turn)*dbinom(n,(S2:S1),1-pmiss)) > 0)	{
#	return(log(sum(dbinom(S1-(S2:S1),S1,freq_turn)*dbinom(n,(S2:S1),1-pmiss))));
	return(logprob_observing_n_shared_given_turnover_and_sampling(n,S1,S2,pmiss,freq_turn));
	} else	{
	return(log(10^-200));
	}
}

logprob_observing_n_shared_given_turnover_and_sampling <- function(n,S1,S2,pmiss,freq_turn)	{
# same as above, but 
# S1 = observed in bin
# S2 = minimum shared in other bin
# n = observed shared in other bin
# freq_turn=expected n=(S1-S2)S1
# pfind = prob finding taxon
if(sum(dbinom(S1-(S2:S1),S1,freq_turn)*dbinom(n,(S2:S1),1-pmiss)) > 0)	{
	return(log(sum(dbinom(S1-(S2:S1),S1,freq_turn)*dbinom(n,(S2:S1),1-pmiss))));
	} else	{
	return(log(10^-200));
	}
}

#rate_info <- c(0.41882725,0.46263251,0.00100000,1.31520510,0.00100000,0.00100000,0.05885117,0.15336673,0.35726634,0.62500000,0.66666667,0.66666667,0.68143587,0.75000000,5.00000000)
#rate_info <- c(0.82010782,0.68499268,0.74587571,0.44251958,0.05692436,0.00100000,0.15934332,0.03558896,0.20185329,0.81783872,0.00100000,0.23175507,0.00100000,0.00100000,0.12260035,0.08378606,0.77601261,0.19505764,0.25069343,0.27851132,0.30632922,0.36196501,0.38978290,0.41760080,0.44541869,0.47323659,0.50105449,0.79166667,0.61132607)
#rate_info <- c(0.71861383,0.45858597,0.50473791,0.28286488,0.28137610,0.24101617,0.04462178,0.01688710,0.23319399,0.11278280,0.00100000,0.19686850,0.33986633,0.00100000,0.00100000,0.07604161,0.77601261,0.19505764,0.28819440,0.27851132,0.30632922,0.36196501,0.38978290,0.40858271,0.44541869,0.47323659,0.50105449,0.75600250,0.61032607)
#rate_info <- c(0.07425486,0.11132716,0.18547175,0.22254405,0.25961635,0.33376095,0.37083324,0.40790554,0.48205014,0.51912243,0.55619473,0.59326703,0.63033933,0.66741163,0.79266667,0.81570082,1.00000000)
turgio_stage_breaks_fraction <- function(divides,min_divides,max_divides)	{
divides <- sort(divides);
too_high <- (1:length(divides))[divides>max_divides];
if (length(too_high)>0)	divides[too_high] <- max_divides[too_high];
too_low <- (1:length(divides))[divides<min_divides];
if (length(too_low)>0)	divides[too_low] <- min_divides[too_low];
return(divides);
}

turgio_stage_breaks_discrete <- function(divides,max_divide=max(divides))	{
divides <- sort(divides);
divides[divides<1] <- 2;	# first divide cannot come before 2nd bin
attempts <- 0;
while (length(unique(divides)) < length(divides) && attempts < 100)	{
	fix_these <- c();
	for (i in 1:length(divides))
		if (sum(divides==divides[i])>1)	fix_these <- c(fix_these,i);
	prob_children <- unique(divides[fix_these]);
	openings <- (2:max_divide)[!(2:max_divide) %in% divides];
	for (i in 1:length(prob_children))	{
		adjacency <- abs(prob_children[i]-openings);
		divides[match(prob_children[i],divides)] <- openings[match(min(adjacency),adjacency)];
		divides <- sort(divides);
		openings <- (2:max_divide)[!(2:max_divide) %in% divides];
		}
	attempts <- attempts+1;
	}
if (attempts==100)
	divides[1:(length(divides)-1)] <- 2:length(divides)
return(divides);
}

log_likelihood_single_diversification_rate <- function(diverse_rate,data_for_p_or_q,interval_durations="")	{
bins <- nrow(data_for_p_or_q);
if (length(interval_durations)==1)	interval_durations <- rep(1,bins);
#print(diverse_rate);
lgl_rate <- c();
for (b in 1:bins)	{
	freq_turn <- Poisson_rate_to_probability(diverse_rate*interval_durations[b]);
	n <- data_for_p_or_q$n[b];
	S1 <- data_for_p_or_q$S1[b];
	S2 <- data_for_p_or_q$S2[b];
	pmiss <- data_for_p_or_q$pmiss[b];
	lgl_rate <- c(lgl_rate,logprob_observing_n_shared_given_turnover_and_sampling(n,S1,S2,pmiss,freq_turn));
	}
return(sum(lgl_rate))
}

reparifarge_rate_info <- function(rate_info, nbins)	{
### routine to take rate_info from an optimization routine and present it normallyy
sep_rates <- ceiling(length(rate_info)/2);
a <- 1+sep_rates;
divides <- sort(rate_info[a:length(rate_info)]);
#1/divides[1]
#nbins <- length(rate_info)/2;
max_divides <- (nbins-((length(divides)-1):0))/nbins;
min_divides <- (0.01+(1:length(divides)))/nbins;
divides <- turgio_stage_breaks_fraction(divides,min_divides,max_divides);

divides <- ceiling(divides*nbins);
#divides <- c(divides,nbins+1);
if (sep_rates==nbins)	{
	divides <- (1:nbins)+1;
	} else	{
#	divides <- turgio_stage_breaks_discrete(divides,max_divide=nbins+1);
	divides <- turgio_stage_breaks_discrete(divides,max_divide=nbins);
	}
divides <- c(1,divides,nbins+1);
div_rates <- vector(length=nbins);
#cbind(divides[1:(sep_rates+1)],divides[2:(sep_rates+2)],rate_info[1:sep_rates])
for (br in 1:sep_rates)	{
	div_rates[divides[br]:(divides[br+1]-1)] <- rate_info[br];
	}
output <- list(div_rates,divides);
names(output) <- c("diversification_rates","divides");
return(output)
}

log_likelihood_multiple_diversification_rates <- function(rate_info,data_for_p_or_q,interval_durations)	{
#print(rate_info);
bins <- nrow(data_for_p_or_q);
if (length(interval_durations)==1)	interval_durations <- rep(1,bins);

key_info <- reparifarge_rate_info(rate_info,nbins=bins);
div_rates <- key_info$diversification_rates;
divides <- key_info$divides;
sep_rates <- ceiling(length(rate_info)/2);
div_rates <- div_rates*interval_durations;
lgl_rate <- c();
#b <- 0;
for (b in 1:bins)	{
#	b <- b+1;
	n <- data_for_p_or_q$n[b];
	S1 <- data_for_p_or_q$S1[b];
	S2 <- data_for_p_or_q$S2[b];
	pmiss <- data_for_p_or_q$pmiss[b];
	freq_turn <- Poisson_rate_to_probability(div_rates[b]);
	lgl_rate <- c(lgl_rate,logprob_observing_n_shared_given_turnover_and_sampling(n,S1,S2,pmiss,freq_turn));
	}
#} else	{
#	lgl_rate <- rep(-10^50,bins);
#	}
#print(c(sum(lgl_rate),rate_info[1:sep_rates],divides));
#print(sum(lgl_rate));
return(sum(lgl_rate));	
}

list_rate_shifts <- function(interval_rates,include_first=FALSE)	{
nb <- length(interval_rates);
rate_shifts <- c();
if (include_first)	rate_shifts <- 1;
for (ir in 2:length(interval_rates))	if (interval_rates[ir]!=interval_rates[ir-1])	rate_shifts <- c(rate_shifts,ir);
return(rate_shifts);
}

count_individual_rates <- function(interval_rates)	{
nb <- length(interval_rates);
rate_shifts <- 1;
for (ir in 2:length(interval_rates))	if (interval_rates[ir]!=interval_rates[ir-1])	rate_shifts <- rate_shifts+1;
return(rate_shifts);
}

accio_temporally_partitioned_diversification_rates_bottom_up <- function(data_for_p_or_q,max_rate,min_rate,ttl_data,interval_durations="")	{
cl <- list(fnscale=-1);
nbins <- nrow(data_for_p_or_q);
if (length(interval_durations)<nbins)
	interval_durations <- rep(1,nbins);

diverse_rate <- sum(saturated_diversification_model[,1]*interval_durations)/sum(interval_durations);
ttl_rates <- 1;
brl <- optim(diverse_rate,fn=log_likelihood_single_diversification_rate,method="L-BFGS-B",data_for_p_or_q=data_for_p_or_q,interval_durations=interval_durations,lower=min_rate,upper=max_rate,control=cl);
#brl <- optim(diverse_rate,fn=log_likelihood_single_diversification_rate,method="L-BFGS-B",data_for_p_or_q=data_for_p_or_q,lower=min_rate,upper=max_rate,control=cl);
diversification_models <- rep(brl$par,nrow(saturated_diversification_model));
diversification_model_loglikelihoods <- brl$value;
diversification_model_AICcs <- modified_AIC(lnL=brl$value,k=ttl_rates,n=ttl_data);
names(diversification_models) <- rownames(data_for_p_or_q);

if (nbins>99)	{
	m_n <- paste("k=00",1,sep="")
	} else if (nbins>9)	{
	m_n <- paste("k=0",1,sep="")
	} else	{
	m_n <- paste("k=",1,sep="")
	}
model_names <- m_n
max_attempts <- 10;x
while (ttl_rates < (nbins-1))	{
	try <- no_improvement <- 0;
	mxlgl <- -10^100;
	ttl_rates <- ttl_rates+1;
	rate_divs <- ttl_rates-1;
	max_divides <- (nbins-((rate_divs-1):0))/nbins;
	min_divides <- (0.01+(1:rate_divs))/nbins;
	max_rates <- c(rep(max_rate,ttl_rates),max_divides);
	min_rates <- c(rep(min_rate,ttl_rates),min_divides);
	if (ttl_rates > 2)
		last_splits <- splits;
	while (no_improvement<10)	{
		if (ttl_rates>2 && (try %% 3)==0)	{
			poss_splits <- (2:nbins)[!(2:nbins) %in% last_splits];
			pspl <- length(poss_splits);
			if (max_attempts > (pspl*3))
				max_attempts <- pspl*3;
			x <- poss_splits[order(runif(pspl))]
			init_divides <- sort(c(last_splits[2:rate_divs],x[1]));
			#### START HERE PETE!!!
			} else {
			init_divides <- (2:nbins)[order(runif(nbins-1))];
			}
		discrete_init_divides <- c(1,sort(init_divides[1:rate_divs]),nbins+1);
		divides <- sort(init_divides[1:rate_divs])/nbins;
	
		diverse_rates <- c();
		for (rr in 1:ttl_rates)	{
			a <- discrete_init_divides[rr];
			z <-  discrete_init_divides[rr+1]-1;
			diverse_rates <- c(diverse_rates,sum(saturated_diversification_model[a:z,1]*interval_durations[a:z])/sum(interval_durations[a:z]));
			}
		rate_info <- c(diverse_rates,divides);
		brl <- optim(rate_info,fn=log_likelihood_multiple_diversification_rates,method="L-BFGS-B",data_for_p_or_q=data_for_p_or_q,interval_durations=interval_durations,lower=rep(min_rate,ttl_rates),upper=rep(max_rate,ttl_rates),control=cl);
		
		if (brl$value > mxlgl)	{
			mbrl <- brl;
			mxlgl <- mbrl$value;
			no_improvement <- 0;
			} else	{
			no_improvement <- no_improvement+1;
			}
		}

#	if (mbrl$value < diversification_model_loglikelihoods[ttl_rates])	{
#		print("help!")
#		}
	
	best_params <- reparifarge_rate_info(rate_info=mbrl$par,nbins);
	this_model <- best_params$diversification_rates;
	splits <- best_params$divides;

	diversification_models <- cbind(diversification_models,this_model);
	diversification_model_loglikelihoods <- c(diversification_model_loglikelihoods,mbrl$value);
	diversification_model_AICcs <- c(diversification_model_AICcs,modified_AIC(lnL=mbrl$value,k=ttl_rates,n=ttl_data));
#	examined_splits_org <- rbind(examined_splits_org,c(splits,rep(0,nbins-ttl_rates)));
	if (nbins>99)	{
		if (ttl_rates<10)	{
			m_n <- paste("k=00",ttl_rates,sep="");
			} else if (ttl_rates<100)	{
			m_n <- paste("k=0",ttl_rates,sep="");
			} else	{
			m_n <- paste("k=",ttl_rates,sep="");
			}
		} else if (nbins>9)	{
		if (ttl_rates<10)	{
			m_n <- paste("k=0",ttl_rates,sep="");
			} else	{
			m_n <- paste("k=",ttl_rates,sep="");
			}
		} else	{
		m_n <- paste("k=",ttl_rates,sep="");
		}
	model_names <- c(model_names,m_n);
	colnames(diversification_models) <- names(diversification_model_loglikelihoods) <- names(diversification_model_AICcs) <- model_names;
	print(c(ttl_rates,round(diversification_model_AICcs,1)));
#	diversification_models
	}

saturated_diversification_model <- c();
for (stg in 1:nbins)	{
	S1 <- data_for_p_or_q$S1[stg];			# synoptic richness
	S2 <- data_for_p_or_q$S2[stg];			# minimum holdovers (= Nt = Nbt - NbL or Nb = Nbt - NFt)
	n <- data_for_p_or_q$n[stg];			# observed holdovers
	pmiss <- data_for_p_or_q$pmiss[stg];	# prob. of missing
	diverse_rate <- boundary_crosser_turnover(Xbt=S1,XFtbL = max(1,S2))/interval_durations[stg];
	br_max <- optim(diverse_rate,fn=loglike_turnover_given_observed_n_and_sampling,method="L-BFGS-B",S1=S1,S2=S2,n=n,pmiss=pmiss,interval_duration=interval_durations[stg],lower=min_rate,upper=max_rate,control=cl);
	saturated_diversification_model <- rbind(saturated_diversification_model,c(br_max$par,br_max$value));
	}
sat_splits <- 1;
for (stg in 2:nbins)	{
	if (saturated_diversification_model[stg,1]!=saturated_diversification_model[stg-1,1])
		sat_splits <- c(sat_splits,stg)
	}

if (length(sat_splits)==nbins)	{
	model_names <- c(model_names,paste("k=",nbins,sep=""));
	diversification_models <- cbind(diversification_models,saturated_diversification_model[,1]);
	diversification_model_loglikelihoods <- c(diversification_model_loglikelihoods,sum(saturated_diversification_model[,2]));
	diversification_model_AICcs <- c(diversification_model_AICcs,modified_AIC(lnL=diversification_model_loglikelihoods[nbins],k=ttl_rates,n=ttl_data));
	colnames(diversification_models) <- names(diversification_model_loglikelihoods) <- names(diversification_model_AICcs) <- model_names;
	} else	{
	max_complexity <- length(sat_splits);
	diversification_models[,max_complexity] <- saturated_diversification_model[,1]
	diversification_model_loglikelihoods[max_complexity] <- sum(saturated_diversification_model[,2])
	diversification_model_AICcs[max_complexity] <- modified_AIC(lnL=(sum(saturated_diversification_model[,2])),k=max_complexity,n=ttl_data);
	if (ncol(diversification_models)>max_complexity)	{
		diversification_models <- diversification_models[,1:max_complexity];
		diversification_model_loglikelihoods <- diversification_model_loglikelihoods[1:max_complexity];
		diversification_model_AICcs <- diversification_model_AICcs[1:max_complexity];
		}
	}

output <- list(diversification_models,diversification_model_loglikelihoods,diversification_model_AICcs);
names(output) <- c("models","loglikelihoods","AICc");
return(output);
}

accio_temporally_partitioned_diversification_rates_top_down <- function(data_for_p_or_q,max_rate,min_rate,ttl_data,interval_durations="")	{
cl <- list(fnscale=-1);
nbins <- nrow(data_for_p_or_q);
if (length(interval_durations)<nbins)
	interval_durations <- rep(1,nbins);

saturated_diversification_model <- c();
for (stg in 1:nbins)	{
	S1 <- data_for_p_or_q$S1[stg];			# synoptic richness
	S2 <- data_for_p_or_q$S2[stg];			# minimum holdovers (= Nt = Nbt - NbL or Nb = Nbt - NFt)
	n <- data_for_p_or_q$n[stg];			# observed holdovers
	pmiss <- data_for_p_or_q$pmiss[stg];	# prob. of missing
	diverse_rate <- boundary_crosser_turnover(Xbt=S1,XFtbL = max(1,S2))/interval_durations[stg];
	br_max <- optim(diverse_rate,fn=loglike_turnover_given_observed_n_and_sampling,method="L-BFGS-B",S1=S1,S2=S2,n=n,pmiss=pmiss,interval_duration=interval_durations[stg],lower=min_rate,upper=max_rate,control=cl);
	saturated_diversification_model <- rbind(saturated_diversification_model,c(br_max$par,br_max$value));
	}
rownames(saturated_diversification_model) <- rownames(data_for_p_or_q);
model_names <- paste("k=",nbins)
#for (i in 1:10)	
#	print(loglike_turnover_given_observed_n_and_sampling(i/10,S1=S1,S2=S2,n=n,pmiss=pmiss,interval_duration=interval_durations[stg]));

splits <- list_rate_shifts(interval_rates = saturated_diversification_model[,1],include_first=FALSE);
ttl_rates <- count_individual_rates(interval_rates = saturated_diversification_model[,1]);
last_splits <- init_splits <- splits;
model_names <- paste("k=",ttl_rates,sep="");

diversification_models <- saturated_diversification_model[,1];
diversification_model_loglikelihoods <- sum(saturated_diversification_model[,2]);
diversification_model_AICcs <- modified_AIC(lnL=diversification_model_loglikelihoods,k=ttl_rates,n=ttl_data);

max_attempts <- 10;

while (ttl_rates > 2)	{
	try <- no_improvement <- 0;
	mxlgl <- -10^100;
#	old_split_range <- 2:(ttl_rates+1);
#	new_split_range <- 2:ttl_rates;
	ttl_rates <- ttl_rates-1;
	rate_divs <- ttl_rates-1;
	max_divides <- (nbins-((rate_divs-1):0))/nbins;
	min_divides <- (0.01+(1:rate_divs))/nbins;
	max_rates <- c(rep(max_rate,ttl_rates),max_divides);
	min_rates <- c(rep(min_rate,ttl_rates),min_divides);
	k_rate_lnls <- c();
	while (no_improvement<20)	{
		poss_split_removed <- last_splits[last_splits %in% (2:nbins)];
		x <- poss_split_removed[order(runif(ttl_rates))];
		xx <- match(x[1],last_splits);
		init_divides <- last_splits[last_splits!=x[1]];
		divides <- init_divides/nbins;
		### START HERE PETE!!!!	
		rate_runs <- cbind(c(1,init_divides),c(init_divides,nbins+1)-1)
#		init_divides <- c(init_divides,nbins+1);
		diverse_rates <- c();
		for (rr in 1:nrow(rate_runs))	{
#			rr <- rr+1;
			a <- rate_runs[rr,1];
			z <- rate_runs[rr,2];
			diverse_rates <- c(diverse_rates,sum(saturated_diversification_model[a:z,1]*interval_durations[a:z])/sum(interval_durations[a:z]));
#			print(c(divides[rr],diverse_rates[rr]));
			}
#		rate_runs <- cbind(rate_runs,diverse_rates);
		rate_info <- c(diverse_rates,divides);
		brl <- optim(rate_info,fn=log_likelihood_multiple_diversification_rates,method="L-BFGS-B",data_for_p_or_q=data_for_p_or_q,interval_durations=interval_durations,lower=rep(min_rate,ttl_rates),upper=rep(max_rate,ttl_rates),control=cl);
		
		k_rate_lnls <- c(k_rate_lnls,brl$value)
		if (brl$value > mxlgl)	{
			mbrl <- brl;
			mxlgl <- mbrl$value;
			no_improvement <- 0;
			} else	{
			no_improvement <- no_improvement+1;
			}
		}
	best_params <- reparifarge_rate_info(rate_info=mbrl$par,nbins);
	this_model <- best_params$diversification_rates;
	last_splits <- splits <- best_params$divides[best_params$divides %in% (2:nbins)];	# save only the splits that we vary

	diversification_models <- cbind(this_model,diversification_models);
	diversification_model_loglikelihoods <- c(mbrl$value,diversification_model_loglikelihoods);
	diversification_model_AICcs <- c(modified_AIC(lnL=mbrl$value,k=ttl_rates,n=ttl_data),diversification_model_AICcs);
#	examined_splits_org <- rbind(examined_splits_org,c(splits,rep(0,nbins-ttl_rates)));
	if (nbins>99)	{
		if (ttl_rates<10)	{
			m_n <- paste("k=00",ttl_rates,sep="");
			} else if (ttl_rates<100)	{
			m_n <- paste("k=0",ttl_rates,sep="");
			} else	{
			m_n <- paste("k=",ttl_rates,sep="");
			}
		} else if (nbins>9)	{
		if (ttl_rates<10)	{
			m_n <- paste("k=0",ttl_rates,sep="");
			} else	{
			m_n <- paste("k=",ttl_rates,sep="");
			}
		} else	{
		m_n <- paste("k=",ttl_rates,sep="");
		}
	model_names <- c(m_n,model_names);
	colnames(diversification_models) <- names(diversification_model_loglikelihoods) <- names(diversification_model_AICcs) <- model_names;
	print(c(ttl_rates,round(diversification_model_AICcs[length(diversification_model_AICcs):1],1)));
#	diversification_models;
	}

diverse_rate <- sum(saturated_diversification_model[,1]*interval_durations)/sum(interval_durations);
ttl_rates <- 1;
brl <- optim(diverse_rate,fn=log_likelihood_single_diversification_rate,method="L-BFGS-B",data_for_p_or_q=data_for_p_or_q,interval_durations=interval_durations,lower=min_rate,upper=max_rate,control=cl);
#brl <- optim(diverse_rate,fn=log_likelihood_single_diversification_rate,method="L-BFGS-B",data_for_p_or_q=data_for_p_or_q,lower=min_rate,upper=max_rate,control=cl);
diversification_models <- cbind(rep(brl$par,nrow(saturated_diversification_model)),diversification_models);
diversification_model_loglikelihoods <- c(brl$value,diversification_model_loglikelihoods);
diversification_model_AICcs <- c(modified_AIC(lnL=brl$value,k=ttl_rates,n=ttl_data),diversification_model_AICcs);
#rownames(diversification_models) <- rownames(data_for_p_or_q);
if (nbins>99)	{
	m_n <- paste("k=00",1,sep="")
	} else if (nbins>9)	{
	m_n <- paste("k=0",1,sep="")
	} else	{
	m_n <- paste("k=",1,sep="")
	}
model_names <- c(m_n,model_names);
colnames(diversification_models) <- names(diversification_model_loglikelihoods) <- names(diversification_model_AICcs) <- model_names;

output <- list(diversification_models,diversification_model_loglikelihoods,diversification_model_AICcs);
names(output) <- c("models","loglikelihoods","AICc");
return(output);
}

#binfixing_detritus <- function()	{
#while (length(unique(divides)) < length(divides) && attempts < 100)	{
#	fix_these <- c();
#	for (i in 1:length(divides))
#		if (sum(divides==divides[i])>1)	fix_these <- c(fix_these,i);
#	fix_these <- (1:sep_rates)[divides %in% unique(divides)];
#	mid_fix <- round(length(fix_these)/2);
#	alterations <- fix_these-fix_these[mid_fix];
#	while (min(divides[fix_these]+alterations)<1)	alterations <- alterations+1;
#	divides[fix_these] <- divides[fix_these]+alterations;
#	too_high <- fix_these[divides[fix_these]>=max_divide];
#	if (length(too_high)>0)	{
#		divides[too_high] <- divides[too_high]-(1+divides[max(too_high)]-max_divide);
#		}
#	divides <- sort(divides);
#	attempts <- attempts+1;
#	divides;
#	}
#}

#p_found_lineage_myr(mu,psi,chrons)	{
#p_duration <- 1;
#return(exp(-mu*chrons)*(1-exp(-psi*chrons)));
#}

p_samples_per_lineage_given_extinction_and_sampling <- function(mu,psi)	{
p_n <- 1-(psi/(mu+psi));	# prob. of missing a species
return(c(p_n,dpois((1:10),-1*log(p_n))));	# prob. of 1:10 finds given lambda & psi
#sum(individual_lineage_myr(q=lambda,tm=(1:1000)/1000))/1000
#integrate(individual_lineage_myr,lower=0,upper=1,q)
#dpois(0,psi/lambda)
#exp(-lambda);
#exp(-psi);
#integrate(individual_lineage_myr,lower=0,upper=2,q=mu)$value
}


