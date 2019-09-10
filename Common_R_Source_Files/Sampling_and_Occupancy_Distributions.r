##### UNITE AND/OR SEPARATE THIS FROM Distributions.r
#dev.off()
library(stringr)
library(gdata)
library(combinat)
library(gtools)
library(sads)
library(subplex)
#library(untb);
ZERO <- 5e-324
MINEXPN <- 10^-10
MINNO <- 5e-324

# routine to count collectins per sampling bin
count_collections_per_bin <- function(bin,locality_info)	{
return(sum(locality_info$Bin==bin))
}

# routine to count collectins per sampling bin for each taxon
accio_collections_per_taxon_per_bin <- function(taxon,occurrences,collections)	{
taxon_finds <- subset(occurrences,occurrences$Taxon==taxon)
taxon_localities <- collections[(collections$Locality %in% taxon_finds$Locality),]
b <- (1:max(collections$Bin))
return(sapply(b,count_collections_per_bin,locality_info=taxon_localities))
}

# first-order AIC
AIC <- function(lnL,k)	{
# L: log-likelihood; # k: parameters; # n: data points
return((-2*lnL)+(2*k))
}

# second-order (modified) AIC
modified_AIC <- function(lnL,k,n)	{
# L: log-likelihood; # k: parameters; # n: data points
#if (n==k)	n <- n+2
if (is.na(n / (n - k - 1)) || (n / (n - k - 1)<0) || (n / (n - k - 1))==Inf)	{
	aic_c <- AIC(lnL,k);
	}	else	aic_c <- (-2*lnL) + (2*k)*(n / (n - k - 1));
return(aic_c);
}

# convert abundances to "Fisher plot" giving # taxa with 1…N Finds
fisher_plot <- function(finds)	{
unique_finds <- sort(unique(finds),decreasing=FALSE)
observed <- vector(length=max(unique_finds))
for (s in 1:length(unique_finds))	observed[unique_finds[s]] <- length(finds[finds==unique_finds[s]])
return(observed)
}

# Get Chao2 Richness estimate from vector giving # taxa with 1…N finds
Chao2_Fisher <- function(observed)	{
ntaxa <- sum(observed)
S <- round(ntaxa + ((observed[1]*observed[1])/(2*(observed[2]+1))) - ((observed[1]*observed[2])/((2*(observed[2]+1)*(observed[2]+1)))))
return(S)
}

# Get Jackknife 2 Richness estimate from vector giving # taxa with 1…N finds
jack2_Fisher <- function(observed)	{
ntaxa <- sum(observed)
ss <- 0
for (i in 1:length(observed))	ss <- ss+(i*observed[i])
S <- round(ntaxa + (observed[1]*(((2*ss)-3)/ss)) - observed[2]*(((ss-2)*(ss-2))/(ss*(ss-1))))
return(S)
}

# Get Jackknife 5 Richness estimate from vector giving # taxa with 1…N finds
jack5_Fisher <- function(observed)	{
ntaxa <- sum(observed)
ss <- 0
for (i in 1:length(observed))	ss <- ss+(i*observed[i])
S <- round(ntaxa + (observed[1]*(((5*ss)-15)/ss)) - observed[2]*((10*(ss*ss)-(70*ss)+125)/(ss*(ss-1))) + observed[3]*(((10*(ss^3))-(120*(ss^2))+(485*ss)-660))/(ss*(ss-1)*(ss-2)) - observed[4]*((ss-4)^4)/(ss*(ss-1)*(ss-2)*(ss-3)) + observed[5]*((ss-5)^5)/(ss*(ss-1)*(ss-2)*(ss-3)*(ss-4)))
return(S)
}

# get expected occurrences given hypothesized occupancy distribution and numbers of collections
#	this generates an "expected" analog to fisher_plot above
expected_occurrences <- function(rocd, ncoll, S)	{
# find expectations of this gamma at this sample size
# exp will give the expected number of taxa found 1…ncoll times
#b <- seq(1,ncoll,by=1);
b <- 1:as.numeric(ncoll);
expected <- dbinom(b,ncoll,rocd[1])
for (j in 2:S)	{
	expected <- expected+dbinom(b,ncoll,rocd[j]);
	}
min <- 1.0e-323
for (i in 1:ncoll)	if (expected[i]==0)	expected[i] <- min
return(expected)
}

# get multinomial lognormal including P[observed total | expected total] as substitute for P[0]
distribution_loglikelihood_mul <- function(observed,expected,oS,hS)	{
#print(c(oS,hS))
mxfind <- length(observed)							# maximum finds observed
prop_expected <- expected[1:mxfind]/sum(expected)		# convert expected species to proportions
prop_expected[prop_expected==0] <- MINNO
lnlo <- observed*log(prop_expected)	# exact probability of observing # taxa with 1, 2, 3, etc. finds
lnlo[is.na(lnlo)] <- 0
#for (i in 1:mxfind)	if (is.na(lnlo[i]))	lnlo[i] <- 0
sobs <- sum(lnlo)
# log probability of observing X taxa from hypothesized hS taxa given expected # taxa with 1, 2, 3, etc. finds
eS <- sum(expected)									# get expected sampled species
if (eS==hS)	eS <- 0.99999999999*hS				# this is a kluge to get around rounding error of eS->hS
if (eS<hS)	{
	lnls <- lfactorial(hS)-(lfactorial(hS-oS)+lfactorial(oS))+(oS*log(eS/hS))+((hS-oS)*log((hS-eS)/hS))
	}	else if (oS<hS)	{
	lnls <- oS*log(MINNO)
	}	else	{
	lnls <- 0	# if we expect to observe all of the species and oS==hS, then P=1.0
	}
return(sobs+lnls)
}

# routine to get loglikelihood of a uniform rate of sampling among hS taxa given finds and possible finds
loglikelihood_uniform_occupancy_rates <- function(sc, hS, observed, oS, ncoll)	{
rocd <- rep(sc,hS)
#print(sc)
if (rocd[1]<=1)	{
	expected <- expected_occurrences(rocd,ncoll,hS)
		# log likelihood
	lnl <- distribution_loglikelihood_mul(observed,expected,oS,hS)
#	print(c(round(sc,20),round(lnl,3)))
	if (is.na(lnl))	lnl <- oS*log(MINNO)
	}	else {
	lnl <- oS*log(MINNO)
	}
return(lnl)
}

# routine to estimate most likely uniform rate of sampling/occupancy for hS taxa given finds & possible finds
optimize_uniform_occupancy_for_hS <- function(hS,observed,oS,ncoll)	{
#print(hS)	# for debugging
cl <- list(fnscale=-1)	# no idea what this means.....
sc <- (sum((1:length(observed))*observed)/hS)/ncoll;
minsc <- sc/100;
maxsc <- length(observed)/ncoll;
w <- optim(sc,fn=loglikelihood_uniform_occupancy_rates,method="L-BFGS-B",hS=hS,oS=oS,observed=observed,ncoll=ncoll,control=cl,lower=minsc,upper=maxsc);
bH <- c(w$par,hS,w$value)
names(bH) <- c("scale","richness","loglikelihood")
return(bH)
}

# routine to estimate most likely uniform rate of sampling/occupancy given finds & possible finds
optimize_uniform_occupancy <- function(finds,ncoll)	{
observed <- fisher_plot(finds)
oS <- length(finds[finds>0])	# observed taxa
minS <- stS <- length(finds)
pa <- 1
pz <- span <- 4		# span of richnesses to consider
enS <- stS+((span-1)*minS)	# highest richness to consider
incr <- floor(enS/span)
cl <- list(fnscale=-1)
peak <- 0
last_hS <- c(0,1,2,3)
while (incr>0)	{
	hS <- round(seq(stS,enS,by=incr),0)
	if (pz>length(hS))	pz <- length(hS)
	if (sum(hS %in% last_hS) == length(hS))	break
	results <- sapply(hS[pa:pz],optimize_uniform_occupancy_for_hS,observed=observed,oS=oS,ncoll=ncoll)
#	results <- sapply(hS[pa:pz],optimize_uniform_occupancy_for_hS,observed=observed,oS=oS,ncoll=ncoll)
	if (pa==2)
		results <-cbind(old_result_1,results)
	if (pz<length(hS))
		results <-cbind(results,old_result_2)
	mlnl <- max(results[3,])
	mlSc <- match(mlnl,results[3,])
	mlS <- hS[mlSc]
	bH <- results[,mlSc]
	spn <- length(results[3,])
	if (incr>1)	{
		# if runs so far are still producing higher likelihoods at higher richnesses:
		if (mlSc==spn)	{
			if (peak==0) {
				# if we have not yet found a peak, then keep increasing richness
				pa <- 2		# make note of best hypothesis
				old_result_1 <- results[,mlSc]	# set aside stats so as to avoid recalculating
				stS <- hS[spn]
				enS <- stS+((spn-1)*minS)
				}	else	{
				pa <- 2				# make note of 2nd best hypothesis
				pz <- (spn-1)		# make note of best hypothesis
				old_result_1 <- results[,spn-1]	# set aside stats so as to avoid recalculating
				old_result_2 <- results[,spn]	# set aside stats so as to avoid recalculating
				enS <- hS[spn]
				stS <- hS[spn-1]
				if (span<incr)	{
					incr <- (enS-stS)/(span-1)
					}	else {
					span <- (enS-stS)
					pz <- span-1
					incr <- 1
					}
#				stS <- hS[spn]-(span*incr)
				}	# end case where last number is best, but this is after finding a peak.
			} else if (mlSc==1)	{
			# if first richness is the best
			peak <- 1
			stS <- hS[mlSc]
			enS <- hS[mlSc+1]
			old_result_1 <- results[,1]
			pa <- 2
			old_result_2 <- results[,2]
			pz <- length(hS)-1
			if (incr==1)	incr <- 0
			if (span<incr)	{
				incr <- (enS-stS)/(span-1)
				}	else {
				span <- enS-stS
				incr <- 1
				}
			} else {
			# if a richness in the middle is best & we still can find a better one
			peak <- 1
			hS2 <- c(mlS-1,mlS+1)
			results2 <- sapply(hS2,optimize_uniform_occupancy_for_hS,observed=observed,oS=oS,ncoll=ncoll)
			if (results2[3,1]>mlnl && results2[3,1]>results2[3,2])	{
			# start just above 2nd best richness and go up to best
				stS <- hS[mlSc-1]
				enS <- hS2[1]
				old_result_1 <- results[,mlSc-1]
				old_result_2 <- results2[,1]
				if (span<incr)	{
					incr <- (enS-stS)/(span-1)
					}	else {
					span <- enS-stS
					incr <- 1
					}
				pz <- length(seq(stS,enS,by=incr))-1
				if (pz < pa)	incr <- 1
				}	else if (results2[3,2]>mlnl && results2[3,2]>results2[3,1]) {
			# start at best richness and go just below 2nd best
				stS <- hS2[2]
				enS <- hS[mlSc+1]
				old_result_1 <- results2[,2]
				old_result_2 <- results[,mlSc+1]
				if (incr>span)	{
					incr <- (enS-stS)/(span-1)
					}	else {
					span <- enS-stS
					incr <- 1
					}
				pz <- length(seq(stS,enS,by=incr))-1
				if (pz < pa)	incr <- 1
				}	else {
				# we already had the best, so just end it
				incr <- 0
				}
			}
		# end case where we have a better richness in middle somewhere.  
		}	else	{
		incr <- 0
		}
	last_hS <- hS
	}
bH[3]<-round(bH[3],2)
uniform_AICc <- round(modified_AIC(bH[3],2,sum(finds)),2)
names(uniform_AICc) <- "AICc"
bH <- c(bH,uniform_AICc)
return(bH)
}

accio_richness_with_prop_greater_than_x_given_exponential <- function(decay,x)	{
d <- min(decay,1/decay)
return(round(1+((log(x)-log(1-d))/log(d))))
}

# exponential distribution with decay rate 1/ev; note that this is truncated at some minimum number
# generate geometric (= exponential) distribution
exponential_distribution <- function(decay)	{
#print(decay);
d <- min(decay,1/decay)
S <- accio_richness_with_prop_greater_than_x_given_exponential(decay,MINEXPN)
ranks <- (1:S)
prop <- (1-d)*(d^(ranks-1))
return(prop)
}

# exponential distribution with decay rate that is rescaled to some "average" rate; note that this is truncated at some minimum number
scaled_exponential_distribution <- function(sc,decay)	{
return(sc*exponential_distribution(decay))
}

# find the decays from one rank to the next
accio_prop_drops <- function(finds)	{
oS <- length(finds)
s1 <- (1:(oS-1))
s2 <- (2:oS)
return(finds[s1]/finds[s2])
}

# upper and lower bounds for plausible exponentials
accio_param_bounds_exponential <- function(finds,observed,ncoll)	{
# step 1: get numbers of species with 1, 2, 3, etc. finds
prop_drops <- unique(accio_prop_drops(finds))
oS <- length(finds)
decay_min <- 1+((min(subset(prop_drops,prop_drops>1))-1)/5)
decay_max <- min(max(prop_drops),1/exp(log(MINEXPN)/oS))
# get the exponential distributions from the two extreme slopes
dmin <- exponential_distribution(decay_min)
dmax <- exponential_distribution(decay_max)
# set scaling so that the most common taxon will be within one half to twice the prop
#	of the observed
scale_min <- (finds[1]/ncoll)/dmax[1]
scale_max <- (finds[1]/ncoll)/dmin[1]
return(c(scale_min,scale_max,decay_min,decay_max))
}

# rough initial estimate of an exponential distribution to seed searches
accio_initial_occupancy_exponential <- function(finds,observed,ncoll)	{
p0 <- vector(length=2)	# pO[1]=rate adjuster; p0[2]=decay
names(p0) <- c("scale", "decay")
# step 1: estimate average shifts
prop_drops <- accio_prop_drops(finds[finds>0]);
p0[2] <- mean(prop_drops)
d <- exponential_distribution(p0[2])
p0[1] <- (finds[1]/ncoll)/d[1]
names(p0) <- c("scale","decay")
return(p0)
}

# get loglikelihood of exponential distribution of sampling rates given finds and collections
loglikelihood_exponential_occupancy_rates <- function(p0, observed, oS, ncoll)	{
r <- p0[1]
decay <- p0[2]
#print(p0)
# p0[1]: r	# p0[2]: ev	# p0[3]: S
rocd <- scaled_exponential_distribution(r,decay)		# basic exponential distribution
if (length(rocd)<oS)	{
	hS <- Chao2_Fisher(observed)
	dummy <- vector(length=(hS-length(rocd)))
	dummy[1] <- rocd[length(rocd)]/decay
	if (length(dummy)>1)		for (i in 2:length(dummy))	dummy[i] <- dummy[i-1]/decay
	rocd <- c(rocd,dummy)
#	rocd <- scaled_exponential_distribution_min_rich(r,decay,hS)
	}	else	{
	minp <- ((10^-7)/ncoll)
	hS <- min(Chao2_Fisher(observed),length(subset(rocd,rocd>minp)))
	}
#rocd <- rocd*max(finds)/ncoll
if (rocd[1]<=1 && hS>=oS)	{
	expected <- expected_occurrences(rocd,ncoll,hS)[1:length(observed)]
		# log likelihood
	lnl <- distribution_loglikelihood_mul(observed,expected,oS,hS)
	}	else {
	lnl <- oS*log(MINNO)
	}
#print(c(r,decay,lnl))
return(round(lnl,2))
}

# optimize exponential distribution
optimize_exponential_occupancy <- function(finds,ncoll)	{
observed <- fisher_plot(finds);
oS <- length(finds);	# observed taxa
finds <- finds[finds>0];
#prop_finds <- finds/ncoll
p0 <- accio_initial_occupancy_exponential(finds,observed,ncoll);	# pO[1]=; p0[3]=richness
bnds <- accio_param_bounds_exponential(finds,observed,ncoll);
cl <- list(fnscale=-1)
w <- optim(p0,fn=loglikelihood_exponential_occupancy_rates,method="L-BFGS-B",oS=oS,observed=observed,ncoll=ncoll,lower=c(bnds[1],bnds[3]),upper=c(bnds[2],bnds[4]),control=cl)

d <- 1/w$par[2]
hS <- 1+round((log((10^-7)/(1-d)))/log(d),0)
bH <- c(w$par,hS,w$value,round(modified_AIC(w$value,2,sum(finds)),2))
names(bH) <- c("scale","decay","richness","loglikelihood","AICc")
return(bH)
}

# generate lognormal distribuiton for S entities with log-variance = ev
lognormal_distribution <- function(mag, S)	{
fi <- seq(1/(S+1),S/(S+1),by=1/(S+1))
prop <- mag^(qnorm(fi,0,1))/sum(mag^(qnorm(fi,0,1)))
return(prop)
}

# generate lognormal distribuiton for S entities with log-variance = mag and rescaled by r
scaled_lognormal_distribution <- function(sc, mag, S)	{
return(sort(sc*lognormal_distribution(mag, S),decreasing=TRUE))
}

# get min-max scale for lognormal
accio_min_max_lognormal_scale <- function(finds,observed,oS,ncoll,p0)	{
d <- lognormal_distribution(p0[2],p0[3])
sd <- finds/(ncoll*d[1:oS])
return(c(min(sd),max(sd)))
}

# get min-max logvariance for lognormal
accio_min_max_lognormal_mag <- function(finds,observed,oS,ncoll)	{
smax <- max(Chao2_Fisher(observed),jack2_Fisher(observed))
### get minimum and maximum parameters
maxmidfnds <- finds[round(oS/2)]
minmidfnds <- finds[oS]/finds[round(oS/2)]
	
benchers <- round(oS/10)
# for max, get the difference between minimum & benchers using oS;
# for min, get the difference between the max & benchers using hS
mnev <- mxev <- 0
for (i in 1:benchers)	{
	s1 <- smax-(i-1)
	s2 <- oS-(i-1)
	x1 <- qnorm(s1/(s1+1),0,1)			# get q-value associated with dominant taxon
	x2 <- qnorm(s2/(s2+1),0,1)			# get q-value associated with dominant taxon
	mnev <- mnev+exp(log(finds[i]/maxmidfnds)/x1)/benchers
	mxev <- mxev+exp(log(finds[i]/minmidfnds)/x2)/benchers
	}	# get rough estimate of magnitude parameter based on Top 10 taxa
mnev <- 1+(mnev-1)/2
mxev <- 1+2*(mxev-1)
return(c(mnev,mxev))
}

# get rough estimate of lognormal given variance in log occupancy 
# get rough estimate of lognormal given variance in log occupancy 
accio_initial_occupancy_lognormal <- function(finds,observed,oS,minS,ncoll)	{
# editted 2019-08-15
p0 <- vector(length=3)	# pO[1]=; p0[3]=richness
names(p0) <- c("scale", "magnitude", "richness")
	# step 1: get rough estimate of lognormal based on variance above the median
p0[3] <- min(Chao2_Fisher(observed),jack2_Fisher(observed))	# lowest richness estimator we will consider
if (p0[3]<=minS)
	p0[3] <- minS+(p0[3]-oS);

if (p0[3]<=oS)	{
	mid <- round(p0[3]/2);
	for (i in 1:min(p0[3],10))	{			# do not let number exceed richness!
		s <- p0[3]-(i-1)
		x1 <- qnorm(s/(s+1),0,1)			# get q-value associated with dominant taxon
		p0[2] <- p0[2]+exp(log(finds[i]/finds[mid])/x1)/10
		}	# get rough estimate of magnitude parameter based on Top 10 taxa
	d <- lognormal_distribution(p0[2],p0[3])
	p0[1] <- (finds[mid]/ncoll)/d[mid]	# set median rate to equal rate of extrapolated median
	}	else	{
	dc <- length(subset(observed,observed>0))
	k <- 1
	uniqfnds <- sort(unique(finds),decreasing=TRUE)
	uo <- observed[observed>0]
	uniqobs <- uo[dc:1]
	midpts <- vector(length=dc)
	midpts[1] <- (1+uniqobs[1])/2
	for (k in 2:dc)	{
		midpts[k] <- midpts[k-1]+((1+uniqobs[k])/2)
		}
	quants <- qnorm((p0[3]-(midpts-1))/(p0[3]+1))	# mag^quant[i] is how many times more common a species is than the median
	### each quant[i] gives the quantile of the median ranked taxon with uniqfnds[i] specimens
	for (s1 in 1:(dc-1))	{
		f1 <- uniqfnds[s1]
		for (s2 in (s1+1):dc)	{
			f2 <- uniqfnds[s2]
			if (s1==1 && s2==2)	{
				avem <- (f1/f2)^(1/(quants[s1]-quants[s2]))
				}	else	{
				avem <- c(avem,(f1/f2)^(1/(quants[s1]-quants[s2])))
				}
			}
		}
	# use the minimum of the median, arithmetic mean and geometric mean.
	p0[2] <- min(median(avem),mean(avem),(exp(sum(log(avem))))^(1/length(avem)))
	d <- sort(lognormal_distribution(p0[2],p0[3]),decreasing=TRUE)
	p0[1] <- median((finds[1:oS]/ncoll)/d[1:oS]);
	if ((p0[1]*d[1])>1)	p0[1] <- 0.99/d[1]
	}
return(p0)
}

# get lower and upper bounds for lognormal magnitude and rescaling
accio_min_max_lognormal_occupancy_params <- function(finds,observed,oS,ncoll)	{
mmm <- finds[1]/finds[oS]
# use a maximum richness estimate to get maximum magnitude
#	this works because the Nth taxon is a higher overall rank given high richness
#	we therefore go from obs. min -> obs. max over the smallest shift in rel. ranks
smax <- max(Chao2_Fisher(observed),jack2_Fisher(observed))
smb <- 1-(oS/(smax+1))
s1b <- 1-(1/(smax+1))
qmb <- qnorm(smb,0,1)
q1b <- qnorm(s1b,0,1)

mag_max <- 2*exp(log(mmm)/(q1b-qmb))

### use observed richness to get lowest estimate.
smin <- oS
sma <- 1-(oS/(smin+1))
s1a <- 1-(1/(smin+1))
qma <- qnorm(sma,0,1)
q1a <- qnorm(s1a,0,1)
mag_min <- exp(log(mmm)/(q1a-qma))/2

# use fewest taxa and biggest magnitude shift to get the most common possible #1
ru <- sort(lognormal_distribution(mag_max,smin),decreasing=TRUE)
# use most taxa and smallest magnitude shift to get the least common possible #1
rl <- sort(lognormal_distribution(mag_min,smax),decreasing=TRUE)

sc_raw <- finds[1]/ncoll
sc_min <- sc_raw/ru[1]
sc_max <- sc_raw/rl[1]
#ru2 <- sort(lognormal_distribution(mag_max,smax),decreasing=TRUE)
#print(c(rl1[1],rl2[1],ru1[1],ru2[1]))
boundary_params <- c(mag_min,mag_max,sc_min,sc_max)
names(boundary_params) <- c("magn_min","magn_max","scale_min","scale_max")
return(boundary_params)
}

# calculate loglikelihood of lognormal for occupancy
loglikelihood_lognormal_occupancy_rates <- function(p0, hS, observed, oS, ncoll)	{
# p0[1]: r	# p0[2]: ev	# p0[3]: S		print(p0)
r <- p0[1]
mag <- p0[2]
rocd <- sort(scaled_lognormal_distribution(r,mag,hS),decreasing=TRUE)		# basic lognormal distribution
#print(c(r,mag))
if (rocd[1]<=1 && length(rocd[rocd>0])==length(rocd))	{
	# expected number of taxa with 1, 2, 3, etc. finds
	#	expected_prop <- expected_prop(rocd,ncoll,p0[3])
	expected <- expected_occurrences(rocd,ncoll,hS)
	# log likelihood
	lnl <- distribution_loglikelihood_mul(observed,expected,oS,hS)
	}	else {
	lnl <- 10000*-oS
	}
return(lnl)
}

# find best lognormal occupancy distribution for particular richness
optimize_lognormal_occupancy_given_hS <- function(hS,oS,observed,ncoll,p0,sc,mm,debug=0)	{
if (debug==1)	print(hS)
#print(hS)
cl <- list(fnscale=-1)
w <- optim(c(br=p0[1],mag=p0[2]),fn=loglikelihood_lognormal_occupancy_rates,method="L-BFGS-B",oS=oS,hS=hS,observed=observed,ncoll=ncoll,lower=c(min(sc),min(mm)),upper=c(max(sc),max(mm)),control=cl)
bH <- c(w$par,hS,w$value)
names(bH) <- c("scale","mag_var","richness","loglikelihood")
return(bH)
}

# find best lognormal occupancy distribution
optimize_lognormal_occupancy <- function(finds,ncoll)	{
observed <- fisher_plot(finds);
oS <- length(finds[finds>0])					# observed + inferred taxa
maxS <- max(oS,(2*Chao2_Fisher(observed)));		# observed + inferred taxa
minS <- length(finds)							# observed + inferred taxa
#observed <- fisher_plot(finds)
p0 <- accio_initial_occupancy_lognormal(finds,observed,oS,minS,ncoll)	# pO[1]=r; p0[1]=magnitude; p0[3]=richness
b0 <- accio_min_max_lognormal_occupancy_params (finds,observed,oS,ncoll)
mm <- b0[1:2];
sc <- b0[3:4];
# make sure that lower & upper bounds contain starting values!
if (mm[1] > p0[2])	{
	mm[1] <- 1+((p0[2]-1)/2)
	}	else if (mm[2] < p0[2])	{
	mm[2] <- 2*p0[2]
	}
if (sc[1] > p0[1])	{
	sc[1] <- p0[1]/2
	}	else if (sc[2] < p0[1])	{
	sc[2] <- p0[1]
	}
pa <- 1
pz <- span <- 4;		# span of richnesses to consider
stS <- minS;
enS <- stS+((span-1)*minS);
incr <- floor(enS/span);
cl <- list(fnscale=-1)
peak <- 0;
last_hS <- c(0,1,2,3);
while (incr>0)	{
	hS <- round(seq(stS,enS,by=incr),0);
	if (max(hS)>maxS)	{
		incr <- (maxS-min(hS))/(length(hS)-1);
		hS <- round(seq(min(hS),maxS,by=incr),0);
		}
	if (sum(hS %in% last_hS) == length(hS))	break;
	if (pz>length(hS))	pz <- length(hS);
	results <- sapply(hS[pa:pz],optimize_lognormal_occupancy_given_hS,oS=oS,observed=observed,ncoll=ncoll,p0=p0,sc=sc,mm=mm,debug=0)
	if (pa==2)
		results <-cbind(old_result_1,results)
	if (pz<length(hS))
		results <-cbind(results,old_result_2)
	mlnl <- max(results[4,])
	mlSc <- match(mlnl,results[4,])
	mlS <- hS[mlSc]
	bH <- results[,mlSc]
	spn <- length(results[1,])
	if (incr>1)	{
		# if runs so far are still producing higher likelihoods at higher richnesses:
		if (mlSc==spn)	{
			if (peak==0) {
				# if we have not yet found a peak, then keep increasing richness
				pa <- 2		# make note of best hypothesis
				old_result_1 <- results[,mlSc]	# set aside stats so as to avoid recalculating
				stS <- hS[spn]
				enS <- stS+((spn-1)*minS)
				}	else	{
				stS <- hS[mlSc-1]
				enS <- hS[mlSc]
				old_result_1 <- results[,mlSc-1]	# set aside stats so as to avoid recalculating
				old_result_2 <- results[,mlSc]	# set aside stats so as to avoid recalculating
				if (span<incr)	{
					incr <- (enS-stS)/(span-1)
					}	else {
					span <- enS-stS
					incr <- 1
					}
				pa <- 2
				pz <- length(seq(stS,enS,by=incr))-1
#				stS <- hS[spn]-(span*incr)
				}	# end case where last number is best, but this is after finding a peak.
			} else if (mlSc==1)	{
			# if first richness is the best
			peak <- 1
			stS <- hS[mlSc]
			enS <- hS[mlSc+1]
			old_result_1 <- results[,1]
			pa <- 2
			old_result_2 <- results[,2]
			pz <- length(seq(stS,enS,by=incr))-1
			if (incr==1)	incr <- 0
			if (span<incr)	{
				incr <- (enS-stS)/(span-1)
				}	else {
				span <- enS-stS
				incr <- 1
				}
			} else {
			# if a richness in the middle is best & we still can find a better one
			peak <- 1
			hS2 <- c(mlS-1,mlS+1)
			results2 <- sapply(hS2,optimize_lognormal_occupancy_given_hS,oS=oS,observed=observed,ncoll=ncoll,p0=p0,sc=sc,mm=mm,debug=0)
			if (results2[4,1]>mlnl && results2[4,1]>results2[4,2])	{
			# start just above 2nd best richness and go up to best
				stS <- hS[mlSc-1]
				enS <- hS2[1]
				old_result_1 <- results[,mlSc-1]
				old_result_2 <- results2[,1]
				if (span<incr)	{
					incr <- (enS-stS)/(span-1)
					}	else {
					span <- enS-stS
					incr <- 1
					}
				pa <- 2
				pz <- length(seq(stS,enS,by=incr))-1
				if (pz < pa)	incr <- 1
				}	else if (results2[4,2]>mlnl && results2[4,2]>results2[4,1]) {
			# start at best richness and go just below 2nd best
				stS <- hS2[2]
				enS <- hS[mlSc+1]
				old_result_1 <- results2[,2]
				old_result_2 <- results[,mlSc+1]
				if (incr>span)	{
					incr <- (enS-stS)/(span-1)
					}	else {
					span <- enS-stS
					incr <- 1
					}
				pa <- 2
				pz <- length(seq(stS,enS,by=incr))-1
				if (pz < pa)	incr <- 1
				}	else {
				# we already had the best, so just end it
				incr <- 0
				}
			pa <- 2
			pz <- spn-1
			}
		# end case where we have a better richness in middle somewhere.  
		last_hS <- hS
		}	else	{
		incr <- 0
		}
#	print(hS)		# for debugging
#	if (incr>0)	print(round(seq(stS,enS,by=incr),0))	# for debugging
	}
lognormal_AICc <- round(modified_AIC(mlnl,3,sum(finds)),2)
bH <- c(bH[1],bH[2],bH[3],round(bH[4],2),lognormal_AICc)
names(bH)[5] <- "AICc"
return(bH)
}

# get beta distribution
beta_distribution <- function(shape1, shape2, S)	{
# NOTE:I get some wonky assed results using dbeta & qbeta
ranks <- (1:S)/(S+1)
numer <- (ranks^(shape1-1))*((1-ranks)^(shape2-1))
#numer[numer==0] <- ZERO
# use log gammas, as gamma does not like numbers much past 100
Bab <- exp((lgamma(shape1)+lgamma(shape2))-lgamma(shape1+shape2))
cdf <- numer/Bab
pdf <- cdf/sum(cdf)
return(pdf)
}

accio_initial_occupancy_beta <- function(obs_mean, obs_var)	{
shape1 <- 0.1
shape2 <- (shape1/obs_mean)-shape1
combo_var <- (shape1*shape2)/(((shape1+shape2)^2)*(shape1+shape2+1))

if (combo_var>obs_var)	{
	while (combo_var>obs_var)	{
		shape1 <- 1.01*shape1
		shape2 <- (shape1/obs_mean)-shape1
		combo_var <- (shape1*shape2)/(((shape1+shape2)^2)*(shape1+shape2+1))
		}
	}	else if (combo_var<obs_var)	{
	while (combo_var<obs_var)	{
		shape1 <- shape1/1.01
		shape2 <- (shape1/obs_mean)-shape1
		combo_var <- (shape1*shape2)/(((shape1+shape2)^2)*(shape1+shape2+1))
		}
	}
return(c(shape1,shape2))
}

# calculate loglikelihood of beta for occupancy
loglikelihood_beta_occupancy_rates <- function(p0, hS, observed, oS, ncoll)	{
# p0[1]: shape1	# p0[2]: beta	
shape1 <- p0[1]
shape2 <- p0[2]
#print(c(shape1,shape2,hS,oS))
mxobs <- length(observed)
if (shape2>=shape1)	{
	rocd <- beta_distribution(shape1, shape2, hS)
	#print(c(r,mag))
	if (rocd[1]<=1 && length(rocd[rocd>0])==length(rocd))	{
	# expected number of taxa with 1, 2, 3, etc. finds
	#	expected_prop <- expected_prop(rocd,ncoll,p0[3])
		expected <- expected_occurrences(rocd,ncoll,hS)
		# log likelihood
		lnl <- distribution_loglikelihood_mul(observed,expected,oS,hS)
		}	else {
		lnl <- 10000*-oS
		}
	}	else lnl <- 10000*-oS
return(round(lnl,3))
}

# find best beta occupancy distribution
optimize_beta_occupancy_given_hS <- function(hS,oS,finds,observed,ncoll)	{
#print(hS)
obs_mean <- mean(c(finds,rep(0,hS-oS)))/ncoll
obs_var <- var(c(finds/ncoll,rep(0,hS-oS)))
p0 <- accio_initial_occupancy_beta(obs_mean,obs_var)
names(p0) <- c("beta_alpha","beta_beta")
cl <- list(fnscale=-1)
w <- optim(p0,fn=loglikelihood_beta_occupancy_rates,method="L-BFGS-B",oS=oS,hS=hS,observed=observed,ncoll=ncoll,lower=p0/100, upper=100*p0,control=cl)
bH <- c(w$par,hS,w$value);
names(bH) <- c("alpha", "beta", "richness","loglikelihood");
return(bH);
}

# find best beta occupancy distribution
optimize_beta_occupancy <- function(finds,ncoll)	{
observed <- fisher_plot(finds)
oS <- length(finds[finds>0])		# observed taxa
minS <- length(finds)				# observed + gap taxa
#bS <- min(Chao2_Fisher(observed),jack2_Fisher(observed))	# lowest richness estimator we will consider
pa <- 1
pz <- span <- 4
stS <- minS
enS <- stS+((span-1)*minS)
incr <- floor(enS/span)
cl <- list(fnscale=-1)
peak <- 0
last_hS <- c(0,1,2,3)
while (incr>0)	{
	hS <- round(seq(stS,enS,by=incr),0)
	if (pz>length(hS))	pz <- length(hS)
	if (sum(hS %in% last_hS) == length(hS))	break
	results <- sapply(hS[pa:pz],optimize_beta_occupancy_given_hS,oS=oS,finds=finds,observed=observed,ncoll=ncoll)
	if (pa==2)
		results <-cbind(old_result_1,results)
	if (pz<length(hS))
		results <-cbind(results,old_result_2)
	mlnl <- max(results[4,])
	mlSc <- match(mlnl,results[4,])
	mlS <- hS[mlSc]
	bH <- results[,mlSc]
	spn <- length(results[1,])
	if (incr>1)	{
		# if runs so far are still producing higher likelihoods at higher richnesses:
		if (mlSc==spn)	{
			if (peak==0) {
				# if we have not yet found a peak, then keep increasing richness
				pa <- 2		# make note of best hypothesis
				old_result_1 <- results[,mlSc]	# set aside stats so as to avoid recalculating
				stS <- hS[spn]
				enS <- stS+((spn-1)*minS)
				}	else	{
				stS <- hS[mlSc-1]
				enS <- hS[mlSc]
				old_result_1 <- results[,mlSc-1]	# set aside stats so as to avoid recalculating
				old_result_2 <- results[,mlSc]	# set aside stats so as to avoid recalculating
				if (span<incr)	{
					incr <- (enS-stS)/(span-1)
					}	else {
					span <- enS-stS
					incr <- 1
					}
				pa <- 2
				pz <- length(seq(stS,enS,by=incr))-1
#				stS <- hS[spn]-(span*incr)
				}	# end case where last number is best, but this is after finding a peak.
			} else if (mlSc==1)	{
			# if first richness is the best
			peak <- 1
			stS <- hS[mlSc]
			enS <- hS[mlSc+1]
			old_result_1 <- results[,1]
			pa <- 2
			old_result_2 <- results[,2]
			pz <- length(hS)-1
			if (incr==1)	incr <- 0
			if (span<incr)	{
				incr <- (enS-stS)/(span-1)
				}	else {
				span <- enS-stS
				incr <- 1
				}
			} else {
			# if a richness in the middle is best & we still can find a better one
			peak <- 1
			hS2 <- c(mlS-1,mlS+1)
			results2 <- sapply(hS2,optimize_beta_occupancy_given_hS,oS=oS,finds=finds,observed=observed,ncoll=ncoll)
			if (results2[4,1]>mlnl && results2[4,1]>results2[4,2])	{
			# start just above 2nd best richness and go up to best
				stS <- hS[mlSc-1]
				enS <- hS2[1]
				old_result_1 <- results[,mlSc-1]
				old_result_2 <- results2[,1]
				if (span<incr)	{
					incr <- (enS-stS)/(span-1)
					}	else {
					span <- enS-stS
					incr <- 1
					}
				pa <- 2
				pz <- length(seq(stS,enS,by=incr))-1
				if (pz < pa)	incr <- 1
				}	else if (results2[4,2]>mlnl && results2[4,2]>results2[4,1]) {
			# start at best richness and go just below 2nd best
				stS <- hS2[2]
				enS <- hS[mlSc+1]
				old_result_1 <- results2[,2]
				old_result_2 <- results[,mlSc+1]
				if (incr>span)	{
					incr <- (enS-stS)/(span-1)
					}	else {
					span <- enS-stS
					incr <- 1
					}
				pa <- 2
				pz <- length(seq(stS,enS,by=incr))-1
				if (pz < pa)	incr <- 1
				}	else {
				# we already had the best, so just end it
				incr <- 0
				}
			}
		# end case where we have a better richness in middle somewhere.  
		}	else	{
		incr <- 0
		}
#	print(hS)		# for debugging
#	if (incr>0)	print(round(seq(stS,enS,by=incr),0))	# for debugging
	last_hS <- hS
	}
beta_AICc <- round(modified_AIC(mlnl,3,sum(finds)),2)
bH <- c(bH[1],bH[2],bH[3],round(bH[4],2),beta_AICc)
names(bH)[5] <- "AICc"
return(bH)
}

# calculate the likelihood and AICc of a model predicting exactly the numbers of taxa with 1…N finds
#	Note: this ad hoc model has one parameter for each taxon
saturated_occupancy_model <- function(finds)	{
observed <- fisher_plot(finds)
oS <- sum(observed)
prop_expected <- observed/sum(observed)
sat_lnl <- sum(observed[observed>0]*log(prop_expected[prop_expected>0]))
sat_AICc <- modified_AIC(sat_lnl,oS,sum(finds))
bH <- c(round(sat_lnl,2),round(sat_AICc,2))
names(bH) <- c("Saturated_log-likelihood","Saturated_AICc")
return(bH)
}

# function to actually fit occurrence distributions in a particular time interval
test_occurrence_distributions <- function(stg,bins,taxon_finds_per_bin,collections_per_bin,tests,bin_labels)	{
# stg: the interval ("stage")
# bins: a vector giving the numbers of the bins, to find the actual number
#	corresponding to stg
bn <- bins[stg]
finds <- sort(taxon_finds_per_bin[bn,],decreasing=TRUE)[sort(taxon_finds_per_bin[bn,],decreasing=TRUE)>0]
#finds <- sort(taxon_finds_per_bin[bn,taxon_finds_per_bin[bn,]>0],decreasing=TRUE)
ncoll <- collections_per_bin[stg]
output_summary_bin <- c(bn,as.character(bin_labels[bn]),length(finds),sum(finds),ncoll)
names(output_summary_bin) <- c("Bin_No","Bin","Taxa","Occurrences","Collections")
ct <- 0
if (!is.na(match("Uniform",tests)))	{
	update <- paste("Doing best Uniform distribution for",as.character(bin_labels[bn]),date(),sep=" ")
	print(update)
	bst_uni <- optimize_uniform_occupancy(finds,ncoll)
	output_summary_bin <- c(output_summary_bin,bst_uni)
	ct <- ct+1
	aiccs <- bst_uni[4]
	if (bst_uni[4]==min(aiccs))	winner <- "Uniform"
	}
if (!is.na(match("Exponential",tests)))	{
	update <- paste("Doing best Exponential distribution for",as.character(bin_labels[bn]),date(),sep=" ")
	print(update)
	bst_exp <- optimize_exponential_occupancy(finds,ncoll)
	output_summary_bin <- c(output_summary_bin,bst_exp)
	ct <- ct+1
	if (ct==1)	{
		aiccs <- bst_exp[5]
		}	else {
		aiccs <- c(aiccs,bst_exp[5])
		}
	if (bst_exp[5]==min(aiccs))	winner <- "Exponential"
	}
if (!is.na(match("Lognormal",tests)))	{
	update <- paste("Doing best Lognormal distribution for",as.character(bin_labels[bn]),date(),sep=" ")
	print(update)
	bst_lgn <- optimize_lognormal_occupancy(finds,ncoll)
	output_summary_bin <- c(output_summary_bin,bst_lgn)
	ct <- ct+1
	if (ct==1)	{
		aiccs <- bst_lgn[5]
		}	else {
		aiccs <- c(aiccs,bst_lgn[5])
		}
	if (bst_lgn[5]==min(aiccs))	winner <- "LogNormal"
	}
if (!is.na(match("Beta",tests)))	{
	update <- paste("Doing best Beta distribution for",as.character(bin_labels[bn]),date(),sep=" ")
	print(update)
	bst_bta <- optimize_beta_occupancy(finds,ncoll)
	output_summary_bin <- c(output_summary_bin,bst_bta)
	ct <- ct+1
	if (ct==1)	{
		aiccs <- bst_bta[5]
		}	else {
		aiccs <- c(aiccs,bst_bta[5])
		}
	if (bst_bta[5]==min(aiccs))	winner <- "Beta"
	}
# get best possible likelihood (i.e., predicts exact distribution with N parameters = No. of taxa)
bst_sat <- saturated_occupancy_model(finds)
ct <- ct+1
if (ct==1)	{
	aiccs <- bst_sat[2]
	}	else {
	aiccs <- c(aiccs,bst_sat[2])
	}
if (bst_sat[2]==min(aiccs))	winner <- "Saturated"
output_summary_bin <- c(output_summary_bin,bst_sat)
names(winner) <- "Best_Model"
output_summary_bin <- c(output_summary_bin,winner)
### add something to state which hypothesis is best
return(output_summary_bin)
}

# Main program using functions above.
# 	occurrences_file specifies tab delimited text file giving collection # and taxon #
# 	collections_file specifies tab delimited text file giving collection # and sampling bin (e.g., stratigraphic stages, geographic units, geographic units by stages, etc.)
# 	bin_labels_file specifies single column text file giving names of sampling units
# 	tests specify distributions to test
# 	min_coll specifies minimum number of collections to test
# 	headers: if the occurrences and collections files have headers such as "Locality" or "Taxon", then TRUE
#		NOTE: these are renamed in the function.
estimate_occupancy_distributions <- function(occurrences_file="Occurrence_Information.txt",collections_file="Locality_Information.txt",study_name="",bin_labels_file="Bin_Labels.txt",tests=c("Uniform","Exponential","Lognormal","Beta"),min_coll=15,headers=TRUE)	{
occurrences <- read.table(file=occurrences_file, header=headers, stringsAsFactors=FALSE, sep="\t")
colnames(occurrences) <- c("Locality","Taxon")
# get taxon data
taxa <- sort(unique(occurrences$Taxon))
# get locality information
collections <- read.table(file=collections_file, header=headers, stringsAsFactors=FALSE, sep="\t")
colnames(collections) <- c("Locality","Bin")
locales <- length(collections$Locality)		# number of localities
bins <- sort(unique(collections$Bin))			# vector giving stratigraphic/geographic/etc. bins
ttl_bins <- length(bins)						# number of bins
# routine to count collections/localities per bin
collections_per_bin <- sapply(bins,count_collections_per_bin,locality_info=collections)

# this routine works on some computers but not others: I can find no rhyme or reason to this!
if (!is.na(bin_labels_file)  && bin_labels_file!="")	{
	bin_labels <- simplify2array(read.table(file=bin_labels_file, header=FALSE, stringsAsFactors=TRUE, sep="\t"))
	names(bins) <- names(collections_per_bin) <- bin_labels[bins]
	}	else {
	bin_labels <- as.character(bins)
	}
#bin_labels <- simplify2array(read.table(file=bin_labels_file, header=FALSE, stringsAsFactors=TRUE, sep="\t"))
names(collections_per_bin) <- bin_labels[bins]

# get finds per taxon per sampling bin
taxon_finds_per_bin <- sapply(taxa,accio_collections_per_taxon_per_bin,occurrences=occurrences,collections=collections)
rownames(taxon_finds_per_bin) <- bin_labels

stg <- (1:ttl_bins)[collections_per_bin>=min_coll]
output_summary <- sapply(stg,test_occurrence_distributions,bins,taxon_finds_per_bin,collections_per_bin,tests,bin_labels)
colnames(output_summary) <- names(collections_per_bin[collections_per_bin>=min_coll])
if (study_name=="")	{
	output_file <- "Occupancy_Distribution_Output.xls"
	}	else	{
	output_file <- paste(study_name,"_Occupancy_Distribution_Output.xls",sep="")	
	}
write.table(t(output_summary),output_file,sep="\t",col.names=TRUE,row.names=FALSE)
return(output_summary)
}

prob_missing_given_lognormal <- function(p,median_sampling,stdev_mag,ncoll)	{
#r <- median_sampling*stdev_mag^qnorm(p)
#print(p)		for debugging
if (median_sampling*stdev_mag^qnorm(p)>1)	{
	return(0)
	}	else	{
	return((1-median_sampling*stdev_mag^qnorm(p))^ncoll)
	}
#return((1-(median_sampling*stdev_mag^qnorm(p)))^ncoll)
}

#rescale <- sampling_lgn_coll[st,1]
#stdev_mag <- sampling_lgn_coll[st,2]
#S <- sampling_lgn_coll[st,3]
rieman_prob_missing_given_lognormal <- function(rescale,stdev_mag,S,ncoll,min_rho=0)	{
lgnd <- scaled_lognormal_distribution(rescale,stdev_mag,S);
#lgnd <- lgnd[lgnd>min_rho];
Sg <- sum(lgnd>min_rho);
median_sampling <- median(lgnd)
p <- seq((1/(S+1)),(S/(S+1)),by=(1/(S+1)))[(S-Sg+1):S];
xx <- sapply(p,prob_missing_given_lognormal,median_sampling,stdev_mag,ncoll)
return(mean(xx))
#return(mean((1-lgnd)^ncoll))
#for (i in 1:S)	prob_missing_given_lognormal()
#prob_missing_given_lognormal(p,median_sampling,stdev_mag,ncoll)
#lgnd <- 1-lgnd
#lgnd <- lgnd^ncoll
#return(mean(lgnd))
}

zipf_distribution <- function(ev, S)	{
ranks <- seq(1,S,by=1)
zpf <- ranks^-ev
prop <- zpf/sum(zpf)
return(prop)
}

scaled_zipf_distribution <- function(r,ev, S)	{
#	rescale zipf to base rate r
return(r*zipf_distribution(ev, S))
}

best_zipf_occupany_rates_set_S <- function(p0, hS, observed, ncoll)	{
# p0[1]: median rate p0[2]: log-log decay	p0[3]: richness
oS <- sum(observed);
rad <- scaled_zipf_distribution(p0[1],p0[2],hS)		# basic lognormal distribution
if (rad[1]<=1)	{
	# expected number of taxa with 1, 2, 3, etc. finds
	#	expected_prop <- expected_prop(rad,ncoll,p0[3])
	expected <- expected_occurrences(rad,ncoll,hS)
	# log likelihood
	lnl <- distribution_loglikelihood_mul(observed,expected,oS,hS)
	}	else {
	lnl <- 1000*-oS
	#lnl <- -1e+308
	}
return(lnl)
}

get_min_max_occupancy_zipf_params <- function(finds,observed,ncoll)	{
oS <- sum(finds>0);
rnk <- length(subset(observed,observed>0))
categoricals <- subset(observed,observed>0)[rnk:1]
categ_ranks2 <- categ_ranks <- vector(length=length(categoricals))
for (g in 1:categoricals[1])	categ_ranks2[1] <- categ_ranks2[1]+g
for (f in 2:length(categoricals))	{
	categ_ranks[f] <- categ_ranks[f-1]+categoricals[f-1]
	for (g in (categ_ranks[f]+1):(categ_ranks[f]+categoricals[f]))	categ_ranks2[f] <- categ_ranks2[f]+g
	}
for (f in 1:length(categoricals))	categ_ranks2[f] <- categ_ranks2[f]/categoricals[f]
slopes <- vector(length=(length(categoricals)-1))
categ_ranks <- categ_ranks+1
for (f in 1:(length(slopes)))	{
	slopes[f] <- log(finds[categ_ranks[f]]/finds[categ_ranks[f+1]])/log(categ_ranks2[f+1]/categ_ranks2[f])
	}

min_ev <- min(slopes)
max_ev <- max(slopes)
d1 <- zipf_distribution(min_ev,max(Chao2_Fisher(observed),jack2_Fisher(observed)))
d2 <- zipf_distribution(max_ev,oS)
m1 <- m2 <- 1
for (s in 1:ceiling(oS/20))	{
	m1 <- m1*(finds[s]/ncoll)/d1[s]
	m2 <- m2*(finds[s]/ncoll)/d2[s]
	}
minmidfnds <- min(m1^(1/ceiling(oS/20)),m2^(1/ceiling(oS/20)))
maxmidfnds <- max(m1^(1/ceiling(oS/20)),m2^(1/ceiling(oS/20)))

return(list(c(minmidfnds,min_ev),c(maxmidfnds,max_ev)));
}

get_initial_occupancy_zipf <- function(finds,observed,ncoll)	{
oS <- sum(finds>0);
p0 <- vector(length=3)
p0[3] <- smin <- min(Chao2_Fisher(observed),jack2_Fisher(observed))	# lowest richness estimator we will consider
rnk <- length(subset(observed,observed>0))
categoricals <- subset(observed,observed>0)[rnk:1]
categ_ranks2 <- categ_ranks <- vector(length=length(categoricals))
for (g in 1:categoricals[1])	categ_ranks2[1] <- categ_ranks2[1]+g
for (f in 2:length(categoricals))	{
	categ_ranks[f] <- categ_ranks[f-1]+categoricals[f-1]
	for (g in (categ_ranks[f]+1):(categ_ranks[f]+categoricals[f]))	categ_ranks2[f] <- categ_ranks2[f]+g
	}
for (f in 1:length(categoricals))	categ_ranks2[f] <- categ_ranks2[f]/categoricals[f]
slopes <- vector(length=(length(categoricals)-1))
categ_ranks <- categ_ranks+1
for (f in 1:(length(slopes)))	{
	slopes[f] <- log(finds[categ_ranks[f]]/finds[categ_ranks[f+1]])/log(categ_ranks2[f+1]/categ_ranks2[f])
	}

p0[2] <- prod(slopes)^(1/length(slopes))
d1 <- zipf_distribution(prod(slopes)^(1/length(slopes)),p0[3])
m1 <- 1
for (s in 1:ceiling(oS/20))	{
	m1 <- m1*(finds[s]/ncoll)/d1[s]
	}
p0[1] <- m1^(1/ceiling(oS/20))
names(p0) <- c("scale", "magnitude", "richness")
return(p0)
}

optimize_zipf_occupancy <- function(finds,ncolls)	{
# pO[1]=median occupancy rate; p0[2]= log-log decay; p0[3]=richness
observed <- fisher_plot(finds)
sO <- length(finds)	# observed taxa
oS <- sum(finds>0);
p0 <- get_initial_occupancy_zipf(finds,observed,ncoll)
# step 1: get rough estimate of lognormal based on variance above the median
smin <- hS <- p0[3]
dd <- get_min_max_occupancy_zipf_params(finds,observed,ncoll)
zipf_lo <- dd[[1]]
zipf_hi <- dd[[2]]
if (p0[1]<zipf_lo[1])	{
	zipf_lo[1] <- p0[1]
	} else if (p0[1]>zipf_hi[1]) {
	zipf_hi[1] <- p0[1]
	}
if (p0[2]<zipf_lo[2])	{
	zipf_lo[2] <- p0[2]
	} else if (p0[2]>zipf_hi[2]) {
	zipf_hi[2] <- p0[2]
	}

cl <- list(fnscale=-1)
w <- optim(p0[1:2],fn=best_zipf_occupany_rates_set_S,method="L-BFGS-B",hS=hS,observed=observed,ncoll=ncoll,lower=zipf_lo,upper=zipf_hi,control=cl)
lblnl <- blnl <- w$value
bH <- c(w$par,hS,w$value)
lgnHs <- c(w$par,w$value)
inc <- 1
while (lblnl==blnl && hS>oS)	{
	hS <- hS+inc
	w <- optim(p0[1:2],fn=best_zipf_occupany_rates_set_S,method="L-BFGS-B",hS=hS,observed=observed,ncoll=ncoll,lower=zipf_lo,control=cl)
	#lgnHs <- rbind(lgnHs,c(w$par,hS,w$value))
	if (w$value>blnl)	{
		bH <- c(w$par,hS,w$value)
		lblnl <- blnl <- w$value
		}	else if (hS==(smin+1))	{
		hS <- smin
		inc <- inc/-1
		lblnl=blnl
		}	else {
		lblnl <- w$value
		}
	}
bH <- c(bH,modified_AIC(blnl,3,oS))
names(bH) <- c("scale", "magnitude", "richness","loglikelihood","AICc")
return(bH)
}
