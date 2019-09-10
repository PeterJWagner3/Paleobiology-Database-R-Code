#install.packages("combinat", dependencies=TRUE)
#install.packages("gtools", dependencies=TRUE)
#install.packages("sads", dependencies=TRUE)
#install.packages("untb", dependencies=TRUE)
library(combinat)
library(gtools)
library(sads)
library(untb)
ZERO <- 1e-323
MINEXPN <- 10^-10
MINNO <- 5e-324
MAXNO <- 1.797693e+308
MAX_J <- 10^7.25

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
	aic_c <- AIC(lnL,k)
	}	else	aic_c <- (-2*lnL) + (2*k)*(n / (n - k - 1))
return(aic_c)
}

# convert abundances to "Fisher plot" giving # taxa with 1…N Finds
fisher_plot <- function(finds)	{
return(hist(finds,breaks=(0:max(finds)),plot=FALSE)$counts)
}
#unique_finds <- sort(unique(finds),decreasing=FALSE)
#observed <- vector(length=max(unique_finds))
#for (s in 1:length(unique_finds))	observed[unique_finds[s]] <- length(finds[finds==unique_finds[s]])
#return(observed)
#}

# Get Chao2 Richness estimate from vector giving # taxa with 1…N finds
Chao2 <- function(finds)	{
o1 <- length(subset(finds,finds==1))
o2 <- length(subset(finds,finds==2))
S <- round(length(finds) + ((o1*o1)/(2*(o2+1))) - ((o1*o2)/((2*(o2+1)*(o2+1)))))
return(S)
}

Chao2_Fisher <- function(observed)	{
ntaxa <- sum(observed)
S <- round(ntaxa + ((observed[1]*observed[1])/(2*(observed[2]+1))) - ((observed[1]*observed[2])/((2*(observed[2]+1)*(observed[2]+1)))))
return(S)
}

jackknife2 <- function(finds)	{
ntaxa <- length(finds)
#o1 <- length(subset(finds,finds==1))
#o2 <- length(subset(finds,finds==2))
ss <- 0
observed <- hist(finds,breaks=0:max(finds),plot=FALSE)$counts
for (i in 1:length(observed))	ss <- ss+(i*observed[i])
S <- round(ntaxa + (observed[1]*(((2*ss)-3)/ss)) - observed[2]*(((ss-2)*(ss-2))/(ss*(ss-1))))
return(S)
}

jackknife5 <- function(finds)	{
ntaxa <- length(finds)
#o1 <- length(subset(finds,finds==1))
#o2 <- length(subset(finds,finds==2))
ss <- 0
observed <- hist(finds,breaks=0:max(finds),plot=FALSE)$counts
for (i in 1:length(observed))	ss <- ss+(i*observed[i])
S <- round(ntaxa + (observed[1]*(((5*ss)-15)/ss)) - observed[2]*((10*(ss*ss)-(70*ss)+125)/(ss*(ss-1))) + observed[3]*(((10*(ss^3))-(120*(ss^2))+(485*ss)-660))/(ss*(ss-1)*(ss-2)) - observed[4]*((ss-4)^4)/(ss*(ss-1)*(ss-2)*(ss-3)) + observed[5]*((ss-5)^5)/(ss*(ss-1)*(ss-2)*(ss-3)*(ss-4)))
return(S)
}

jack2_Fisher <- function(observed)	{
ntaxa <- sum(observed)
ss <- 0
for (i in 1:length(observed))	ss <- ss+(i*observed[i])
S <- round(ntaxa + (observed[1]*(((2*ss)-3)/ss)) - observed[2]*(((ss-2)*(ss-2))/(ss*(ss-1))))
return(S)
}

jack5_Fisher <- function(observed)	{
ntaxa <- sum(observed)
ss <- 0
for (i in 1:length(observed))	ss <- ss+(i*observed[i])
S <- round(ntaxa + (observed[1]*(((5*ss)-15)/ss)) - observed[2]*((10*(ss*ss)-(70*ss)+125)/(ss*(ss-1))) + observed[3]*(((10*(ss^3))-(120*(ss^2))+(485*ss)-660))/(ss*(ss-1)*(ss-2)) - observed[4]*((ss-4)^4)/(ss*(ss-1)*(ss-2)*(ss-3)) + observed[5]*((ss-5)^5)/(ss*(ss-1)*(ss-2)*(ss-3)*(ss-4)))
return(S)
}

count_finds_in_assemblage <- function(assemblage,abundance_data)	{
return(sum(abundance_data$Counts[abundance_data$Assemblage==assemblage]))
}

finds_per_assemblage <- function(abundance_data)	{
assemblage <- sort(unique(abundance_data$Assemblage))
return(sapply(assemblage,count_finds_in_assemblage,abundance_data))
}

count_taxa_in_assemblage <- function(assemblage,abundance_data)	{
return(sum(abundance_data$Assemblage==assemblage))
}

taxa_per_assemblage <- function(abundance_data)	{
assemblage <- sort(unique(abundance_data$Assemblage))
return(sapply(assemblage,count_taxa_in_assemblage,abundance_data))
}

# find expectations of this gamma at this sample size
# exp will give the expected number of taxa found 1…ncoll times
expected_abundances <- function(rel_ab_dist, nspec, S)	{
expected <- vector(length=nspec)
for (t in 1:S)	{
	for (i in 1:nspec)	{
		expected[i] <- expected[i]+dbinom(i,nspec,rel_ab_dist[t])
		}
	}
return(expected)
}

# get the likelihood of observed numbers of taxa with 1…N finds given expected numbers of taxa with 1…N finds
distribution_loglikelihood_mul <- function(observed,expected,oS,hS)	{
#print(c(oS,hS))		# for debugging
mxfind <- length(observed)							# maximum finds observed
prop_expected <- expected[1:mxfind]/sum(expected)		# convert expected species to proportions
prop_expected[prop_expected==0] <- MINNO
lnlo <- observed*log(prop_expected)	# exact probability of observing # taxa with 1, 2, 3, etc. finds
lnlo[is.na(lnlo)] <- 0
#for (i in 1:mxfind)	if (is.na(lnlo[i]))	lnlo[i] <- 0
sobs <- sum(lnlo)
# log probability of observing X taxa from hypothesized hS taxa given expected # taxa with 1, 2, 3, etc. finds
eS <- sum(expected)									# get expected sampled species
while (eS==hS)	eS <- 0.99999*eS
	# this is a kluge to get around rounding error of eS->hS

if (hS>=oS && round(eS,5)<round(hS,5))	{
	# get the probability of observing oS of hS species given that we expected to observe eS of hS species
	lnls <- lfactorial(hS)-(lfactorial(hS-oS)+lfactorial(oS))+(oS*log(eS/hS))+((hS-oS)*log((hS-eS)/hS))
	# this should be the norm: more species than observed hypothesized
	}	else	{
	lnls <- oS*log(MINNO)
	# if hypothesized true richness is less than observed, then this is impossible
	}
return(sobs+lnls)
}

# given a "Preston Plot" with X taxa with abundance 1…N, return a rank abundance distribution
#	Good for converting log-series & zero sum to RADs
convert_expected_to_relative_abundance <- function(taxa_w_n_finds)	{
#cum_expect <- floor(cumsum(taxa_w_n_finds))	# get the cumulative # species with n finds
cum_expect <- round(cumsum(taxa_w_n_finds),0)	# get the cumulative # species with n finds
uniq_expect <- c(unique(cum_expect))			# get the unique values of n finds
crit_counts <- match(uniq_expect,cum_expect)	# actual abudances
raw_abund <- vector(length=max(uniq_expect))
raw_abund[1:uniq_expect[1]] <- crit_counts[1]
if (length(crit_counts)>1)
	for (i in 2:length(crit_counts))
		raw_abund[(uniq_expect[i-1]+1):uniq_expect[i]] <- crit_counts[i]
return(sort(raw_abund/sum(raw_abund),decreasing=TRUE))
}

# when given a vector 2 4 5 6 6 7, return that there are 2 2 1 1 0 1 entities in each cell
entities_in_rank  <- function(entities)	{
dummy <- c(0,entities)
members <- entities
for (i in 1:length(entities))	{
	members[i] <- entities[i]-dummy[i]
	}
return(members)
}

# estimate how many specimens an sad should have given P[finds] and S taxa
N_given_S_and_Prob_N <- function(prob_n_finds,S)	{
dope <- floor(0.5+cumsum(S*prob_n_finds))
uniq_dope <- unique(dope)[unique(dope)>0]
dope_finds <- match(uniq_dope,dope)
dope_mems <- entities_in_rank(uniq_dope)
return(sum(dope_mems*dope_finds))
}

optimize_uniform_abundance <- function(counts)	{
# written 2017-01-28
minS <- length(counts)
nspec <- sum(counts)
oS <- hS <- minS
observed <- fisher_plot(counts)
mxlnl <- lnl <- -1*MAXNO
while (lnl == mxlnl)	{
	rel_ab_dist <- rep(1/hS,hS)
	raw_expected <- expected_abundances(rel_ab_dist,nspec,hS)
	lnl <- distribution_loglikelihood_mul(observed,expected=raw_expected,oS,hS)
	if (mxlnl < lnl)	mxlnl <- lnl
	hS <- hS+1
	}
bH <- c(hS-1,round(mxlnl,3),round(modified_AIC(mxlnl,1,nspec),3))
names(bH) <- c("Uniform_S","Uniform_log-likelihood","Uniform_AICc")
return(bH)
}

# generate geometric (= exponential) distribution
geometric_distribution <- function(decay)	{
if (decay>1)	decay <- (1/decay)
S <- round(1+((log(10^-9)-log(1-decay))/log(decay)))
ranks <- (1:S)
prop <- (1-decay)*(decay^(ranks-1))
return(prop)
}

# get log-likelihood of geometric given decay
loglikelihood_geometric_rad <- function(decay,nspec,oS,observed)	{
# p0[1]: r	# p0[2]: decay	# p0[3]: S
rel_ab_dist <- geometric_distribution(decay)		# basic exponential distribution
hS <- length(rel_ab_dist)
#print(c(round(decay,5),hS))		# for debugging
if (hS>=oS)	{
	raw_expected <- expected_abundances(rel_ab_dist,nspec,hS)
		# log likelihood
	lnl <- distribution_loglikelihood_mul(observed,expected=raw_expected,oS,hS)
	}	else	{
	lnl <- oS*log(MINNO)
	}
#print(c(round(decay,5),round(lnl,10)))		# for debugging
return(lnl)
}

# get best-fit geometric
optimize_geometric_abundance <- function(counts)	{
oS <- length(counts)				# observed taxa
observed <- fisher_plot(counts)
decay <- (min(counts)/max(counts))^(1/(length(counts)-1))
max_decay <- exp(log(decay)/2)
min_decay <- exp(log(decay)*2)

cl <- list(fnscale=-1)
nspec <- sum(counts)
w <- optim(decay,fn=loglikelihood_geometric_rad,method="L-BFGS-B",nspec=nspec,oS=oS,observed=observed,lower=min_decay,upper=max_decay,control=cl)
bH <- c(round(w$par,6),round(w$value,3),round(modified_AIC(w$value,1,nspec),3))
names(bH) <- c("Geometric_decay","Geometric_log-likelihood","Geometric_AICc")
return(bH)
}

# estimate log-series distribution given alpha and N
logseries_distribution <- function(alpha,J)	{
# J is the true (incompletely sampled) population size
# alpha = Fisher's Biodiversity Index
prob_n_finds <- dls((1:J),J,alpha)	# dls from sads
N <- 0
S <- 2
while (N<J)	{
	S <- S+1
	N <- N_given_S_and_Prob_N(prob_n_finds,S)
	}
taxa_w_n_finds <- S*prob_n_finds
sad <- convert_expected_to_relative_abundance(taxa_w_n_finds)
return (sad)
}

# get log-likelihood of log-series given alpha & J
loglikelihood_logseries_rad <- function(oS,nspec,alpha,J,observed)	{
#print(round(alpha,10))		# for debugging
rel_ab_dist <- logseries_distribution(alpha=alpha,J=J)		# 
hS <- length(rel_ab_dist)
if (hS>=oS)	{
	raw_expected <- expected_abundances(rel_ab_dist,nspec,hS)
#	lnl <- abundance_distribution_loglikelihood_suf(observed,expected=raw_expected,oS=oS,hS=hS)
	lnl <- distribution_loglikelihood_mul(observed,expected=raw_expected,oS,hS)
	}	else {
	lnl <- oS*log(MINNO)
	}
return(lnl)
}

# get minimum alpha for log-series given observed taxa oS & population size J
accio_minimum_logseries_alpha_given_S_and_J <- function(oS,J)	{
alpha <- 1
prob_n_finds <- dls((1:J),J,alpha)	# dls from sads
prob_n_finds <- prob_n_finds/sum(prob_n_finds)
N <- N_given_S_and_Prob_N(prob_n_finds,oS)
if (N>J)	{
	while (N>J)	{
		alpha <- alpha*1.1
		prob_n_finds <- dls((1:J),J,alpha)	# dls from sads
		prob_n_finds <- prob_n_finds/sum(prob_n_finds)
		N <- N_given_S_and_Prob_N(prob_n_finds,oS)
		}
	}	else	{
	last_N <- N
	run_N <- 0
	while (N<J && run_N<10)	{
		alpha <- alpha/1.1
		prob_n_finds <- dls((1:J),J,alpha)	# dls from sads
		prob_n_finds <- prob_n_finds/sum(prob_n_finds)
		N <- N_given_S_and_Prob_N(prob_n_finds,oS)
		if (N==last_N)	{
			run_N <- run_N+1
			}	else	{
			last_N <- N
			run_N <- 0
			}
#		print(c(alpha,N))		# for debugging
		}
	}
return(alpha)
}

# get the best alpha given a total population size
optimize_logseries_abundance_given_J <- function(J,observed,nspec,oS,init_alpha,apprise=FALSE)	{
# add something to get minimum alpha
if (apprise)	print(paste("Log-Series J = ",J," ",date(),sep=""))		# for updating, as logseries is slow!
#if (apprise)	print(paste("Log-Series J = ",J,sep=""))			# for updating, as logseries is slow!
#min_alpha <- accio_minimum_logseries_alpha_given_S_and_J(oS,J)
cl <- list(fnscale=-1)
alpha <- init_alpha
w <- optim(alpha,fn=loglikelihood_logseries_rad,method="L-BFGS-B",oS=oS,nspec=nspec,observed=observed,J=J,lower=init_alpha/10,upper=10*init_alpha,control=cl)
lsj_output <- c(round(w$par,5),J,round(w$value,3))
names(lsj_output) <- c("alpha","J","lnL")
return(lsj_output)
}

# get best-fit log-series
#powers <- c(4.531478917,5.031477138,5.531478917,6.031478754)
optimize_logseries_abundance <- function(counts,byte_trail=FALSE,span=4, max_J=MAX_J,apprise=FALSE)	{
oS <- length(counts)
observed <- fisher_plot(counts)
base <- 10
peak <- 0
nspec <- sum(counts)
exp_st <- log10(nspec)
incr <- 0.5
powers <- exp_st+incr*(0:span)
exp_en <- max(powers)
pa <- 1
pz <- length(powers)
init_alpha <- fishers.alpha(sum(counts),oS)
if (byte_trail)	byte_trail_file <- paste("Log_Series_Trail_oS=",oS,"_N=",nspec,".txt",sep="")
while (incr>0.001) {
	J <- round(base^powers)
	results <- sapply(J[pa:pz],optimize_logseries_abundance_given_J,oS=oS,nspec=nspec,observed=observed,init_alpha=init_alpha,apprise=apprise)
	if (byte_trail)	{
		if (pa==2)	{
			byte_trail_results <- rbind(byte_trail_results,t(results))
			}	else {
			byte_trail_results <- t(results)
			}
		write.table(byte_trail_results,byte_trail_file,sep="\t",col.names=TRUE,row.names=FALSE)
		}
	if (pa==2)	
		results <-cbind(old_result_1,results)
	if (pz<length(powers))
		results <-cbind(results,old_result_2)
	mlln <- max(results[3,])
	mlJc <- match(mlln,results[3,])
	mlJ <- results[2,mlJc]
	bH <- results[,mlJc]
	if (mlJc==length(results[3,]))	{
		if (peak==0)	{
			lower <- powers[mlJc]
			upper <- min(log10(max_J),(lower + (incr*span)))
			incr <- (upper-lower)/span
			powers <- lower+incr*(0:span)
			old_result_1 <- results[,mlJc]
			pa <- 2
			# we are still improving by increasing J
			}	else	{
			# already found a peak J, but the higest in this search is the best
			lower <- powers[mlJc-1]
			upper <- powers[mlJc]
			incr <- (upper-lower)/span
			powers <- lower+incr*(0:span)
			old_result_1 <- results[,mlJc-1]
			old_result_2 <- results[,mlJc]
			pa <- 2
			pz <- length(powers) - 1
			}
		# end case where highest J is the best
		}	else if (mlJc==1)	{
			# the lowest examined J is best
		if (peak==0)	{
			upper <- powers[mlJc+1]
			lower <- powers[mlJc]
			old_result_1 <- results[,mlJc]
			old_result_2 <- results[,mlJc+1]
			pa <- 2
			pz <- length(powers) - 1
			incr <- (upper-lower)/span
			powers <- lower+incr*(0:span)
			peak <- 1
			}	else {
			peak <- 1
			lower <- powers[mlJc]
			upper <- powers[mlJc+1]
			incr <- (upper-lower)/span
			powers <- lower+incr*(0:span)
			old_result_1 <- results[,mlJc]
			old_result_2 <- results[,mlJc+1]
			pa <- 2
			pz <- length(powers) - 1
			}
		# end case where lowest J is the best
		}	else	{
		# some J in the middle is the best
		peak <- 1
		J2 <- round(c(base^(powers[mlJc]-(incr/10)),base^(powers[mlJc]+(incr/10))),0)
		results2 <- sapply(J2,optimize_logseries_abundance_given_J,oS=oS,nspec=nspec,observed=observed,init_alpha=init_alpha,apprise=apprise)
		if (byte_trail)	{
			if (pa==2)	{
				byte_trail_results <- rbind(byte_trail_results,t(results))
				}	else {
				byte_trail_results <- t(results)
				}
			write.table(byte_trail_results,byte_trail_file,sep="\t",col.names=TRUE,row.names=FALSE)
			}
		if (results2[3,1]>results2[3,2])	{
			# lower J improves likelihood
			lower <- powers[mlJc-1]
			old_result_1 <- results[,mlJc-1]
			if (results2[3,1]>results[3,mlJc])	{
				# if lower J is better than bH, then reset
				upper <- powers[mlJc]-(incr/10)
				bH <- old_result_2 <- results2[,1]
				}	else {
				upper <- powers[mlJc]
				old_result_2 <- results[,mlJc]
				}
			incr <- (upper-lower)/span
			powers <- lower+incr*(0:span)
			} else	{
			# higher J improves likelihood
			upper <- powers[mlJc+1]
			old_result_2 <- results[,mlJc+1]
			if (results2[3,2]>results[3,mlJc])	{
				# if higher J is better than bH, then reset
				lower <- powers[mlJc]+(incr/10)
				bH <- old_result_1 <- results2[,2]
				}	else {
				lower <- powers[mlJc]
				old_result_1 <- results[,mlJc]
				}
			incr <- (upper-lower)/span
			powers <- lower+incr*(0:span)
			}
		pa <- 2
		pz <- length(powers)-1
		}
	}
bH <- c(bH,modified_AIC(bH[3],2,nspec))
names(bH) <- c("Log-Series_Alpha","Log-Series_J", "Log-Series_loglikelhood", "Log-Series_AICc")
return(bH)
}

# get best zero sum multinomial for all three variables
optimo_logseries_abundance <- function(hubbells,observed,nspec,oS,min_alpha,max_alpha,min_J,max_J)	{
test_alpha <- min_alpha+(hubbells[1]*(max_alpha-min_alpha))
test_J <- round((min_J+(hubbells[2]*(max_J-min_J))),0)
rel_ab_dist <- logseries_distribution(test_J,test_alpha)
hS <- length(rel_ab_dist)
if (hS>=oS)	{
	raw_expected <- expected_abundances(rel_ab_dist,nspec,S=hS)
	lnl <- distribution_loglikelihood_mul(observed,expected=raw_expected,oS=oS,hS=hS)
	}	else {
	lnl <- oS*log(MINNO)
	}
return(lnl)
}

# calcualate lognormal for S entities that increase in magnitude by mag every standard deviation
lognormal_distribution <- function(mag, S)	{
fi <- (S:1)/(S+1)
prop <- mag^(qnorm(fi,0,1))/sum(mag^(qnorm(fi,0,1)))
return(prop)
}

# get log-likelihood for lognormal
loglikelihood_lognormal_rad <- function(oS,nspec,mag,hS,observed)	{
#print(c(round(mag,9)))
rel_ab_dist <- lognormal_distribution(mag=mag,S=hS)		# basic lognormal distribution
raw_expected <- expected_abundances(rel_ab_dist,nspec,S=hS)
#lnl <- abundance_distribution_loglikelihood_suf(observed,expected=raw_expected,oS=oS,hS=hS)
lnl <- distribution_loglikelihood_mul(observed,expected=raw_expected,oS,hS)
return(round(lnl,2))
}

# get log-likelihood for lognormal
loglikelihood_lognormal_rad_for_optim <- function(oS,nspec,rand_no,hS,observed,min_mag,max_mag)	{
mag <- min_mag + rand_no*(max_mag - min_mag)
#print(c(round(mag,9)))
rel_ab_dist <- lognormal_distribution(mag=mag,S=hS)		# basic lognormal distribution
raw_expected <- expected_abundances(rel_ab_dist,nspec,S=hS)
#lnl <- abundance_distribution_loglikelihood_suf(observed,expected=raw_expected,oS=oS,hS=hS)
lnl <- distribution_loglikelihood_mul(observed,expected=raw_expected,oS,hS)
return(round(lnl,2))
}

# get minimum and maximum magnitude of increase for a lognormal given hS
accio_min_and_max_lognormal_mag_given_hS <- function(observed,counts,hS)	{
oS <- sum(observed)				# observed taxa
mxfinds <- length(observed)
mnfinds <- min((1:mxfinds)[!observed %in% 0])
rank_shifts <- sort(qnorm((1:hS)/(hS+1)),decreasing=TRUE)[1:oS]
local_m <- mag_shifts <- vector(length=(oS-1))
for (i in 1:(oS-1))		{
	mag_shifts[i] <- counts[i]/counts[i+1]
	local_m[i] <- exp(log(mag_shifts[i])/(rank_shifts[i]-rank_shifts[i+1]))
	}
if (max(observed)==1)	{
	min_mag <- min(local_m)
	}	else {
	min_mag <- 1.01
	}
max_mag <- prod(local_m)^(1/(length(local_m)))
if (max_mag < min_mag)	{
	m <- min_mag
	min_mag <- max_mag
	max_mag <- m
	}
return(c(min_mag,max_mag))
}

# get most likely magnitude of increase for a lognormal distribution of hS entities
optimize_lognormal_abundance_given_hS <- function(observed,counts,hS)	{
#print(hS)		# for debugging
oS <- sum(observed)				# observed taxa
nspec <- sum((1:length(observed))*observed)
mxfinds <- length(observed)
mnfinds <- min((1:mxfinds)[!observed %in% 0])
mm <- accio_min_and_max_lognormal_mag_given_hS(observed,counts,hS)
mag <- prod(mm)^0.5
#inev <- mag <- exp(log(mxfinds/mnfinds)/(qnorm((hS-1)/(hS+1))-qnorm((hS-oS)/(hS+1))))
cl <- list(fnscale=-1)
rand_no <- (mag - mm[1])/(mm[2] - mm[1])
w <- optim(rand_no,loglikelihood_lognormal_rad_for_optim,,method="L-BFGS-B",oS=oS,nspec=nspec,hS=hS,observed=observed,min_mag=mm[1],max_mag=mm[2],lower=0,upper=1,control=cl)
w$par <- mm[1] + (w$par*(mm[2] - mm[1]))
#w <- optim(mag,fn=loglikelihood_lognormal_rad,method="L-BFGS-B",oS=oS,nspec=nspec,hS=hS,observed=observed,lower=mm[1],upper=mm[2],control=cl)
bH <- c(w$par,hS,w$value)
names(bH) <- c("Lognormal_magnitude", "Lognormal_S","Lognormal_log-likelihood")
return(bH)
}

# get most likely lognormal distribution for a population of some sort
optimize_lognormal_abundance <- function(counts)	{
stS <- oS <- length(counts)				# observed taxa
observed <- fisher_plot(counts)
nspec <- sum(counts)
span <- 5
enS <- stS+((span-1)*oS)
incr <- floor(enS/span)
cl <- list(fnscale=-1)
peak <- 0
while (incr>0)	{
	hS <- seq(stS,enS,by=incr)
	results <- sapply(hS,optimize_lognormal_abundance_given_hS,observed=observed,counts=counts)
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
				stS <- hS[spn]
				enS <- stS+((spn-1)*oS)
				}	else	{
				enS <- hS[spn]
				if (span<incr)	{
					incr <- floor((enS-hS[spn-1])/span)
					}	else {
					span <- enS-hS[spn-1]
					incr <- 1
					}
				stS <- hS[spn]-(span*incr)
				}	# end case where last number is best, but this is after finding a peak.
			} else if (mlSc==1)	{
			# if first richness is the best
			peak <- 1
			stS <- hS[mlSc]
			enS <- hS[mlSc+1]-1
			if (incr==1)	incr <- 0
			if (span<incr)	{
				incr <- floor((enS-stS)/span)
				}	else {
				span <- enS-stS
				incr <- 1
				}
			} else {
			# if a richness in the middle is best & we still can find a better one
			peak <- 1
			hS2 <- c(mlS-1,mlS+1)
			results2 <- sapply(hS2,optimize_lognormal_abundance_given_hS,observed=observed,counts=counts)
			if (results2[3,1]>mlnl && results2[3,1]>results2[3,2])	{
			# start just above 2nd best richness and go up to best
				stS <- hS[mlSc-1]+1
				enS <- hS[mlSc]
				if (span<incr)	{
					incr <- floor((enS-stS)/span)
					}	else {
					span <- enS-stS
					incr <- 1
					}
				}	else if (results2[3,2]>mlnl && results2[3,2]>results2[3,1]) {
			# start at best richness and go just below 2nd best
				stS <- hS[mlSc]
				enS <- hS[mlSc+1]-1
				if (incr>span)	{
					incr <- floor((enS-stS)/span)
					}	else {
					span <- enS-stS
					incr <- 1
					}
				}	else {
				# we already had the best, so just end it
				incr <- 0
				}
			}
		# end case where we have a better richness in middle somewhere.  
		}	else	{
		incr <- 0
		}
	}
bH <- c(round(bH[1],6),bH[2],bH[3],modified_AIC(mlnl,2,nspec))
names(bH)[4] <- "Lognormal_AICc"
return(bH)
}

# generate Zipf distribution for S individuals with log-decay = zg
zipf_distribution <- function(zg, S)	{
return(((1:S)^-zg)/sum((1:S)^-zg))
#zpf <- ranks^-zg
#prop <- zpf/sum(zpf)
#return(prop)
}

# get the likelihood of a Zipf model with log-log decay = zg & hS entities given oS observed with nspec individuals
loglikelihood_zipf_rad <- function(oS,nspec,zg,hS,observed)	{
#print(zg)
rel_ab_dist <- zipf_distribution(zg=zg,S=hS)		# basic zipf distribution
raw_expected <- expected_abundances(rel_ab_dist,nspec,S=hS)
lnl <- distribution_loglikelihood_mul(observed,expected=raw_expected,oS=oS,hS=hS)
return(round(lnl,3))
}

# get the most likely Zipf model for a particular richness
optimize_zipf_abundance_given_hS <- function(counts,observed,max_zipf,hS)	{
#print(hS)		# for debugging
oS <- sum(observed)				# observed taxa
nspec <- sum((1:length(observed))*observed)
mxfinds <- length(observed)
mnfinds <- min((1:mxfinds)[!observed %in% 0])
inzg <- zg <- log(mxfinds/mnfinds)/log(oS-1)
cl <- list(fnscale=-1)
w <- optim(zg,fn=loglikelihood_zipf_rad,method="L-BFGS-B",oS=oS,nspec=nspec,hS=hS,observed=observed,lower=0,upper=max_zipf,control=cl)
bH <- c(w$par,hS,w$value)
names(bH) <- c("Zipf_log-log_decay","Zipf_hypothesized_richness","Zipf_log-likelihood")
return(bH)
}

# get the maximum zipf decay implied by the distribution of data points
accio_max_zipf_decay <- function(counts,oS,nspec)	{
xx <- log(1:oS)
yy <- log(counts/nspec)-log(counts[1]/nspec)
zz <- -1*(yy/xx)[2:length(xx)]
return(max(zz))
}

# get the most likely Zipf model
optimize_zipf_abundance <- function(counts,span=5)	{
stS <- oS <- length(counts)				# observed taxa
observed <- fisher_plot(counts)
enS <- stS+((span-1)*oS)
incr <- floor(enS/span)
cl <- list(fnscale=-1)
peak <- 0
nspec <- sum(counts)
max_zipf <- 2*accio_max_zipf_decay(counts,oS,nspec)
pa <- 1
pz <- span
while (incr>0)	{
	hS <- seq(stS,enS,by=incr)
	results <- sapply(hS[pa:pz],optimize_zipf_abundance_given_hS,counts=counts,observed=observed,max_zipf=max_zipf)
	if (pa==2)	{
		results <- cbind(bH,results)
		} else if (pz<span)	{
		results <- cbind(results,bH)
		}
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
				stS <- hS[spn]
				pa <- 2
				pz <- span
				enS <- stS+((spn-1)*oS)
				}	else	{
				enS <- hS[spn]
				if (span<incr)	{
					incr <- floor((enS-hS[spn-1])/span)
					}	else {
					span <- enS-hS[spn-1]	# check this! TDay
					incr <- 1
					}
				stS <- hS[spn]-(span*incr)
				pa <- 1
				pz <- span-1
				}	# end case where last number is best, but this is after finding a peak.
			} else if (mlSc==1)	{
			# if first richness is the best
			peak <- 1
			stS <- hS[mlSc]
			enS <- hS[mlSc+1]-1
			if (incr==1)	incr <- 0
			if (span<incr)	{
				incr <- floor((enS-stS)/span)
				}	else {
				span <- enS-stS
				incr <- 1
				}
			pa <- 2
			pz <- span
			} else {
			# if a richness in the middle is best & we still can find a better one
			peak <- 1
			hS2 <- c(mlS-1,mlS+1)
			results2 <- sapply(hS2,optimize_zipf_abundance_given_hS,counts=counts,observed=observed,max_zipf=max_zipf)
			if (results2[3,1]>mlnl && results2[3,1]>results2[3,2])	{
			# start just above 2nd best richness and go up to best
				bH <- results2[,1]
				stS <- hS[mlSc-1]+1
				enS <- hS2[1]
				if (span<incr)	{
					incr <- floor((enS-stS)/span)
					}	else {
					span <- 1+enS-stS
					incr <- 1
					}
				pa <- 2
				pz <- span
				}	else if (results2[3,2]>mlnl && results2[3,2]>results2[3,1]) {
			# start at best richness and go just below 2nd best
				stS <- hS2[2]
				enS <- hS[mlSc+1]-1
				bH <- results2[,2]
				if (incr>span)	{
					incr <- floor((enS-stS)/span)
					}	else {
					span <- 1+enS-stS
					incr <- 1
					}
				pa <- 1
				pz <- span-1
				}	else {
				# we already had the best, so just end it
				incr <- 0
				}
			}
		# end case where we have a better richness in middle somewhere.  
		}	else	{
		incr <- 0
		}
	}
bH <- c(round(bH[1],6),bH[2],bH[3],modified_AIC(mlnl,k=2,n=nspec))
names(bH)[4] <- "Zipf_AICc"
return(bH)
}

# generate Zipf distribution for S individuals with log-decay = zg
zipf_mandelbrot_distribution <- function(zg,zb,S)	{
return((1/((1:S)+zb)^zg)/sum(1/((1:S)+zb)^zg))
}

# get the likelihood of a Zipf model with log-log decay = zg & hS entities given oS observed with nspec individuals
loglikelihood_zipf_mandelbrot_rad <- function(oS,nspec,zg,zb,hS,observed)	{
if (zg<0 || zb<0)	{
	return(-1*MAXNO)
	}	else	{
	if (zb>0 && zg>0)	{
		rel_ab_dist <- zipf_mandelbrot_distribution(zg=zg,zb=zb,S=hS)	# basic zipf-mandelbrot distribution
		} else if (zb==0 && zg>0)	{
		rel_ab_dist <- zipf_distribution(zg=zg,S=hS)			# basic zipf distribution
		} else if (zg==0)	{
		rel_ab_dist <- rel_ab_dist <- rep(1/hS,hS)
		}
	raw_expected <- expected_abundances(rel_ab_dist,nspec,S=hS)
	lnl <- distribution_loglikelihood_mul(observed,expected=raw_expected,oS=oS,hS=hS)
#	print(c(zg,zb,hS,lnl))
	return(round(lnl,3))
	}
}

loglikelihood_zipf_mandelbrot_rad_both <- function(zm,oS,nspec,hS,observed)	{
#zm: vector giving gamma and beta parameters
zg <- zm[1]	# gamma parameter (log-log decay)
zb <- zm[2]	# beta parameter
#print(c(zg,zb,hS))
if (zg<0 || zb<0)	{
	return(-1*MAXNO)
	}	else	{
	if (zb>0 && zg>0)	{
		rel_ab_dist <- zipf_mandelbrot_distribution(zg=zg,zb=zb,S=hS)	# basic zipf-mandelbrot distribution
		} else if (zb==0 && zg>0)	{
		rel_ab_dist <- zipf_distribution(zg=zg,S=hS)			# basic zipf distribution
		} else if (zg==0)	{
		rel_ab_dist <- rel_ab_dist <- rep(1/hS,hS)
		}
	raw_expected <- expected_abundances(rel_ab_dist,nspec,S=hS)
	lnl <- distribution_loglikelihood_mul(observed,expected=raw_expected,oS=oS,hS=hS)
#	print(c(zg,zb,hS,lnl))
	return(round(lnl,3))
	}
}

loglikelihood_zipf_mandelbrot_rad_zb <- function(zb,oS,nspec,zg,hS,observed)	{
if (zg<0 || zb<0)	{
	return(-1*MAXNO)
	}	else	{
	if (zb>0 && zg>0)	{
		rel_ab_dist <- zipf_mandelbrot_distribution(zg=zg,zb=zb,S=hS)	# basic zipf-mandelbrot distribution
		} else if (zb==0 && zg>0)	{
		rel_ab_dist <- zipf_distribution(zg=zg,S=hS)			# basic zipf distribution
		} else if (zg==0)	{
		rel_ab_dist <- rel_ab_dist <- rep(1/hS,hS)
		}
	raw_expected <- expected_abundances(rel_ab_dist,nspec,S=hS)
	lnl <- distribution_loglikelihood_mul(observed,expected=raw_expected,oS=oS,hS=hS)
#	print(c(zg,zb,hS,lnl))
	return(round(lnl,3))
	}
}

optimize_zipf_mandelbrot_abundance_given_hS_and_zm <- function(zg,hS,observed,init_beta=-1,max_beta=100)	{
#print(zg)
oS <- sum(observed)				# observed taxa
nspec <- sum((1:length(observed))*observed)
incr <- 0.1
if (init_beta==-1)	{
	zb <- seq(0,2.0,by=incr)
	} else	{
	zb <- seq(init_beta-1,init_beta+1.0,by=incr)
	}
zb_lnl <- sapply(zb,loglikelihood_zipf_mandelbrot_rad_zb,oS,nspec,zg,hS,observed)
ml_zb_c <- match(max(zb_lnl),zb_lnl)
ml_zb <- max(zb_lnl)

improve <- 1
while (ml_zb_c==length(zb_lnl) && (improve>0.01 && min(zb) < max_beta))	{
	a <- zb[length(zb_lnl)-1]
	b <- max(zb)+(max(zb)-min(zb))
	zb <- seq(a,b,by=incr)
	zb_lnl <- sapply(zb,loglikelihood_zipf_mandelbrot_rad_zb,oS,nspec,zg,hS,observed)
	ml_zb_c <- match(max(zb_lnl),zb_lnl)
	if (max(zb_lnl) < ml_zb)	{
		if (ml_zb_c==1)	{
			ml_zb_c <- 2
			} else if (ml_zb_c==length(zb))	{
			ml_zb_c <- length(zb)-1
			}
		improve <- min(abs(xx[ml_zb_c,2]-xx[ml_zb_c-1,2]),abs(xx[ml_zb_c,2]-xx[ml_zb_c+1,2]))
#		improve <- max(zb_lnl) - ml_zb
		ml_zb <- max(zb_lnl)
		}
	}

#if (ml_zb_c==1)
#	ml_zb_c <- 2

bzb <- zb[ml_zb_c]
while (improve>0.01 && min(zb) < max_beta)	{
	lml_zb_l <- lml_zb
	lml_zb <- ml_zb
	incr <- incr/10
	zb <- seq(zb[ml_zb_c]-(10*incr),zb[ml_zb_c]+(10*incr),by=incr)
	zb_lnl <- sapply(zb,loglikelihood_zipf_mandelbrot_rad_zb,oS,nspec,zg,hS,observed)
	if (ml_zb < max(zb_lnl))	{
		ml_zb_c <- match(max(zb_lnl),zb_lnl)
		bzb <- zb[ml_zb_c]
		ml_zb <- max(zb_lnl)
		if (ml_zb_c==1)	{
			ml_zb_c <- 2
			} else if (ml_zb_c==length(zb))	{
			ml_zb_c <- length(zb)-1
			}
		improve <- min(abs(zb_lnl[ml_zb_c]-zb_lnl[ml_zb_c-1]),abs(zb_lnl[ml_zb_c]-zb_lnl[ml_zb_c+1]))
		}	else	{
		improve <- 0
		}
#	if (ml_zb_c==1)	ml_zb_c <- 2
	}
#print(c(zg,bzb,ml_zb))
#return(ml_zb)
names(bzb) <- "Z-M_beta"
names(ml_zb) <- "lnl_Z-M"
return(c(bzb,ml_zb))
}

# get the most likely Zipf-Mandelbrot model for a particular richness
optimize_zipf_mandelbrot_abundance_given_hS <- function(hS,counts,observed,init_zg=-1)	{
# hS: hypothesized true richness
# counts: # finds per observed species
# observed: vector of length N giving # taxa with 1…N finds
# init_zg: log-log decay from best Zipf hypothesis
oS <- sum(observed)			# observed taxa
mxfinds <- length(observed)
mnfinds <- min((1:mxfinds)[!observed %in% 0])
if (init_zg==-1)	{
	inzg <- log(mxfinds/mnfinds)/log(oS-1)
	} else	{
	inzg <- init_zg
	}
incrm <- 1
zg <- seq(inzg-floor(inzg),inzg+floor(inzg),by=incrm)
while (length(zg)<3)	{
	incrm <- incrm/2
	zg <- seq(inzg-floor(inzg),inzg+floor(inzg),by=incrm)
	}
xx <- base::t(sapply(zg,optimize_zipf_mandelbrot_abundance_given_hS_and_zm,hS=hS,observed=observed))
lnl_ml_zm <- max(xx[,2])
improve <- 1
while (match(max(xx[,2]),xx[,2])==nrow(xx) && improve>=0.01)	{
	init_beta <-xx[nrow(xx),1]
	o_zm <- zg[(length(zg)-1):length(zg)]
	o_xx <- xx[((length(zg)-1):length(zg)),]
	a <- zg[length(zg)]
	b <- a+(max(zg)-min(zg))
	zg <- seq((a+incrm),b,by=incrm)
	xx <- base::t(sapply(zg,optimize_zipf_mandelbrot_abundance_given_hS_and_zm,hS=hS,observed=observed,init_beta=init_beta))
	xx <- rbind(o_xx,xx)
	zg <- c(o_zm,zg)
	mlsc <- match(max(xx[,2]),xx[,2])
	if (mlsc==1)	{
		mlsc <- 2
		} else if (mlsc==length(zg))	{
		mlsc <- length(zg)-1
		}
	if ((xx[mlsc,2]-xx[mlsc-1,2]) < improve)
		improve <- min(abs(xx[mlsc,2]-xx[mlsc-1,2]),abs(xx[mlsc,2]-xx[mlsc+1,2]))
	if (lnl_ml_zm < max(xx[,2]))	{
		lnl_ml_zm <- max(xx[,2])			# log-likelihood of ML (so far) model
		ml_zm_c <- match(lnl_ml_zm,xx[,2])	# cell with ML (so far) values
		ml_zm_g <- zg[ml_zm_c]				# ml log-log decay (gamma)
		}
	}

while (improve>=0.01)	{
	if (ml_zm_c > 1)	{
		incrm <- (zg[ml_zm_c]-zg[ml_zm_c-1])/2
		} else	{
		incrm <- (zg[ml_zm_c+1]-zg[ml_zm_c])/2
		}
	if (ml_zm_c<length(zg))	{
		zg <- seq(zg[ml_zm_c]-(2*incrm),zg[ml_zm_c]+(2*incrm),by=incrm)
		} else	{
		zg <- seq(zg[ml_zm_c]-(2*incrm),zg[ml_zm_c]+(2*incrm),by=incrm)
		#zg <- seq(zg[ml_zm_c-1],zg[ml_zm_c]+(2*incrm),by=incrm)
		}
	xx <- base::t(sapply(zg,optimize_zipf_mandelbrot_abundance_given_hS_and_zm,hS=hS,observed=observed))
	if (lnl_ml_zm <= max(xx[,2]))	{
#		improve <- max(xx[,2])-lnl_ml_zm
		mlsc <- match(max(xx[,2]),xx[,2])	# get cell with best value
		if (mlsc==1)	{
			mlsc <- 2
			} else if (mlsc==length(zg))	{
			mlsc <- length(zg)-1
			}
		improve <- min(abs(xx[mlsc,2]-xx[mlsc-1,2]),abs(xx[mlsc,2]-xx[mlsc+1,2]))
		lnl_ml_zm <- max(xx[,2])
		ml_zm_c <- match(lnl_ml_zm,xx[,2])
		ml_zm_g <- zg[ml_zm_c]
		} else	{
		improve <- 0
		}
	}
bH <- c(ml_zm_g,xx[ml_zm_c,1],hS,xx[ml_zm_c,2])
names(bH) <- c("Zipf_Mandel_log-log_decay","Zipf_Mandel_beta","Zipf_Mandel_hypothesized_richness","Zipf_Mandel_log-likelihood")
print(bH)
return(bH)
}

optimize_zipf_mandelbrot_abundance_given_hS_alt <- function(hS,counts,observed,init_zm,max_zb,max_zg)	{
# hS: hypothesized true richness
# counts: # finds per observed species
# observed: vector of length N giving # taxa with 1…N finds
# init_zg: log-log decay from best Zipf hypothesis
oS <- sum(observed)			# observed taxa
nspec <- sum((1:length(observed))*observed)
cl <- list(fnscale=-1)
zm <- init_zm
#b1 <- loglikelihood_zipf_rad(oS,nspec,zg=zm[1],hS=hS,observed)
xx <- optim(zm,fn=loglikelihood_zipf_mandelbrot_rad_both,method="L-BFGS-B",oS=oS,nspec=nspec,hS=hS,observed=observed,lower=c(0,0),upper=c(max_zg,max_zb),control=cl)
#if (xx$value < b1)	{
#	zm[2] <- 0.1
#	xx <- optim(zm,fn=loglikelihood_zipf_mandelbrot_rad_both,method="L-BFGS-B",oS=oS,nspec=nspec,hS=hS,observed=observed,lower=c(0,0),upper=c(2*max_zg,max_zb),control=cl)
#	}
#try <- 0
#while (xx$value < b1 && try<10)	{
#	zm[2] <- xx$par[2]/2
#	zm[1] <- xx$par[1]
#	xx <- optim(zm,fn=loglikelihood_zipf_mandelbrot_rad_both,method="L-BFGS-B",oS=oS,nspec=nspec,hS=hS,observed=observed,lower=c(0,0),upper=c(2*max_zg,max_zb),control=cl)
#	try <- try+1
#	}
bH <- c(xx$par[1],xx$par[2],hS,xx$val)
names(bH) <- c("Zipf_Mandel_log-log_decay","Zipf_Mandel_beta","Zipf_Mandel_hypothesized_richness","Zipf_Mandel_log-likelihood")
print(bH)
return(bH)
}

# get the most likely Zipf-Mandelbrot model
optimize_zipf_mandelbrot_abundance <- function(counts,init_zg=-1,init_hS=-1,span=5)	{
# counts: vector of length S with each element the number of finds for species 1…S
oS <- length(counts)				# observed taxa
if (init_hS==-1)	{
	stS <- oS
	}	else	{
	stS <- init_hS
	}
observed <- fisher_plot(counts)
if (init_zg==-1)	{
	mxfinds <- length(observed)
	mnfinds <- min((1:mxfinds)[!observed %in% 0])
	inzg <- log(mxfinds/mnfinds)/log(oS-1)
	} else	{
	inzg <- init_zg
	}

#zipf_dist <- zipf_distribution(zg=init_zg,S=stS)
inzb <- 50#(1/counts[1]/nspec)^(1/inzg)-1
inzb <- 0
init_zm <- c(inzg,inzb)
enS <- stS+((span-1)*oS)
incr <- floor(enS/span)
cl <- list(fnscale=-1)
peak <- 0
nspec <- sum(counts)
#max_zg <- 2*accio_max_zipf_decay(counts,oS,nspec)
max_zg <- max_zb <- 100
ispan <- span
pa <- 1
pz <- span
dhS <- c()
dresults <- c()
improve <- 1
while (incr>0 && improve>0.05)	{
	hS <- seq(stS,enS,by=incr)[!seq(stS,enS,by=incr) %in% dhS]
#	if (max_zg < 2*init_zm[1])	max_zg <- 2*init_zm[1]
	results <- sapply(hS,optimize_zipf_mandelbrot_abundance_given_hS_alt,counts=counts,observed=observed,init_zm=init_zm,max_zg=max_zg,max_zb=100)
#	print(hS)
#	results <- sapply(hS[pa:pz],optimize_zipf_mandelbrot_abundance_given_hS,counts=counts,observed=observed,max_zipf=max_zipf)
#	results <- sapply(hS,optimize_zipf_mandelbrot_abundance_given_hS,counts=counts,observed=observed,init_zg=init_zg)
	dhS <- c(dhS,hS)
	dresults <- cbind(dresults,results)
	if (length(dhS) > ispan)	{
		a <- sum(hS<bH[3]) 
		if (a==length(hS))	{
			results <- cbind(results,bH)
			} else if (a==0)	{
			results <- cbind(bH,results)
			} else	{
			b <- ncol(results)
			results <- cbind(results[,1:a],bH,results[,(a+1):b])
			}
		}
	hS <- results[3,]
	mlnl <- max(results[4,])
	mlSc <- match(mlnl,results[4,])
	mlS <- hS[mlSc]
	bH <- results[,mlSc]
	spn <- ncol(results)
	if (incr>1)	{
		# if runs so far are still producing higher likelihoods at higher richnesses:
		if (mlSc==spn)	{
			if (peak==0) {
				# if we have not yet found a peak, then keep increasing richness
				stS <- hS[spn]
				pa <- 2
				pz <- span
				enS <- stS+((spn-1)*oS)
				}	else	{
				enS <- hS[spn]
				if (span<incr)	{
					incr <- floor((enS-hS[spn-1])/span)
					}	else {
					span <- enS-hS[spn-1]	# check this! TDay
					incr <- 1
					}
				stS <- hS[spn]-(span*incr)
				pa <- 1
				pz <- span-1
				}	# end case where last number is best, but this is after finding a peak.
			} else if (mlSc==1)	{
			# if first richness is the best
			peak <- 1
			stS <- hS[mlSc]
			enS <- hS[mlSc+1]-1
			if (incr==1)	incr <- 0
			if (span<incr)	{
				incr <- floor((enS-stS)/span)
				}	else {
				span <- enS-stS
				incr <- 1
				}
			pa <- 2
			pz <- span
			} else {
			# if a richness in the middle is best & we still can find a better one
			peak <- 1
			hS2 <- c(mlS-1,mlS+1)
			init_zm <- results[1:2,mlSc]
			results2 <- sapply(hS2,optimize_zipf_mandelbrot_abundance_given_hS_alt,counts=counts,observed=observed,init_zm=init_zm,max_zb=max_zb,max_zg=100)
#			results2 <- sapply(hS2,optimize_zipf_mandelbrot_abundance_given_hS,counts=counts,observed=observed,init_zg=init_zg)
			dhS <- c(dhS,hS2)
			dresults <- cbind(dresults,results2)
			if (results2[4,1]>mlnl && results2[4,1]>results2[4,2])	{
			# start just above 2nd best richness and go up to best
				bH <- results2[,1]
				stS <- hS[mlSc-1]+1
				enS <- hS2[1]
				if (span<incr)	{
					incr <- max(1,floor((enS-stS)/span))
					if (incr==1)
						span <- 1+enS-stS
					}	else {
					span <- 1+enS-stS
					incr <- 1
					}
				pa <- 2
				pz <- span
				}	else if (results2[4,2]>mlnl && results2[4,2]>results2[4,1]) {
			# start at best richness and go just below 2nd best
				stS <- hS2[2]
				enS <- hS[mlSc+1]-1
				bH <- results2[,2]
				if (incr>span)	{
					incr <- floor((enS-stS)/span)
					}	else {
					span <- 1+enS-stS
					incr <- 1
					}
				pa <- 1
				pz <- span-1
				}	else {
				# we already had the best, so just end it
				incr <- 0
				}
			}
		# end case where we have a better richness in middle somewhere.  
		}	else	{
		incr <- 0
		}
	sup <- max(dresults[4,])-dresults[4,]
	improve <- min(sup[sup!=0])
	mlov <- match(max(dresults[4,]),dresults[4,])
	init_zm <- dresults[1:2,mlov]
	}
mlnl <- max(dresults[4,])
mlov <- match(max(dresults[4,]),dresults[4,])
if (dresults[2,mlov]>0)	{
	Zipf_Mandel_AICc <- modified_AIC(mlnl,k=3,n=nspec)
	} else	{
	Zipf_Mandel_AICc <- modified_AIC(mlnl,k=2,n=nspec)
	}
bH <- c(dresults[,mlov],Zipf_Mandel_AICc)
#names(bH)[2] <- "Zipf_Mandelbrot_Beta"
names(bH)[5] <- "Zipf_Mandelbrot_AICc"
return(bH)
}

# get gamma distribution
gamma_distribution <- function(a, b, S)	{
p0 <- vector(length=S)
p0[1] <- S/(S+1)
for(i in 2:S) p0[i] <- p0[i-1]-(1/(S+1))
prop <- qgamma(p0,a,b)/sum(qgamma(p0,a,b))
return(prop)
}

# get minimum alpha
accio_minimum_alpha_for_gamma_one <- function(hS)	{
min_alpha <- 1
iS <- hS
while (iS==hS)	{
	rel_ab_dist <- gamma_distribution(a=min_alpha,b=min_alpha,S=hS)		# basic lognormal distribution
	iS <- sum(rel_ab_dist>MINNO)
	if (iS == hS)	min_alpha <- min_alpha/2
	}
return(2*min_alpha)
}

# get log-likelihood for gamma
loglikelihood_gamma_one_rad <- function(oS,nspec,alpha,hS,observed)	{
#print(c(log10(alpha),hS))
rel_ab_dist <- gamma_distribution(a=alpha,b=alpha,S=hS)		# basic lognormal distribution
#print(rel_ab_dist)
if (sum(rel_ab_dist>MINNO) == hS)	{
	raw_expected <- expected_abundances(rel_ab_dist,nspec,S=hS)
	#lnl <- abundance_distribution_loglikelihood_suf(observed,expected=raw_expected,oS=oS,hS=hS)
	lnl <- distribution_loglikelihood_mul(observed,expected=raw_expected,oS,hS)
	}	else	lnl <- -1*MAXNO
return(round(lnl,2))
}

# get most likely magnitude of increase for a lognormal distribution of hS entities
optimize_gamma_one_abundance_given_hS <- function(observed,counts,hS)	{
#print(hS)		# for debugging
oS <- sum(observed)				# observed taxa
nspec <- sum((1:length(observed))*observed)
cl <- list(fnscale=-1)
alpha <- 2
min_alpha <- accio_minimum_alpha_for_gamma_one(hS)
w <- optim(alpha,loglikelihood_gamma_one_rad,method="L-BFGS-B",oS=oS,nspec=nspec,hS=hS,observed=observed,lower=min_alpha,upper=10000,control=cl)
#w <- optim(mag,fn=loglikelihood_lognormal_rad,method="L-BFGS-B",oS=oS,nspec=nspec,hS=hS,observed=observed,lower=mm[1],upper=mm[2],control=cl)
bH <- c(w$par,hS,w$value)
names(bH) <- c("Gamma_magnitude", "Gamma_S","Gamma_log-likelihood")
return(bH)
}

# get most likely gamma distribution for a population of some sort
optimize_gamma_one_abundance <- function(counts)	{
stS <- oS <- length(counts)				# observed taxa
observed <- fisher_plot(counts)
nspec <- sum(counts)
span <- 5
enS <- stS+((span-1)*oS)
incr <- floor(enS/span)
cl <- list(fnscale=-1)
peak <- 0
while (incr>0)	{
	hS <- seq(stS,enS,by=incr)
	results <- sapply(hS,optimize_gamma_one_abundance_given_hS,observed=observed,counts=counts)
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
				stS <- hS[spn]
				enS <- stS+((spn-1)*oS)
				}	else	{
				enS <- hS[spn]
				if (span<incr)	{
					incr <- floor((enS-hS[spn-1])/span)
					}	else {
					span <- enS-hS[spn-1]
					incr <- 1
					}
				stS <- hS[spn]-(span*incr)
				}	# end case where last number is best, but this is after finding a peak.
			} else if (mlSc==1)	{
			# if first richness is the best
			peak <- 1
			stS <- hS[mlSc]
			enS <- hS[mlSc+1]-1
			if (incr==1)	incr <- 0
			if (span<incr)	{
				incr <- floor((enS-stS)/span)
				}	else {
				span <- enS-stS
				incr <- 1
				}
			} else {
			# if a richness in the middle is best & we still can find a better one
			peak <- 1
			hS2 <- c(mlS-1,mlS+1)
			results2 <- sapply(hS2,optimize_lognormal_abundance_given_hS,observed=observed,counts=counts)
			if (results2[3,1]>mlnl && results2[3,1]>results2[3,2])	{
			# start just above 2nd best richness and go up to best
				stS <- hS[mlSc-1]+1
				enS <- hS[mlSc]
				if (span<incr)	{
					incr <- floor((enS-stS)/span)
					}	else {
					span <- enS-stS
					incr <- 1
					}
				}	else if (results2[3,2]>mlnl && results2[3,2]>results2[3,1]) {
			# start at best richness and go just below 2nd best
				stS <- hS[mlSc]
				enS <- hS[mlSc+1]-1
				if (incr>span)	{
					incr <- floor((enS-stS)/span)
					}	else {
					span <- enS-stS
					incr <- 1
					}
				}	else {
				# we already had the best, so just end it
				incr <- 0
				}
			}
		# end case where we have a better richness in middle somewhere.  
		}	else	{
		incr <- 0
		}
	}
bH <- c(round(bH[1],6),bH[2],bH[3],modified_AIC(mlnl,2,nspec))
names(bH)[4] <- "Gamma_AICc"
return(bH)
}

# get the expected number of taxa with 1…J individuals given theta & m (new species appearance rates)
#	This is modified from the untb volkov function, which chokes on both very low and very high numbers
#	My "kluge" solution is to replace numbers rounded to zero or infinity with the minimum or maximum numbers that R can handle
#	Note that I cut all of the binning stuff: you can do that on your own with minimal effort
Volkov_distribution <- function (J, theta, m) {
gam <- m * (J - 1)/(1 - m)

volkov_integrand <- function(y, n) {
#	print(y)			# for debugging
	log_mult <- (lgamma(J + 1) - lgamma(n + 1) - lgamma(J -n + 1) + lgamma(gam) - lgamma(J + gam) + lgamma(n + y) + lgamma(J - n + gam - y) - lgamma(1 + y) - lgamma(gam - y) - y * theta/gam)
	log_mult[log_mult<log(MINNO)] <- log(MINNO)
	accio <- theta * exp(log_mult)
	accio[accio==Inf] <- MAXNO
	1*accio
    }

integrate_volkov_for_n <- function(n) {
#	print(n)			# for debugging
	integrate(volkov_integrand, lower = 0, upper = gam, n = n, stop.on.error=FALSE)$value
    }

nn <- (1:J)	# vector of possible numbers of individuals
out <- sapply((1:J), integrate_volkov_for_n)
return(out)
}

# calculate Zero Sum Multinomial for total population J, new species rate m and theta
zero_sum_multinomial_distribution <- function(J,theta,m)	{
# J is the true (incompletely sampled) population size
#	NOTE: We do not want to use counted specimens for J
#	Low J assumes a "new" community for which immigration is still important!
#	Particularly with time averaging, we are sampling "old" communities in part
# theta = Hubbell's Fundamental Biodiversity Number
# m = migration rate
# Volkov_distribution function is modified heavily from untb
init_abund <- Volkov_distribution(J, theta, m)	# this gives probability of initial abundance
init_abund[is.na(init_abund)] <- MINNO
cum_init_abund <- round(cumsum(init_abund),0)
sp_no_up <- unique(cum_init_abund)
abunds <- match(sp_no_up,cum_init_abund)
unq_fnds <- length(sp_no_up)
S <- max(sp_no_up)
sp_no_up <- c(0,sp_no_up)
ab_dist <- vector(length=S)
for (i in 1:unq_fnds)	ab_dist[(sp_no_up[i]+1):sp_no_up[i+1]] <- abunds[i]
return (rel_ab_dist=sort(ab_dist/sum(ab_dist),decreasing=TRUE))
}

# find the minimum theta that will 
accio_minimum_theta_and_maximum_m_given_J <- function(J,S,max_m=0.999,dbug=FALSE)	{
last_theta <- min_theta <- 1
ad <- Volkov_distribution(J,theta=min_theta,m=max_m)
new_S <- sum(ad)
while (new_S < 1)	{
	last_S <- new_S
	last_m <- max_m
	max_m <- 1 - ((1-max_m)*2)
	ad <- Volkov_distribution(J,theta=min_theta,m=max_m)
	new_S <- sum(ad)
	if (dbug)	print(paste("m =",round(max_m,5),"S =",round(new_S,2),date(),sep=" "))
	}

if (new_S < S)	{
	while(new_S < S)	{
		last_theta <- min_theta
		min_theta <- min_theta*1.5
		ad <- Volkov_distribution(J,theta=min_theta,m=max_m)
		last_S <- new_S
		new_S <- sum(ad)
		if (dbug)	print(paste("θ =",round(min_theta,2),"S =",round(new_S,2),date(),sep=" "))
		}
	bst_theta <- last_theta + (min_theta - last_theta)*(1-((new_S-S)/(new_S-last_S)))
	}	else if (new_S >= S)	{
	while(new_S >= S)	{
		last_theta <- min_theta
		min_theta <- min_theta/1.1
		ad <- Volkov_distribution(J,theta=min_theta,m=max_m)
		last_S <- new_S
		new_S <- sum(ad)
		if (dbug)	print(paste("θ =",round(min_theta,2),"S =",round(new_S,2),date(),sep=" "))
		}
	bst_theta <- last_theta - (last_theta - min_theta)*(1-((S-new_S)/(last_S-new_S)))
	}
baseline <- c(bst_theta,max_m)
names(baseline) <- c("min_theta","max_m")
return(baseline)
}

# get the minimum plausible m (immigration rate) given some J & theta that yields S species
accio_minimum_m_and_maximum_theta_given_J <- function(J,S,dbug=FALSE)	{
max_theta <- 100*log10(J)
last_m <- next_m <- 0.5
ad <- Volkov_distribution(J,max_theta,m=next_m)
max_S <- new_S <- round(sum(ad),0)
if (new_S < S)	{
	# if the distributions have too few species, increase m
#	ad_l <- Volkov_distribution(J,max_theta,m=(next_m-0.05))
#	new_Sl <- round(sum(ad_l),0)
#	ad_u <- Volkov_distribution(J,max_theta,m=(next_m+0.05))
#	new_Su <- round(sum(ad_u),0)
#	new_S <- max(new_Sl,new_Su)
#	if (new_Sl>new_Su)	dm <- -0.05
	while (round(new_S,0) < 1)	{
		last_theta <- max_theta
		max_theta <- max_theta / 1.25
		ad <- Volkov_distribution(J,max_theta,m=next_m)
		last_S <- new_S
		new_S <- round(sum(ad),0)
		if (dbug)	print(c(round(max_theta,3),round(new_S,2),date()))
		}
	ad_l <- Volkov_distribution(J,max_theta,m=(next_m-0.05))
	new_Sl <- sum(ad_l)
	ad_u <- Volkov_distribution(J,max_theta,m=(next_m+0.05))
	new_Su <- sum(ad_u)
	if (new_Sl > new_Su)	{
		up <- FALSE
		}	else up <- TRUE
	m_mod <- 1.1
	theta_mod <- 1.05
	last_theta <- next_theta <- max_theta
	overshot <- 0
	while (new_S < S)	{
		last_S <- new_S
		if (overshot==0)	{
			penult_m <- last_m
			last_m <- next_m
			if (up)	{
				next_m <- 1 - ((1-last_m)/m_mod)
				}	else {
				next_m <- last_m/m_mod
				}
			}	else	{
			last_theta <- next_theta
			next_theta <- last_theta * theta_mod
			}
		ad <- Volkov_distribution(J,next_theta,m=next_m)
#		adl <- Volkov_distribution(J,max_theta/1.1,m=next_m)
#		adu <- Volkov_distribution(J,max_theta*1.1,m=next_m)
		new_S <- sum(ad)
#		new_Sl <- sum(adl)
#		new_Su <- sum(adu)
		if (dbug)	print(c(round(next_m,4),round(next_theta,2),round(new_S,3),date()))
		if (new_S > max_S)	{
			max_S <- new_S
			}	else	{
			overshot <- 1
			next_m <- last_m
			new_S <- max_S
			}
		}
	min_m <- next_m+(last_m-next_m)*((new_S-S)/(new_S-last_S))
	max_theta <- last_theta+(next_theta-last_theta)*((new_S-S)/(new_S-last_S))
	}	else {
	while (new_S > S && next_m > MINNO)	{
		last_S <- new_S
		last_m <- next_m
		next_m <- last_m/1.1
		ad <- Volkov_distribution(J,max_theta,m=next_m)
		new_S <- sum(ad)
		if (dbug)
			print(paste(round(next_m,4),round(next_theta,2),round(new_S,3),date(),sep=" "))
#		if (dbug)	print(c(next_m,new_S,date()))
		}
	min_m <- next_m-(next_m-last_m)*((S-new_S)/(last_S-new_S))
	}
baseline <- c(max_theta,min_m)
names(baseline) <- c("max_theta","min_m")
return(baseline)
}

# rough approximation of theta given population size and observed richness
accio_init_theta <- function(counts)	{
oS <- length(counts)
nspec <- sum(counts)
theta <- 1
S <- 1 + (theta * log(1+(nspec-1)/theta))
if (S > oS)	{
	while (S > oS)	{
		theta <- theta/1.05
		S <- 1 + (theta * log(1+(nspec-1)/theta))
		}
	}	else if (S < oS)	{
	while (S < oS)	{
		theta <- theta*1.05
		S <- 1 + (theta * log(1+(nspec-1)/theta))
		}
	}
return (theta)
}

# calculate Zero Sum RAD for a theta & m given J
loglikelihood_zero_sum_multinomial_rad <- function(hub,min_theta,max_theta,J,nspec,observed,oS,dbug=FALSE)	{
# hub: includes theta (hub[1]) & m (hub[2])
# J: True Population Size
# nspec: number pf specimens
# observed: a linear Fisher plot giving # species with X specimens
# oS observed species
theta <- min_theta+(hub[1]*(max_theta-min_theta))
m <- hub[2]
#if (dbug==2)	print(c(round(theta,8),round(m,8)))				# for debugging
rel_ab_dist <- zero_sum_multinomial_distribution(J=J,theta=theta,m=m)		# get expected RAD given zero sum parameters
hS <- length(rel_ab_dist)
#if (hS>=oS)	{
if (length(rel_ab_dist)>=oS)	{
	raw_expected <- expected_abundances(rel_ab_dist,nspec,S=hS)
	lnl <- distribution_loglikelihood_mul(observed,expected=raw_expected,oS=oS,hS=hS)
	}	else {
	lnl <- oS*log(MINNO)
	}
if (dbug)	print(paste(round(theta,3),round(m,3),round(lnl,2),J,sep=" "))
return(lnl)
}

# get the best zero sum multinomial at a given J
optimize_zero_sum_multinomial_abundance_given_J <- function(J,observed,init_theta,init_m,apprise=FALSE)	{
oS <- sum(observed)				# observed taxa
nspec <- sum((1:length(observed))*observed)
#if (apprise)	print(paste("Zero Sum Multinomial J = ",J,sep=""))		# for updating, as zero-sum is slow
if (apprise)	print(paste("Zero Sum Multinomial J = ",J," ",date(),sep=""))		# for updating, as logseries is slow!
# given J and the minimum theta, there is some minimum m that generates at least oS taxa
#vvv <- accio_minimum_m_and_maximum_theta_given_J(J,S=oS,dbug=TRUE)
#dd <- fitvolkov(counts)$Coefficients
min_theta <- min(init_theta/4,1)
max_theta <- 4*init_theta
hub <- c((init_theta-min_theta)/(max_theta-min_theta),init_m)
cl <- list(fnscale=-1)
w <- optim(hub,fn=loglikelihood_zero_sum_multinomial_rad,method="L-BFGS-B",min_theta=min_theta,max_theta=max_theta,J=J,oS=oS,nspec=nspec,observed=observed,dbug=FALSE,lower=c(0,1e-05),upper=c(1,0.99999),control=cl)
bH <- c(w$par,J,w$value)
bH[1] <- min_theta+(w$par[1]*(max_theta-min_theta))
names(bH) <- c("Zero_Sum_Multinomial_Theta","Zero_Sum_Multinomial_m","Zero_Sum_Multinomial_J","Zero_Sum_Multinomial_log-likelihood")
return(bH)
} # end

accio_zero_sum_multinomial_init_parameters <- function(counts)	{
oS <- length(counts)
init_theta <- accio_init_theta(counts)
init_m <- optimal.prob(counts)
iJ <- sum(counts)
rel_ab_dist <- zero_sum_multinomial_distribution(J=iJ,theta=init_theta,m=init_m)		# get expected RAD given zero sum parameters
iS <- length(rel_ab_dist)
# make sure that parameters create enough richness given observed specimens
while (iS < oS)	{
	init_theta <- init_theta * 1.01
	init_m <- 1-((1-init_m)/1.01)
	rel_ab_dist <- zero_sum_multinomial_distribution(J=iJ,theta=init_theta,m=init_m)		# get expected RAD given zero sum parameters
	iS <- length(rel_ab_dist)
#	print(c(init_theta,init_m,iS))
	}
return(c(init_theta,init_m))
}

# get the best zero sum multinomial
optimize_zero_sum_multinomial_abundance <- function(counts,byte_trail=FALSE,span=4,max_J=MAX_J,apprise=FALSE)	{
oS <- length(counts)
observed <- fisher_plot(counts)
base <- 10
peak <- 0
nspec <- sum(counts)
exp_st <- log10(nspec)
incr <- 0.5
powers <- exp_st+incr*(0:span)
exp_en <- max(powers)
pa <- 1
pz <- length(powers)
if (byte_trail)	byte_trail_file <- paste("Zero_Sum_Trail_oS=",oS,"_N=",nspec,".txt",sep="")
init_params <- accio_zero_sum_multinomial_init_parameters(counts)
init_theta <- init_params[1]
init_m <- init_params[2]
while (incr>0.001) {
	J <- round(base^powers)
	results <- sapply(J[pa:pz],optimize_zero_sum_multinomial_abundance_given_J,observed=observed,init_theta=init_theta,init_m=init_m,apprise=apprise)
	if (byte_trail)	{
		if (pa==1)	{
			byte_trail_results <- t(results)
			}	else {
			byte_trail_results <- rbind(byte_trail_results,t(results))
			}
		write.table(byte_trail_results,byte_trail_file,sep="\t",col.names=TRUE,row.names=FALSE)
		}
	if (pa==2)
		results <-cbind(old_result_1,results)
	if (pz<length(powers))
		results <-cbind(results,old_result_2)
	mlln <- max(results[4,])
	mlJc <- match(mlln,results[4,])
	mlJ <- results[3,mlJc]
	bH <- results[,mlJc]
	if (mlJc==length(results[4,]))	{
		if (peak==0)	{
			lower <- powers[mlJc]
			upper <- min(log10(max_J),(lower + (incr*span)))
			incr <- (upper-lower)/span
			powers <- lower+incr*(0:span)
			old_result_1 <- results[,mlJc]
			pa <- 2
			# we are still improving by increasing J
			}	else	{
			# already found a peak J, but the higest in this search is the best
			lower <- powers[mlJc-1]
			upper <- powers[mlJc]
			incr <- (upper-lower)/span
			old_result_1 <- results[,mlJc-1]
			old_result_2 <- results[,mlJc]
			pa <- 2
			pz <- length(powers) - 1
			}
		powers <- lower+incr*(0:span)
		# end case where highest J is the best
		}	else if (mlJc==1)	{
			# the lowest examined J is best
		if (peak==0)	{
			upper <- powers[mlJc+1]
			lower <- powers[mlJc]
			old_result_1 <- results[,mlJc]
			old_result_2 <- results[,mlJc+1]
			pa <- 2
			pz <- length(powers) - 1
			incr <- (upper-lower)/span
			powers <- lower+incr*(0:span)
			peak <- 1
			}	else {
			peak <- 1
			lower <- powers[mlJc]
			upper <- powers[mlJc+1]
			incr <- (upper-lower)/span
			powers <- lower+incr*(0:span)
			old_result_1 <- results[,mlJc]
			old_result_2 <- results[,mlJc+1]
			pa <- 2
			pz <- length(powers) - 1
			}
		# end case where lowest J is the best
		}	else	{
		# some J in the middle is the best
		peak <- 1
		J2 <- round(c(base^(powers[mlJc]-(incr/10)),base^(powers[mlJc]+(incr/10))),0)
		results2 <- sapply(J2,optimize_zero_sum_multinomial_abundance_given_J,observed=observed,init_theta=init_theta,init_m=init_m,apprise=apprise)
#		results2 <- sapply(J2,optimize_zero_sum_multinomial_abundance_given_J,observed=observed,max_m=0.999)
		if (byte_trail)	{
			byte_trail_results <- rbind(byte_trail_results,t(results2))
			write.table(byte_trail_results,byte_trail_file,sep="\t",col.names=TRUE,row.names=FALSE)
			}
		if (results2[4,1]>results2[4,2])	{
			# lower J improves likelihood
			lower <- powers[mlJc-1]
			old_result_1 <- results2[,1]
			if (results2[4,1]>results[4,mlJc])	{
				# if lower J is better than bH, then reset
				upper <- powers[mlJc]-(incr/10)
				bH <- old_result_2 <- results2[,1]	# the new most likely hypothesis so far
				}	else {
				upper <- powers[mlJc]
				old_result_2 <- results[,mlJc]		# no change in the most likely hypothesis
				}
			incr <- (upper-lower)/span
			powers <- lower+incr*(0:span)	# HERE
			} else	{
			# higher J improves likelihood
			upper <- powers[mlJc+1]
			old_result_2 <- results[,mlJc+1]
			if (results2[4,2]>results[4,mlJc])	{
				# if higher J is better than bH, then reset
				lower <- powers[mlJc]+(incr/10)
				bH <- old_result_1 <- results2[,2]	# the new most likely hypothesis so far
				}	else {
				lower <- powers[mlJc]
				old_result_1 <- results[,mlJc]		# no change in the most likely hypothesis
				}
			incr <- (upper-lower)/span
			powers <- lower+incr*(0:span)
			}
		pa <- 2
		pz <- length(powers)-1
		}
	}
Zero_Sum_Multinomial_AICc <- modified_AIC(bH[4],k=3,n=nspec)
names(Zero_Sum_Multinomial_AICc) <- "Zero_Sum_Multinomial_AICc"
bH <- c(bH,Zero_Sum_Multinomial_AICc)
return(bH)
}

# get best zero sum multinomial for all three variables
optimo_zero_sum_multinomial_abundance <- function(hubbells,observed,min_theta,max_theta,min_m,max_m,min_J,max_J)	{
test_theta <- min_theta+(hubbells[1]*(max_theta-min_theta))
test_m <- min_m+(hubbells[2]*(max_m-min_m))
test_J <- round((min_J+(hubbells[3]*(max_J-min_J))),0)
rel_ab_dist <- zero_sum_multinomial_distribution(test_J,test_theta,test_m)
hS <- length(rel_ab_dist)
if (hS>=oS)	{
	raw_expected <- expected_abundances(rel_ab_dist,nspec,S=hS)
	lnl <- distribution_loglikelihood_mul(observed,expected=raw_expected,oS=oS,hS=hS)
	}	else {
	lnl <- oS*log(MINNO)
	}
return(lnl)
}

# calculate the likelihood and AICc of a model predicting exactly the numbers of taxa with 1…N finds
#	Note: this ad hoc model has one parameter for each taxon
saturated_abundance_model <- function(counts)	{
observed <- fisher_plot(counts)
nspec <- sum(counts)
oS <- sum(observed)
prop_expected <- observed/sum(observed)
sat_lnl <- sum(observed[observed>0]*log(prop_expected[prop_expected>0]))
sat_AICc <- modified_AIC(sat_lnl,oS,n=nspec)
bH <- c(round(sat_lnl,3),round(sat_AICc,2))
names(bH) <- c("Saturated_log-likelihood","Saturated_AICc")
return(bH)
}

# a parametric bootstrap test to find the range of expected best-possible-fit
#	likelihoods given a true model distribution
parametric_bootstrap_test <- function(runner_up,params,oS,nspec,runs=1000)	{
if (runner_up=="Geometric")	{
	trad <- geometric_distribution(params[1])
	}	else if (runner_up=="Log_Series")	{
	trad <- logseries_distribution(params[1],params[2])
	}	else if (runner_up=="Zero_Sum_Multinomial")	{
	trad <- zero_sum_multinomial_distribution(params[1],params[2],params[3])
	}	else if (runner_up=="Lognormal")	{
	trad <- lognormal_distribution(params[1],params[2])
	}	else if (runner_up=="Zipf")	{
	trad <- zipf_distribution(params[1],params[2])
	}

for (i in 1:runs)	{
	if (i==1)	{
		sim_mps <- accio_expected_best_possible_distribution_support(trad,oS,nspec)
		}	else {
		sim_mps <- rbind(sim_mps,accio_expected_best_possible_distribution_support(trad,oS,nspec))
		}
	}
return(sim_mps)
}

# simulate a range of communities following an model distribution to get expected "best-possible" for oS taxa
accio_expected_best_possible_distribution_support <- function(trad,oS,nspec)	{
sS <- 0
overshot <- 0
while (sS != oS)	{
	sim_sad <- simulated_community(trad,nspec)
	sS <- length(sim_sad)
	if (sS > oS)	overshot <- overshot+1
	if (overshot==100)	sS <- oS
#	print(sS)
	}
sim_sat <- saturated_abundance_model(sim_sad[1:oS])
return(sim_sat)
}
	
# simulate a communities following an model distribution
simulated_community <- function(trad,nspec)	{
samples <- sort(runif(nspec))
crad <- cumsum(trad)
srad <- count_simulated_specimens(crad,samples,nspec)
return(srad)
}

# count simulated specimens given a string of uniform random numbers
count_simulated_specimens <- function (crad,samples,nspec)	{
sS <- 1
sspec <- 0 
while (sspec < nspec)	{
	if (sS==1)	{
		srad <- length(subset(samples,samples<crad[sS]))
		}	else {
		srad <- c(srad,length(subset(samples,samples<crad[sS]))-sspec)
		}
	sspec <- sum(srad)
	sS <- sS+1
	}
return(srad)
}

# routine to test requested distributions, set up to use sapply or lapply
test_abundance_distributions <- function(asm,abundance_data,tests,assemblages,assemblage_names,temporary_output=TRUE,byte_trail=FALSE,span=4,max_J=MAX_J,apprise=FALSE,provide_saturated=TRUE,analysis_name="")	{
# get data for this assemblage
counts <- sort(abundance_data$Counts[abundance_data$Assemblage==asm],decreasing=TRUE);
collection_name <- assemblage_names[match(asm,assemblages)];
summary_stats <- c(asm,as.character(collection_name),length(counts),sum(counts))
names(summary_stats) <- c("Assemblage_No","Assemblage_Name","Taxa","Individuals")
ct <- 0
asmno <- assemblages[match(asm,assemblages)]
if (!is.na(match("Geometric",tests)))	{
	status <- paste("Finding best Geometric for",assemblage_names[asmno],date(),sep=" ")
	print(status)
	bst_geom <- optimize_geometric_abundance(counts)
	summary_stats <- c(summary_stats,bst_geom)
	ct <- ct+1
	if (ct==1)	{
		AICcs <- bst_geom[3]
		}	else AICcs <- c(AICcs,bst_geom[3])
	if (bst_geom[3]==min(AICcs))	{
		if (ct>1)	runner_up <- winner
		winner <- "Geometric"
		best_params <- bst_geom[1:2]
		}
	}
if (!is.na(match("Logseries",tests)))	{
	status <- paste("Finding best Logseries for",assemblage_names[asmno],date(),sep=" ")
	print(status)
	bst_lgsr <- optimize_logseries_abundance(counts,byte_trail=byte_trail,span=span,max_J=max_J,apprise=apprise)
	summary_stats <- c(summary_stats,bst_lgsr)
	ct <- ct+1
	if (ct==1)	{
		AICcs <- bst_lgsr[4]
		}	else AICcs <- c(AICcs,bst_lgsr[4])
	if (bst_lgsr[4]==min(AICcs))	{
		if (ct>1)	runner_up <- winner
		winner <- "Log_Series"
		best_params <- bst_lgsr[1:2]
		}
	}
if (!is.na(match("Zero-Sum-Multinomial",tests)))	{
	status <- paste("Finding best Zero Sum Multinomial for",assemblage_names[asmno],date(),sep=" ")
	print(status)
	bst_zero <- optimize_zero_sum_multinomial_abundance(counts,byte_trail=byte_trail,span=span,max_J=max_J,apprise=apprise)
	summary_stats <- c(summary_stats,bst_zero)
	ct <- ct+1
	if (ct==1)	{
		AICcs <- bst_zero[5]
		}	else AICcs <- c(AICcs,bst_zero[5])
	if (bst_zero[5]==min(AICcs))	{
		if (ct>1)	runner_up <- winner
		winner <- "Zero_Sum_Multinomial"
		best_params <- bst_zero[1:3]
		}
	}
if (!is.na(match("Lognormal",tests)))	{
	status <- paste("Finding best Lognormal for",assemblage_names[asmno],date(),sep=" ")
	print(status)
	bst_lgnr <- optimize_lognormal_abundance(counts)
	summary_stats <- c(summary_stats,bst_lgnr)
	ct <- ct+1
	if (ct==1)	{
		AICcs <- bst_lgnr[4]
		}	else AICcs <- c(AICcs,bst_lgnr[4])
	if (bst_lgnr[4]==min(AICcs))	{
		if (ct>1)	runner_up <- winner
		winner <- "Lognormal"
		best_params <- bst_lgnr[1:3]
		}
	}
if (!is.na(match("Zipf",tests)))	{
	status <- paste("Finding best Zipf for",assemblage_names[asmno],date(),sep=" ")
	print(status)
	bst_zipf <- optimize_zipf_abundance(counts)
	summary_stats <- c(summary_stats,bst_zipf)
	ct <- ct+1
	if (ct==1)	{
		AICcs <- bst_zipf[4]
		}	else AICcs <- c(AICcs,bst_zipf[4])
	if (bst_zipf[4]==min(AICcs))	{
		if (ct>1)	runner_up <- winner
		winner <- "Zipf"
		best_params <- bst_zipf[1:3]
		}
	}
# get best possible likelihood (i.e., predicts exact distribution with N parameters = No. of taxa)
if (provide_saturated)	{
	bst_sat <- saturated_abundance_model(counts)
	ct <- ct+1
	if (ct==1)	{
		AICcs <- bst_sat[2]
		}	else AICcs <- c(AICcs,bst_sat[2])
	if (bst_sat[2]==min(AICcs))	{
		if (ct>1)	runner_up <- winner
		winner <- "Saturated"
		}
	summary_stats <- c(summary_stats,bst_sat)
	}
names(winner) <- "Best_Model"
summary_stats <- c(summary_stats,winner)
if (temporary_output)	{
	if (analysis_name=="") {
		temp_file_name <- paste("Collection_",asm,"_",assemblage_names[asmno],".xls",sep="")
		} else {
		temp_file_name <- paste(analysis_name,"_Collection_",asm,"_",assemblage_names[asmno],".xls",sep="")
		}
	write.table(summary_stats,temp_file_name,sep="\t",col.names=TRUE,row.names=TRUE)
	}
return(summary_stats)
}

# Main program using functions above.
# 	occurrences_file specifies tab delimited text file giving collection # and taxon #
# 	collections_file specifies tab delimited text file giving collection # and sampling bin (e.g., stratigraphic stages, geographic units, geographic units by stages, etc.)
# 	tests specify distributions to test
# 	min_taxa minimum number of taxa to test an assemblage
# 	min_counts minimum number of finds to test an assemblage
#	temporary_output: if TRUE, then provide output file for each collection analyzed
#	byte_trail: if TRUE, then most-likely results for particular hypothesized J values for logseries and zero sum
#	maxJ: maximum hypothesized population size for logseries and zero-sum analyses
#	apprise: if TRUE, then provide apprises for slow-moving logseries and zero sum analyses
#	provide_saturated: if TRUE, then get the best possible model. (Really only relevant for N>20 or so)
#	analysis_name: if provided, then append output files with names
find_best_relative_abundance_models <- function(abundances_file,assemblage_names_file,tests,min_counts,temporary_output=FALSE,byte_trail=FALSE,span=4,max_J=MAX_J,apprise=FALSE,provide_saturated=FALSE,analysis_name)	{
abundance_data <- read.table(file=abundances_file, header=TRUE, stringsAsFactors=TRUE,sep="\t")
keep <- (1:ncol(abundance_data))[colnames(abundance_data) %in% c("Assemblage","Taxon","Counts")];
abundance_data <- abundance_data[,keep];
#colnames(abundance_data) <- c("Assemblage","Taxon","Counts")
abundance_data <- subset(abundance_data,abundance_data$Counts>0);
assemblages <- sort(unique(abundance_data$Assemblage));
taxa <- sort(unique(abundance_data$Taxon));

assemblage_finds <- finds_per_assemblage(abundance_data);
assemblage_rich <- taxa_per_assemblage(abundance_data);
assemblage_names <- simplify2array(read.table(file=assemblage_names_file, header=FALSE, stringsAsFactors=TRUE, sep="\t"));

# get assemblages with enough specimens
asm1 <- assemblages[assemblage_finds>=min_counts]
# get assemblages with enough taxa
asm2 <- assemblages[assemblage_rich>=min_taxa]
# get intersection of taxon-rich and high-count assemblages
asm <- asm1[asm1 %in% asm2]
#asm_names <- assemblage_names[asm1 %in% asm2]
output_summary <- sapply(asm,test_abundance_distributions,abundance_data=abundance_data,tests=tests,assemblages=assemblages,assemblage_names=assemblage_names,temporary_output=temporary_output,byte_trail=byte_trail,span=span,max_J=MAX_J,apprise=apprise,provide_saturated)

return(output_summary)
}

find_best_relative_abundance_models_by_group <- function(abundances_file,taxon_file,individual_groups,assemblage_names_file,tests,min_counts,temporary_output=FALSE,byte_trail=FALSE,span=4,max_J=MAX_J,apprise=FALSE,provide_saturated=FALSE)	{
abundance_data <- read.table(file=abundances_file, header=TRUE, stringsAsFactors=FALSE)
colnames(abundance_data) <- c("Assemblage","Taxon","Counts")
abundance_data <- subset(abundance_data,abundance_data$Counts>0)
assemblages <- sort(unique(abundance_data$Assemblage))
taxa <- sort(unique(abundance_data$Taxon))

assemblage_finds <- finds_per_assemblage(abundance_data)
assemblage_rich <- taxa_per_assemblage(abundance_data)
assemblage_names <- simplify2array(read.table(file=assemblage_names_file, header=FALSE, stringsAsFactors=TRUE, sep="\t"))

# get assemblages with enough specimens
asm1 <- assemblages[assemblage_finds>=min_counts]
# get assemblages with enough taxa
asm2 <- assemblages[assemblage_rich>=min_taxa]
# get intersection of taxon-rich and high-count assemblages
asm <- asm1[asm1 %in% asm2]

taxa <- read.table(file=taxon_file, header=TRUE, stringsAsFactors=TRUE,sep="\t")
ttl_taxa <- nrow(taxa)
groups <- length(individual_groups)
ttl_finds <- nrow(abundance_data)

reduced_abundances <- abundance_data[(1:ttl_finds)[abundance_data$Assemblage %in% asm],]
output_summary <- c()
for (tx in 1:groups)	{
	ingroup <- (1:ttl_taxa)[taxa$Laflamme_clade %in% individual_groups[tx]]
	ingroup_abundances <- reduced_abundances[reduced_abundances$Taxon %in% ingroup,]
	ingroup_assemblages <- unique(ingroup_abundances$Assemblage)
	ingroup_assemblage_finds <- finds_per_assemblage(ingroup_abundances)
	ingroup_assemblage_rich <- taxa_per_assemblage(ingroup_abundances)

	# get assemblages with enough specimens
	iasm1 <- ingroup_assemblages[ingroup_assemblage_finds>=min_counts]
	# get assemblages with enough taxa
	iasm2 <- ingroup_assemblages[ingroup_assemblage_rich>=min_taxa]
	# get intersection of taxon-rich and high-count assemblages
	asm <- iasm1[iasm1 %in% iasm2]
	ioutput_summary <- sapply(asm,test_abundance_distributions,abundance_data=ingroup_abundances,tests=tests,assemblages=assemblages,assemblage_names=ingroup_assemblages,temporary_output=temporary_output,byte_trail=byte_trail,span=span,max_J=MAX_J,apprise=apprise,provide_saturated=provide_saturated)
	output_summary <- rbind(output_summary,base::t(ioutput_summary))
	}

return(output_summary)
}
#file_name <- paste(individual_groups[tx],"_Output.xls",sep="")
#write.table(output_summary,file_name,row.names = FALSE,col.names = TRUE,sep="\t")