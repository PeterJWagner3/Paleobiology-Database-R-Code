# some basic stats & probability routines that I wrote for myself because I didn't like the way that R did them...

mann_whitney <- function(data_v,category)	{
ranked <- rank(data_v)
cats <- unique(category)
n1 <- sum(category==cats[1])
n2 <- sum(category==cats[2])
nt <- length(category)
R1 <- sum(ranked[(1:nt)[category==cats[1]]])
U1 <- R1 -((n1*(n1+1))/2)
R2<- sum(ranked[(1:nt)[category==cats[2]]])
U2 <- R2 -((n2*(n2+1))/2)
mU <- (n1*n2)/2
unique_ranks <- sort(unique(ranked))
numr <- 0
for (i in 1:length(unique_ranks))	{
	ti <- sum(ranked %in% unique_ranks[i])
	numr <- numr + ((ti^3)-ti)
	}
numr <- numr/(nt*(nt-1))
sigU <- ((n1*n2/12)*(nt+1-numr))^0.5	# expected standard deviation in ranks corrected for ties
#(n1*n2*(n1+n2+1)/12)^0.5
if (abs(U2-mU)>0)	{
	pval <- 2*integrate(dnorm,lower=abs(U2-mU)/sigU,upper=10)$value;
	} else	{
	pval <- 1.0;  # if we have only the same numbers in both categories, then it's a perfect match
	}
output <- data.frame(U=as.numeric(max(U2,U1)),pval=as.numeric(pval),stringsAsFactors = F);
return(output)
}

Poisson_rate_to_probability <- function(expectation)	{
return(1-exp(-expectation))
}

probability_to_Poisson_rate <- function(pn)	{
return(-1*log(1-pn))
}

binomial_support_bars_ext <- function (n,N,support=1.0)	{
ml_Rate <- n/N;
rates <- seq(0.001,0.999,by=0.001);
rate_loglikelihood <- log(dbinom(n,N,rates));
max_lnL <- max(rate_loglikelihood);
rate_support <- rate_loglikelihood-max_lnL;
ok_rates <- rates[rate_support>=-support];
return(c(min(ok_rates),max(ok_rates)));
}

colMax <- function(data_matrix)	{
col_maxes <- c();
for (i in 1:ncol(data_matrix))	col_maxes <- c(col_maxes,max(data_matrix[,i]));
return(col_maxes)
}
	
