# generates a jackknifed distribution of the statistic given in theta.fun
# theta.fun is the function applied to each pseudosample
# x is either a vector of univariate data or a vector of row indices for multivariate data
# xdata is the matrix of multivariate data that matches x if x is row indices
# jack.gf is grouping factor in case the jackknifed units are not the individual
# (e.g., family, locality, etc.). Default is to jackknife by individual
# ... additional arguments transfered to theta.fun
# CI calculation is from 
# http://www.math.ntu.edu.tw/~hchen/teaching/LargeSample/references/R-bootstrap.pdf 
# and http://darwin.phyloviz.net/ComparingPartitions/index.php?link=Tut9

jackknife <- function(x, theta.fun, jack.gf=as.character(1:length(x)), alpha=c(0.025,0.975), ...) {
	call <- match.call()
	thetahat <- theta.fun(x, ...)
	n <- length(unique(jack.gf))
	jackdist0 <- c()
	for (gfl in unique(jack.gf)) {
		xl <- x[-which(jack.gf==gfl)]
		jackdist0[gfl] <- theta.fun(xl, ...)
		}
	jackdist <- n*thetahat-(n-1)*jackdist0 # corrected for estimation bias
	
	est <- mean(jackdist)
	ci <- est+qt(alpha,n-1)*sqrt(var(jackdist)/n)
	return(list(thetahat=thetahat, est=est, confpoints=ci, jackdist=jackdist, jackdist.or=jackdist0, call=call))}
	
	
# To jackknife bi- and multivariate datasets,
# write theta so that its argument x is the set of observation indices
# and simply pass the vector 1,2,..n as data to jack.
# the multivariate dataset is then passed on as the additional argument xdata
# For example, to bootstrap
# the correlation coefficient from a set of 15 data pairs:
#       xdata <- matrix(rnorm(30),ncol=2)
#       n <- 15
#       theta.fun <- function(x,xdata){cor(xdata[x,1],xdata[x,2])} # x is a vector specifying row indices
#       results <- jack(1:n, theta.fun, xdata, gf)
