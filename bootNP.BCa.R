# Calculates the observed theta, bootstraps the sample, and calculates BCa CI's
# x is a vector of either univariate data or 1:nrow(X) for a multivariate dataset X 
# (see example below), in which case X would be passed on as an additional argument of theta.fun
# For theta that includes comparison of two multivariate datasets X would 
# be rbind(X1,X2) and theta.fun would have an additional argument specifying grouping factor. 
# If replace=FALSE then data is permuted instead of bootstrapped

# Output is a list including: 
# the observed estimate (thetahat), 
# a vector of pseudovalues (thetastar), 
# the confidence interval for the specified alpha, including BCa if BCa=TRUE,
# the BCa parameters, z0, acc, and u, if BCa=TRUE,

### modified from http://www-rohan.sdsu.edu/~babailey/stat672/bootstrap.r
### and http://www.rohan.sdsu.edu/~babailey/stat672/bcanon.r

# To bootstrap bi- and multivariate datasets,
# write theta.fun so that its argument x
# is a set of observation indices
# and simply pass as data to bootstrap the vector 1,2,..n.
# For example, to bootstrap
# the correlation coefficient from a set of 15 data pairs:
#       xdata <- matrix(rnorm(30),ncol=2)
#       n <- 15
#       theta.funCor <- function(x,xdata){cor(xdata[x,1],xdata[x,2])} 
#       results <- bootNP.BCa(x=1:n, nboot=20, theta.fun=theta.funCor, xdata)


bootNP.BCa <- function(x, nboot, alpha=c(0.025,0.975), replace=TRUE, BCa=TRUE, theta.fun, ...){
	call <- match.call()
	thetahat <- theta.fun(x,...) # observed value
	n <- length(x)
	
	# Checking if a grouping factor is supplied for theta.fun (for comparing two multivariate datasets) 
	# and resample accordingly 
	args <- list(...)
	if (!"gf"%in%names(args)) { # no grouping factor, so all one group
		bootsam <- apply(matrix(x, nrow=n, nc=nboot), MARGIN=2, sample, size=n, replace=replace)
		} else {
			gf <- args["gf"] # grouping factor supplied for theta.fun
			ngf <- table(gf) # number of groups in the grouping factor
			bootsam <- c()
			i=1
			for (j in 1:length(ngf)){
				bj <- apply(matrix(i:(i+ngf[j]-1), nrow=ngf[j], nc=nboot), MARGIN=2, sample, size=ngf[j], replace=replace)
				i <- i+ngf[j]
				bootsam <- rbind(bootsam, bj)
				} # ensuring that original sample size is maintained within each group
			}
		
	thetastar <- apply(bootsam, MARGIN=2, theta.fun,...) # pseudovalues
 	confpoints <- NULL; z0 <- NULL; acc <- NULL; u=NULL
 	if (BCa==TRUE) {
 		z0 <- qnorm(sum(thetastar<thetahat, na.rm=TRUE)/(nboot+1))
	   	u <- rep(0,n)
   		for(i in 1:n) {u[i] <- theta.fun(x[-i],...)}  # Jackknife pseudovalues
   		bias <- mean(u, na.rm=TRUE)-u
   		acc <- sum(bias^3)/(6*(sum(bias^2))^1.5) # acceleration
   		zalpha <- qnorm(alpha)
   		cf <- pnorm(z0+(z0+zalpha)/(1-acc*(z0+zalpha))) # correction factors for both alphas
   		confpoints <- quantile(thetastar, probs=cf, na.rm=TRUE)
   		names(confpoints) <- paste("BCaCI", alpha, sep="")
 	} else confpoints <- quantile(thetastar, probs=alpha, na.rm=TRUE)
 	return(list(thetahat=thetahat, thetastar=thetastar, confpoints=confpoints, z0=z0, acc=acc, u=u, call=call, theta.args=names(list(...)), groupingFactor=ifelse(exists("gf"), gf, NA)))}



