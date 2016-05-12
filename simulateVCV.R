# simulates a correlation or a variance-covariance matrix with a given mean correlation
# p is number of parameters, q is mean correlation in the matrix
# cv is the coefficient of variance of the off-diagonal elements
# matrix specifies whether output is correlation ("R") or v/cv ("V") matrix
# if matrix="V" variances are drawn from a chi-square distribution 
simulateVCV <- function(q, p, cv=min((1-q)^2/q, 0.9), matrix=c("R","V")) {
	R <- matrix(1,p,p)
	qz <- 0.5*log((1+q)/(1-q))
	cv <- ifelse((1-q)^2/q<cv, (1-q)^2/q, cv)
	z <- rnorm(p*(p-1)/2, mean=qz, sd=cv*qz)	
	R[lower.tri(R)] <- (exp(z/0.5)-1)/(exp(z/0.5)+1)
	rt <- t(R)
	R[upper.tri(R)] <- rt[upper.tri(rt)]
	if (matrix=="V") {
		sd <- sqrt(rchisq(p,p))
		return(sweep(sweep(R, 1, sd, "*"), 2, sd, "*"))
		} else return(R)	
	}