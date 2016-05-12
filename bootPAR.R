# parametric bootstrap
# S1 and S2 are covariance matrices
# S2 is optional
# func is the function to apply on S1 (and S2)
# n is a vector of 2 or 1 elements indicating sample size to be drawn. 
# if n is 2 elements then two samples are drawn and the function is assumed to require two samples.
# The two samples can be either from the same matrix, if S2 isn't provided, 
# or one from S1 and one from S2, if S2 is provided.
   
 samppa <- function(N, S1, func, S2=NULL, ...) {
		X1 <- rmvnorm(N[1], sigma=S1)
		if (length(N)==2) {
			if (is.null(S2)) {
				X2 <-rmvnorm(N[2], sigma=S1)
				} else {
					X2 <- rmvnorm(N[2], sigma=S2)}
					func(var(X1), var(X2), ...)}
		else {func(var(X1), ...)}
		}

bootPAR <- function(S1, S2=NULL, n, tail=2, n.boot=500, alpha=0.05, func, ...) {
		N <- matrix(n, nr=length(n), nc=n.boot)
		dist <- sort(apply(N, 2, samppa, S1=S1, S2=S2, func=func, ...))
		dist}
