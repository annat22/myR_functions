# Common subspace, or Krzanowski's method for comparing two covariance matrices
# following Zelditch et al., 2006, Evolution & Development, 8, 46-60 and green book p.308
# Returns a distance metric, whereas Krzanowski's original metric 
# is a similarity one (Blows et al., 2004, American Naturalist, 163, 329-340)
# This is the fastest and most stable calculation
# V1, V2 are the covaraince matrices, q is number of dimensions to retain

comsubsp <- function(V1,V2, q=floor(ncol(V1)/2), ...) {
 	E <- eigen(V1, symmetric=TRUE)$vectors
 	Q <- E[,1:q]%*%solve(E)[1:q,]
 	E <- eigen(V2, symmetric=TRUE)$vectors
 	R <- E[,1:q]%*%solve(E)[1:q,]
	j <- try(eigen(Q-R)$values, silent=TRUE)
	if (is.character(j)) NA else sqrt(sum(asin(as.numeric(j[1:q]))^2))
 	} # distance; Zelditch version 
		
