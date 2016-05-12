# shrinks the matrix to avoid singularity
# M is a covariance matrix; tol is shrinking tollerance
# returns the shrunk matrix
# for both correlation and covariance matrices
# Based on Jorjani, H., L. Klei, and U. Emanuelson. 2003. A Simple Method for Weighted Bending of Genetic (Co)variance Matrices. Journal of Dairy Science 86:677-679. 

shrink <- function(M, tol=10^-8) {
    ei <- eigen(M, symmetric=TRUE)
    d <- ei$values
    rtol <- tol * mean(d)
    if (min(d) < rtol) {
    	if (sum(diag(M))==ncol(M)) {
    		di <- d
    		di[d<rtol] <- 2*rtol
    		di <- di*(sum(d)/sum(di))
    		Ms <- ei$vectors %*% (di * t(ei$vectors))
    		Ms <- cov2cor(Ms)
    		} else {
    			d[d<rtol] <- rtol
    			Ms <- ei$vectors %*% (d * t(ei$vectors))
    			}
    	dimnames(Ms) <- dimnames(M)
       	return(Ms)} else {return(M)}
   }
