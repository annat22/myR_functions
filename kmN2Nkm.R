# Converts an array of k x m x N to a matrix of N x k*m.
# Input (XX) is an array of k x m x N 
# (e.g., N=number of specimens, k=number of landmarks, m=number of dimensions). 
# Output is a matrix of N specimens, where each specimen is a row vector arranged as {x1, y1, z1, x2, y2, z2...}

kmN2Nkm <- function(XX) {
	m <- ncol(XX)
	X <- t(apply(XX,3,t))
	if (!is.null(dimnames(XX)[[1]])) {dimnames(X) <- list(dimnames(XX)[[3]], paste(rep(rownames(XX), each=m), rep(colnames(XX),m), sep="_"))}
	X}
