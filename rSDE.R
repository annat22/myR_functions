# calculates the relative standard deviation of eigenvalues based on Van Valen (1974)
# for either covariances or correlation matrices (M)
# Covariance matrices are standardized by the total variance
# For correlation matrices this reduces to rSDE of Pavlicev et al. (2009) 

rSDE <- function(M) {
	d <- eigen(M, symmetric=TRUE)$values
	p <- length(d)
	sqrt(sum((d-mean(d))^2)*p/(sum(d)*sum(d)*(p-1)))}
