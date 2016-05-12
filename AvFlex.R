# Calculates average flexibility 
# Following Rolian 2009 and Marroig et al. 2009
# M is a Variance-covariance matrix
# B is a set of random vectors (skewers)
# B can be either provided or generated within (defult)

AvFlex <- function(M, n.it=1000, B=NULL) {	
	if (is.null(B)) {
		B <- matrix(rnorm (ncol(M)*n.it, mean = 0, sd = 1), ncol(M), n.it)}
	B <- t(t(B)/sqrt(colSums(B^2)))
	Z <- M%*%B
	Z <- t(t(Z)/sqrt(colSums(Z^2)))
	r <- diag(t(Z)%*%B)
	z <- mean(0.5*log((1+r)/(1-r)))
	(exp(z/0.5)-1)/(exp(z/0.5)+1)}

#AvFlex <- function(M, n.it=1000, B) {
#	B <- t(t(B)/sqrt(colSums(B^2)))
#	mean(diag(t(B) %*% M %*% B)/sqrt(diag(t(B) %*% (M %*% M) %*% B)))
#} # using the equation; doesn't include fisher's z transformation but it doesn't matter much