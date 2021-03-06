#### AvCondE
# Calculated average conditional evolvability following Hansen and Haoule 2008 using simulations
# V is a VCV matrix; the set of random vectors (B) can be either determined in advance or generated within (defult)
require(MASS)
AvCondE <- function(V, n.it=1000, B=NULL, ...) {	
	if (is.null(B)) {
		B <- matrix(rnorm (ncol(V)*n.it, mean = 0, sd = 1), ncol(V), n.it)
		B <- t(t(B)/sqrt(colSums(B^2))) # standardized to unit length
		}
		mean(1/diag(t(B)%*%ginv(V)%*%B))
		}
	
