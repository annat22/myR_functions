#calculates simmilarity between two covariance matrices (V1 and V2) 
# using the random skewer method,following Cheverud 1996 J.Evol.Biol 9:5-42 
# as explained in Cheverud and Marroig 2007 Genet.Mol.Biol. 30(2):461-469 
# and Marroig et al 2009 Evol.Biol 36:136-148; 
# skewers are drawn from a normal distribution following Marroig et al. 2012 Evo.Bio. 38:225-241

randsk <- function(V1, V2, n.it=1000, B=NULL) {	
	if (is.null(B)) {B <- matrix(rnorm (ncol(V1)*n.it, mean = 0, sd = 1), ncol(V1), n.it)} # generating a sample of selection vectors
	B <- t(t(B)/sqrt(colSums(B^2))) # scaling them to unit length
	Z1 <- V1%*%B # response vectors of first VCV matrix
	Z2 <- V2%*%B # response vectors of second VCV matrix
	Z1 <- t(t(Z1)/sqrt(colSums(Z1^2))) # scaling response vectors of V1
	Z2 <- t(t(Z2)/sqrt(colSums(Z2^2))) # scaling response vectors of V2
	r <- diag(t(Z1)%*%Z2) # caluculating their dot-product
	z <- mean(0.5*log((1+r)/(1-r))) # taking their mean using fisher's transformation
	round((exp(z/0.5)-1)/(exp(z/0.5)+1), 2) # re-transforming the mean to vary between -1 and 1
	}


