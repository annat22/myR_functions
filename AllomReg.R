# Allometric regression following Rolian 2009 and Lleonart et al. 2000
# using logged geometric mean of the raw variables (arythmetric mean of the logged variables)
# Input (E) is the raw (not log-transformed) interlandmark distances. 
# Output is the traits theoretical values at the group average size (mean GM).
# If scale=TRUE (default) then the data is standardized by the group average size so that all individuals are scaled to a theoretical body size of 1.

AllomReg <- function(E, scale=TRUE) {
	gm <- apply(E, 1, function(x){prod(x)^(1/length(x))})
	GMi <- matrix(gm, nr=nrow(E), nc=ncol(E))
	GMo <- matrix(mean(gm), nr=nrow(E), nc=ncol(E))
	B <- matrix(coef(lm(log(E)~log(gm)))[2,],  nr=nrow(E), nc=ncol(E), byrow=TRUE)
	Er <- (E*(GMo/GMi)^B)
	if (scale==TRUE) {res <- Er/GMo} else {res <- Er}
	res}
