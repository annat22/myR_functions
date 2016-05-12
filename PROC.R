# Partial procrustes superimposition. 
# "PROC" does the actuall superimpostion. The other two functions are auxiliary:
# "scale.centroid" scales a centered configuration by its centroid size 
# "fitX" performs one round of fitting one configuration to another. 
# Missing data is ignored but should designated by NA in the data matrix.
# AA is an array of the k x m x N format (k landmakrs, m dimensions, N specimens).  
# Return the Procrustes coordinates ($rotated) and the optimization statistic ($opt.stat).

scale.centroid <- function(X) {X/sqrt(sum(X^2, na.rm=TRUE))} # X is centered already

fitX <- function(X, Ref) {	
		X[is.na(Ref)] <- NA
		Ref[is.na(X)] <- NA
		M <- t(na.omit(Ref)) %*% na.omit(X)
		SVD <- svd(M)
		V <- SVD$v # eigenvectors scaled to 1 unit length
		U <- SVD$u
		L <- diag(SVD$d)
		S <- ifelse(L<0, -1, L)
		S <- ifelse(L>0, 1, L)
		RM <- V %*% S %*% t(U) # the rotation matrix
		X %*% RM} # rotate X onto the reference


PROC <- function(AA, sizescaled=TRUE) {
	Xc <- array(apply(AA, 3, scale, scale=FALSE), dim=dim(AA)) # centered
	if (sizescaled==TRUE) {Xs <- array(apply(Xc, 3, scale.centroid), dim=dim(AA))} else {Xs <- Xc}
	R <- Xs[,,1]
	d <- 1
	while (d>0.001) {
		Xf <- array(apply(Xs, 3, fitX, Ref=R), dim=dim(AA), dimnames=dimnames(AA))
		MC <- colMeans(aperm(Xf, c(3,1,2)), na.rm=TRUE) # mean configuration
		MC[is.na(R)] <- NA
		d <- sqrt(sum((MC-R)^2, na.rm=T))
		R <- MC
		Xs <- Xf}
	list(rotated=Xf, opt.stat=d)}
