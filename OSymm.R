# This function uses object symmetry to calculate a symmetric configuration 
# and reconstruct missing lateral landmarks, following
# Klingenberg et al. 2002 (Evolution 56:1909-1920)

# This function now works with both 3D and 2D; no need to specify which it is
# X is a matrix of k landmarks by m dimensions of one specimen
# "midline", "right", and "left" are vectors indicating the names of the 
# midline, right and left landmarks, respectively. 
# (i.e., c(midline, right, left) should match rownames(X), not necessarily same order
# missing data should be designated with NA

# The output includes the original configuration ($rec.orig) and the symmetric consensus
# configuration ($symmconf). Both have missing landmarks reconstructed
# based on the other side if available. 
# Missing landmarks that cannot be reconstructed based on bilateral symmetry 
# (e.g.,midline landmarks or both sides are unavailable)
# are simply ignored and returned as NA.
# No further rotation or reflection is done, therefore don't expect the specimens to be 
# oriented in any particualr direction. To achieve that, please use the AMP function.

# NOTE: the resulting configuration will be re-organized to have all the midline 
# landmarks first, then the right lateral ones, and then the left lateral ones. 
# It is therefore HIGHLY RECOMMENDED that landmarks will be designated by names, so
# rownames(X) is a charcter vector rather than NULL, matching c(midline, right, left))

OSymm <- function(X, midline, right, left) {
	ncl <- ncol(X); nrw <- nrow(X)
	Xr <- cbind(X[,-ncl], -X[,ncl]) # the last axis is being reflected: z for 3D, y for 2D
	## reorganizing the original configuration and renaming the reflected one
	Xo.na <- Xo <- X[c(midline, right, left),]
	Xr.na <- Xr <- Xr[c(midline, left, right),]
	rownames(Xr.na) <- rownames(Xr) <- c(midline, right, left)
	## finding the rotation matrix between the original and the reflected:
	# landmarks that are missing for at least one of the sides are ignored 
	# for finding both the centoid and the rotation matrix
	Xo.na[is.na(Xr)] <- NA; Xr.na[is.na(Xo)] <- NA
	# centroids:
	mr <- mo <- matrix(colMeans(X,na.rm=TRUE), byrow=TRUE, nr=nrw, nc=ncl)
	mr[,ncl] <- -mr[,ncl]
	# Finding the rotation matrix:
	SVD <- svd(t(na.omit(Xr.na-mr)) %*% na.omit(Xo.na-mo))
	L <- diag(SVD$d)
	S <- ifelse(L<0, -1, L)
	S <- ifelse(L>0, 1, L)
	RM <- SVD$v %*% S %*% t(SVD$u) # the rotation matrix
	## centering and rotating the original and reflected configurations
	# including all LM's:
	Xo.c.rot <- (Xo-mo) %*% RM                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 
	Xr.c <- Xr-mr
	## symmetric configuration as the mean of the original and the reflected copies
	symm.conf <- apply(array(c(Xr.c, Xo.c.rot), dim=c(nrw,ncl,2)), 1:2, mean, na.rm=TRUE)
	dimnames(symm.conf) <- dimnames(Xo.c.rot)
	## reconstructing the missing lateral landmarks based on the other side
	Xo.c.rot[which(is.na(Xo.c.rot))] <- Xr.c[which(is.na(Xo.c.rot))]	
	list(rec.orig=Xo.c.rot, symm.conf=symm.conf)}
	
	