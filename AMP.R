## Aligns the specimen along the midline plane ($aligned).
## This function now works with both 3D and 2D; no need to specify which it is
## X is a matrix of k landmarks by m dimensions of one specimen;
# "midline", "right", and "left" are vectors indicating the names (appear 
# as row names of X; could also work with indices but not recommended) 
# of the midline, right and left landmarks, respectively.
# Missing data should be designated with NA
# "fillNA" specify the option to fill in missing lateral landmarks based on the other side
# Make sure the left and right landmarks are specified in the same order
# so that the first element in the "left" vector specifies the left equivalent 
# of the first element in the "right" vector, etc.
# If fillNA=TRUE (the default) then it fills in missing data by reflecting 
# from the other side. In this case $aligned will include the reflected landmarks.
# If both sides are missing for a certain landmark, that landmark is simply ignored and returned as NA.
# Missing midline landmarks are simply ignored and returned as NA
# It also provides the average of the left and right ($averaged).
# The resulting aligned specimen will have the midline axis as the last 
# column (z in 3D, y in 2D), the longest non-midline axis as the first column (x in 3D), 
# and the second non-midline axis as y for 3D.
# No further rotation and reflection is done, therefore the specimen may be facing any 
# direction along each of the axes.

# NOTE: the resulting configuration will be re-organized to have all the midline 
# landmarks first, then the right lateral ones, and then the left lateral ones. 
# It is therefore HIGHLY RECOMMENDED that landmarks will be designated by names in X
# (i.e., rownames(X) is a charcter vector rather than NULL)

AMP <- function(X, midline, right, left, fillNA=TRUE) { 
		ncl <- ncol(X) # determining whether it's 2D or 3D
		Xm <- na.omit(X[midline,])
		Xc <- t(t(X)-colMeans(Xm)) # translating all the LM's based on the centroid of the midline
		W <- na.omit(Xc[midline,]) 
		RM <-svd(var(W))$v # finding the rotation matrix
		Xa <- Xc%*%RM # rotating the whole specimen to align along the midline
		Mr <- Xa[right,]
		Ml <- Xa[left,]
		i <- which(is.na(Mr))
		Mlr <- Ml; Mlr[,ncl] <- -Mlr[,ncl]
		Mr[i] <- Mlr[i] # filling NA's in the right side by reflecting the left
		i <- which(is.na(Ml))
		Mrr <- Mr; Mrr[,ncl] <- -Mrr[,ncl]
		Ml[i] <- Mrr[i] # filling NA's in the left side by reflecting the right
		if (fillNA==TRUE) {	Xal <- rbind(Xa[midline,], Mr, Ml)} else {Xal <- Xa[c(midline,right,left),]}
		Mr[,ncl] <- -(Mr[,ncl]) # reflecting the right side onto the left for averaging
		Mav <- (Mr+Ml)/2 # averaged
		list(aligned=Xal, averaged=rbind(Xa[midline,], Mav))}
		
		