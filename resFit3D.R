# AA is an array of k landmarks x 3 dimensions x N specimens assumed to be already superimposed by a least square fit procedure.
# missing landmarks will be ignored throughout the sample, so all the configurations will be superimposed using only landmarks that occur in all of the configurations.
# The function RF3D carries a 3D generalized resistant fit superimposition based on Rholf and Slice 1990 Syst.Zool 39(1):40-59 and Slice 1996 Advances in Morphometrics p.179-199

tripletS <- function(triplet, ref, tar) {
	R <-t(ref[triplet,])
	T <- t(tar[triplet,])
	u1 <- (R[,2]-R[,1])/sqrt((R[,2]-R[,1])%*%(R[,2]-R[,1]))
	v1 <- (T[,2]-T[,1])/sqrt((R[,2]-R[,1])%*%(R[,2]-R[,1]))
	u2a <- (R[,3]-R[,1])-((R[,3]-R[,1])%*%u1)*u1
	u2 <- u2a/sqrt(u2a%*%u2a)
	v2a <- (T[,3]-T[,1])-((T[,3]-T[,1])%*%v1)*v1
	v2 <- v2a/sqrt(v2a%*%v2a)
	u3 <- c(u1[2]*u2[3]-u1[3]*u2[2], u1[3]*u2[1]-u1[1]*u2[3], u1[1]*u2[2]-u1[2]*u2[1])
	v3 <- c(v1[2]*v2[3]-v1[3]*v2[2], v1[3]*v2[1]-v1[1]*v2[3], v1[1]*v2[2]-v1[2]*v2[1])
	H <- t(cbind(u1,u2,u3))%*%cbind(v1,v2,v3)
	solve(H+diag(c(1,1,1)))%*%(H-diag(c(1,1,1)))} # one triplet

fitXR <- function(X, ref) {	
		tar<- X/median(apply(as.matrix(dist(X)), 2, median))
		k <- nrow(ref)
		tau <- median(apply(as.matrix(dist(ref))/as.matrix(dist(tar)), 2, median, na.rm=TRUE))
		tri <- combn(1:k, 3)
		SS <- matrix(apply(tri, 2, tripletS, ref=ref, tar=tar), 9, ncol(tri))
		Shat <- matrix(apply(SS, 1, median), 3,3)
		RM <- (diag(c(1,1,1))+Shat)%*%solve(diag(c(1,1,1))-Shat)
		tau*tar%*%RM+matrix(apply(ref-tau*tar%*%RM, 2, median), byrow=TRUE, nr=k, nc=3)} # fits one configuration onto another

# the generalized procedure:
ersFit3D <- function(XX) {
	i <- unique(u <- which(is.na(XX), arr.ind=TRUE)[,1]) # row number of missing landmarks
	if(length(i)>0) AA<-XX[-i,,] else AA<-XX 
	M <- apply(AA, 1:2, median)
	M <-M/median(apply(as.matrix(dist(M)), 2, median))
	d <- 1
	dd <- 1
	itr <- 0
	while (d>0.005 & itr<3) {
		Xf <- array(apply(AA, 3, fitXR, ref=M), dim=dim(AA), dimnames=dimnames(AA))
		MC <- apply(Xf, 1:2, median)
		d <- median(sqrt(diag((MC-M)%*%t(MC-M))))
		M <- MC
		AA <- Xf
		itr<-itr+1
		dd <- append(dd, d)}
	if(length(i)>0) {
		YY <- array(NA, dim=dim(XX), dimnames=dimnames(XX))
		for (h in 1:dim(YY)[3]) { x <- c(t(AA[,,h]))
				for (j in 1:length(i)) { 
					x <- append(x,rep(NA,ncol(AA)),after=ncol(AA)*(i[j]-1))}
				YY[,,h] <- matrix(x, byrow=TRUE, nr=nrow(XX), nc=ncol(XX))}
		AA <- YY}
	list(rotated=AA, opt.stat=d, itr=itr)}