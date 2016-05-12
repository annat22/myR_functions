# Computes interlandmark distances
# X is a specimen matrix of k landmarks by m dimensions; works for both 2D and 3D
# LM is a matrix of kd x 2 where kd is the number of desired distances; each row in LM specifies the 
# row indices (or landmark names if applicable) in X for the pair of landmarks between which the
# distance is to be calculated. by default computes all pairwise distances
# last updated July 20 2010

ILMD <- function(X, LM=NULL) {
		k <- nrow(X)
		if(is.null(rownames(X))){lnames <- paste("l",1:k,sep=".")} else {(lnames <- rownames(X))}
		D <- as.matrix(dist(X))
		l1 <- lnames[combn(k,2)[1,]]
		l2 <- lnames[combn(k,2)[2,]]		
		Dnm <- paste(l1, l2, sep="-")
		d <- D[lower.tri(D)]; names(d) <- Dnm
		if (is.null(LM)) {d} else {
			a <- b <- cbind(match(LM[,1],lnames),match(LM[,2],lnames))
			o <- t(apply(a,1,order))
			a[which(o[,1]==2),1] <- b[which(o[,1]==2),2]
			a[which(o[,1]==2),2] <- b[which(o[,1]==2),1]
			d[paste(lnames[a[,1]],lnames[a[,2]], sep="-")]
			}}
