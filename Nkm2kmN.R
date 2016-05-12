# Converts a matrix of N x k*m to an array of k x m x N. 
# Input (X) is a matrix of N specimens
# each specimen is a row vector o flandmarks arranged as {x1, y1, z1, x2, y2, z2...}

Nkm2kmN <- function(X, m=3) {
		k <- ncol(X)/m
		N <- nrow(X)
		AA <- array(apply(array(t(X), dim=c(m,k,N)), 3, t), dim=c(k,m,N))
		cn <- colnames(X)[seq(1,ncol(X),m)]
		if (m==2) mm <- c("x","y") else mm <- c("x","y","z")
		#dimnames(AA) <- list(substr(cn,1,nchar(cn)-2), mm, rownames(X))
		dimnames(AA) <- list(paste("LM", 1:k, sep=""), mm, rownames(X))
		AA}
