# caluclates the expected correlations between ILMD's due to their shared LM (absolute values); as well as the observed correlations (absolute values) and the vector signifies which corelation element has a shared landmark (0 if they don't share and 1 if they do)
# input is an array of k landmarks by m dimensions by N specimens with landmark names as rownames 

corILMD <- function(XX) {
	ILMD <- function(X, lm=NULL) {
		k <- nrow(X)
		D <- c()
		Dnm <- c()
		for (i in 1:(k-1)) { 
			for (j in (i+1):k) {
				D <- append(D, sqrt(sum((X[i,]-X[j,])^2)))
				Dnm <- append(Dnm, paste(rownames(X)[i],"-", rownames(X)[j], sep=""))
					}}
			names(D) <- Dnm} 
		D}
	M <- 	log(t(apply(XX, 3, ILMD)))
	R <- S <- cor(M)
	cn <- matrix(unlist(strsplit(colnames(M),"-")),2,ncol(R))
	for (i in 1:(ncol(cn)-1)) { 
		for (j in (i+1):ncol(cn)){
			S[i,j] <- ifelse(any(cn[,i]%in%cn[,j]),1,0)}}
	r.obs <- abs(R[upper.tri(R)])
	s <- S[upper.tri(S)]
	expcor.ild <- function(l1,l2,lshared) {
		c1 <- l1-lshared
		c2 <- l2-lshared
		abs(0.5*(c1%*%c2/(sqrt(c1%*%c1)*sqrt(c2%*%c2))))}
	expcorIL <- function(X, dn) {
		n <- ncol(dn)
		Re <- matrix(nr=n,nc=n)
		for (i in 1:(n-1)) { 
			for (j in (i+1):n){
				l <- as.factor(c(dn[,i],dn[,j]))
				si <- which(tabulate(l)==2)
				if (length(si)>0) {
					ls <- levels(l)[si]
					l1 <- levels(l)[-si][1]
					l2 <- levels(l)[-si][2]
					Re[i,j] <- expcor.ild(X[l1,],X[l2,],X[ls,])} else {Re[i,j] <- 0}
			}}
	Re[upper.tri(Re)]}
	rem <- rowMeans(apply(XX,3,expcorIL, dn=cn))
	list(allILMD=M, r.obs=r.obs, r.exp=rem, shared=s)}