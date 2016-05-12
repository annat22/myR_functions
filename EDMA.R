

# X is either an array of landmark configurations (k landmarks by m dimensions by N specimens) 
# or a matrix of interlandmark distances (N specimens by k(k-1)/2 ILMD's if k is number of landmarks) 
# argument 'type' specifies if X is landmark configurations (="LM") or interlandmark distances (="ILMD")
# m is number of dimensions
# output is the mean form "omega", the matrix of interlandmark distances "ED" 
# and the EDMA variance-covariance matrix "V"
# requires functions "ILMD" and "convEDMA" (below)
EDMA <- function(X, type=c("LM","ILMD"), m=3) {
	if (type=="LM") {ED <- t(apply(X,3,ILMD)); k=nrow(X)}
	if (type=="ILMD") {ED <- X; k=(1+sqrt(1+8*ncol(X)))/2}
	N <- nrow(ED)
	Em <- apply(ED^2, 2, mean)
	S <- apply((t(ED^2)-Em)^2, 1, sum)/N
	if (m==2){om<-(Em^2-S)^0.25}
 	if (m==3){om<-(Em^2-1.5*S)^0.25}
 	A <- array(apply(ED, 1, convEDMA), dim=c(k,m,N))
	Bs<-array(NA, dim=c(k,k,N))
 	for (i in 1:N){
      Cc<-apply(A[,,i],2,mean)
      Ac<-t(t(A[,,i])-Cc)
      Bs[,,i]<-Ac%*%t(Ac)}
  	B<-apply(Bs, 1:2,mean)
  	M<-convEDMA(om)
  	V<-(B-M%*%t(M))/k
	list(omega=om, ED=ED, V=V)}
	
# X is a matrix of k landmarks by m dimesnions (one specimen configuration)
# output is a vector of all pairwise interlandmark distances for that specimen
ILMD <- function(X) {
		k <- nrow(X)
		if(is.null(rownames(X))){lnames <- paste("l",1:k,sep="")} else {(lnames <- rownames(X))}
		D <- as.matrix(dist(X))
		l1 <- lnames[combn(k,2)[1,]]
		l2 <- lnames[combn(k,2)[2,]]		
		Dnm <- paste(l1, l2, sep="-")
		d <- D[lower.tri(D)]; names(d) <- Dnm
		d}


convEDMA <- function(ed, m=3) {
	l <- length(ed)
	k <- (1+sqrt(1+8*l))/2
	mat<-diag(0,k)
	mat[row(mat)>col(mat)]<-ed; mat<-t(mat)
	mat[row(mat)>col(mat)]<-ed
	C1<-diag(k)-1/k*matrix(1,k,k)
  	B<- -0.5*C1%*%mat^2%*%C1
 	eC<-eigen(B)
  	eve<-eC$vectors
  	eva<-eC$values
  	MD<-matrix(NA, k, m)
  	for (i in 1:m) {MD[,i]<-sqrt(eva[i])*eve[,i]}
 	MD}


