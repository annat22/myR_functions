# Calculates the statistic (W) for Levene's test for multivariate data, 
# using Brown-Forstythe modification (median instead of mean)
# X1 and X2 are the datasets whose variances are to be commpared

levinesBFW <- function(X1, X2) {
	Xm1 <- matrix(apply(X1,2,median),byrow=TRUE,nr=nrow(X1),nc=ncol(X1))
	y1 <- sqrt(diag((X1-Xm1)%*%t(X1-Xm1)))
	m1 <- mean(y1)
	Xm2 <- matrix(apply(X2,2,median),byrow=TRUE,nr=nrow(X2),nc=ncol(X2))
	y2 <- sqrt(diag((X2-Xm2)%*%t(X2-Xm2)))
	m2 <- mean(y2)
	wss <- sum((y1-m1)^2)+sum((y2-m2)^2)
	ass <- (m1-mean(c(y1,y2)))^2+(m2-mean(c(y1,y2)))^2
	W <- (ass/wss)*(Nt1+Nt2-2)
	W}
