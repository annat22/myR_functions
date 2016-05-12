# scales a contingency table (Xob) by its chi-square deviations

chisq.dev <- function(Xob, exp=FALSE) {
	Mr <- diag(rowSums(Xob)) # row margins
	Mc <- diag(colSums(Xob)) # col margins
	T <- 1/matrix(sum(Xob), nrow(Xob), ncol(Xob))
	Xexp <- Mr%*%T%*%Mc
	st.dev <- (Xob-Xexp)/sqrt(Xexp)
	dimnames(st.dev) <- dimnames(Xexp) <- dimnames(Xob)
	if(exp==FALSE) res <- st.dev else res <- list("st.dev"=st.dev, "expected"=Xexp)
	res}
	
