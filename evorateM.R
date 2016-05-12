# Calculates the evolutionary rate matrix following Revell (2009) and Revell and Collar (2010)
# Input is the data matrix (X) and the tree in ape's phylo format
# rownames(X) should include the tip labels of tree but order does not need to be the same
# both tree and X can contain taxa that do not appear in the other object, these taxa will be ignored
# Output is a symmetric matrix with nrow and ncol equal to ncol(X)

require(ape)
evorateM <- function(X, phy) {
			tx <-  phy$tip.label
			Nt <- length(tx)
			Xi <- X[tx,]
			Se <- vcv.phylo(phy)[tx,tx] # C in Revell (2009)
			invS <- solve(Se)
			one <- matrix(1,Nt,1)
			a <- t(t(one)%*%invS%*%Xi)*sum(invS)^-1
			Re <- t(Xi-one%*%t(a))%*%invS%*%(Xi-one%*%t(a))/(Nt-1)
			Re}