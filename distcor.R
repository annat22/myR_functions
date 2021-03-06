# Generates either a distogram or an autocorrelogram based on 
# Hardy, O. J., and S. Pavoine. 2012. Evolution 66:2614-2621 
# and Diniz-Filho, J. A. F. 2001. Evolution 55:1104-1109
# code is modified from Hardy and Pavoine 2012
# includes a comparison to Brownian motion model (simulations scaled by the original data)
# and to a case of no phylogenetic structure (permutations of the original data)
# X is a data matrix (for multivariate data) or vector (for univariate data) and tree as a phylo object
# dcfun specifies whether the Euclidean distance ("ED") is calculated or Moran's I ("IM")
# nsims specifies the number of simulations for the BM simulations and permutations
# default (plot=TRUE) generates the plot itself (either distogram or autocorrelogram), 
# otherwise you can set plot=FALSE and plot it yourself based on the output
# output is a matrix in which every row is a class of patristic distance over which 
# the distance/similarity is averaged. columns are the mean patristic distance for that class ("pdcm"), 
# the observed value of distance/similarity ("obs"), the mean ("Pm"), lower ("Pl"), and upper ("Pu") 
# bounds of the 95% CI of permutations, and the mean ("BMm"), lower ("BMl"), and upper ("BMu") 
# bounds of the 95% CI of BM simulations

require(abind); require(geiger)

distcor <- function(X, tree, dcfun=c("ED","IM"), nsims = 100, plot=TRUE) {
	X <- as.matrix(X)
	pd <- cophenetic.phylo(tree)
	pd <- pd/max(pd)
	class <-c(0, (1:10)/10)
	if (dcfun=="ED") class <- class^2
	pdcm <- c()
	for(j in 2:length(class)) {pdcm <- c(pdcm, mean(pd[pd>class[j-1] & pd<=class[j]]))}
	XP <- array(apply(array(X, dim=c(nrow(X),ncol(X),nsims)), 3, function(X){X[sample(1:nrow(X),nrow(X)),]}), dim=c(nrow(X),ncol(X),nsims)) #permutations of species to simulate no phylogenetic structure
	XS <- sim.char(tree, var(X), nsims=nsims) # simulations of BM
	XX <- abind(XP, XS, X, along=3)
	if (dcfun=="ED") dfun <- function(X){as.matrix(dist(X))} else if (dcfun=="IM") dfun <- function(X){X%*%t(X)/sum(diag(var(X))) + 1/(nrow(X)-1)}
	cM0 <- c()
	for (h in 1:dim(XX)[3]) {
		X <- as.matrix(XX[,,h])
		D <- dfun(X)
		cm <- c()
		for(j in 2:length(class)) {cm <- c(cm, mean(D[pd>class[j-1] & pd<=class[j]]))}
		cM0 <- cbind(cM0, cm)			
		}
	na <- which(is.na(pdcm))
	cMP <- cM0[-na,1:nsims]
	cMS <- cM0[-na,(nsims+1):(nsims*2)]
	cM <- cbind("pdcm"=pdcm[-na], "obs"=cM0[-na,ncol(cM0)], "Pm"=rowMeans(cMP), "Pl"=apply(cMP, 1, quantile, 0.025), "Pu"=apply(cMP, 1, quantile, 0.975), "BMm"=rowMeans(cMS), "BMl"=apply(cMS, 1, quantile, 0.025), "BMu"=apply(cMS, 1, quantile, 0.975))
	if (plot==TRUE) {
		ltys <- c(NA,1,3,3,3,1,1,1)
		cols <- c(NA,"Black","black","grey","grey","black","grey","grey")
		if(dcfun=="ED") {
			cM[,"pdcm"] <- sqrt(cM[,"pdcm"])
			legpos <- "bottomright"
			yttl <- "Euclidean distance"
			xttl <- "Square root of patristic distance"
			} 
		else if(dcfun=="IM") {
			legpos <- "topright"
			yttl <- "Moran's I (trait similarity)"
			xttl <- "patristic distance"
			}
		plot(cM[,"pdcm"], cM[,"obs"], ylim=c(min(0,min(cM[,-1])),max(cM[,-1])), pch=19, xlab=xttl, ylab=yttl)
		for (j in 2:ncol(cM)) {lines(cM[,"pdcm"], cM[,j], lty=ltys[j], col=cols[j])}
		legend(legpos, legend=c("Observed","Permutated","Brownian"), cex=0.8, lty=c(1,3,1), pch=c(19,NA,NA))
		}
	cM
	}

