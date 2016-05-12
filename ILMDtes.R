# generates a definition matrix for interlandmark distances 
# based on a tesselation representation of a 2D or 3D configuration.
# X is k x m matrix of landmark coordinates of one configuration
# of k landmarks and m dimensions. (presumably mean configuration of a sample)
# if X is symmetric then only one side of the configuration should be included 
# (e.g., only midline and right bilateral landmarks); otherwise symmetry is ignored.
# output is a matrix with two columns. each row is a pair of landmarks that define one interlandmark distance
# use ILMD.R function in order to calculate the interlandmark distances themselves from the definition matrix generated here.

require(geometry)
ILMDtes <- function(X) {
	if (is.null(rownames(X))) rownames(X) <- paste("l",1:nrow(X),sep=".")
	tes <- delaunayn(X)
	if (ncol(X)==3) {ild0 <- rbind(tes[,1:2], tes[,2:3], tes[,3:4], tes[,c(1,4)])}
	if (ncol(X)==2) {ild0 <- rbind(tes[,1:2], tes[,2:3], tes[,c(1,3)])}
	ir <- which(ild0[,1]>ild0[,2])
	ild0[ir,] <- ild0[ir,2:1]
	ild.nm <- rownames(ild0) <- apply(ild0,1,paste, collapse="-")
	uil <- unique(ild.nm)
	ild <- ild0[uil,]
	ild <- ild[order(ild[,1]),]
	cbind(rownames(X)[ild[,1]], rownames(X)[ild[,2]])}