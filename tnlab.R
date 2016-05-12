# Attaches taxon labels to the edges matrix. 
# "phy" is the phylogenetic tree as an ape object
# returns the edges matrix with an additional column containing label for the taxon listed in the second column
# as given by tree$tip.label and tree$node.label 

tnlab <- function(phy) {
	ntip <- length(phy$tip.label)
	nnd <- phy$Nnode
	if (is.null(phy$node.label)) {phy$node.label <- as.character((ntip+1):(ntip+nnd))}
	tnlab <- character(nrow(phy$edge))
	tnlab[which(phy$edge[,2]<=ntip)] <- phy$tip.label
	tnlab[which(phy$edge[,2]>ntip)] <- phy$node.label[-1]
	cbind(phy$edge, tnlab)}
