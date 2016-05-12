# calculates disparity for a specific node in a tree
# Disparity is calculated as either the trace of the matrix ("trace")
# average square distance between all pairwise squared distances ("ave.sq")

tip.disp <- function (node, data, tree, disp=c("trace","avg.sq")) {
    Ntip <- length(tree$tip.label)
    nb.node <- which(tree$node.label==node)
    Xt <- data[tips(tree, Ntip+nb.node), ]
    if (disp=="trace") {d <- sum(diag(var(Xt)))}
    if (disp=="avg.sq") {d <- mean(dist(Xt, method = "euclidean")^2)}   
    d} 
