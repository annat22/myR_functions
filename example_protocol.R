## An example for using functions unifyVD, AMP, and OSymm
## This protocol and the functions included can be used with both 2D and 3D data
# input is in the k x m*N format; k landmakrs, m dimensions, N specimens; 
# the first row is specimens ID's at the top of each x-coordinate and is not read into X, the first column is LM's names
# in the input example file there are three specimens. Landmarks can be in any order regarding dorsal, ventral, midline etc. 
# as long as their names (not indices) are specified in the proper places below (see annotations)
# missing data should be designated with NA
# NOTE: all three functions reorganize the landmarks in their output which is why it's recommended to use the landmark names 
# rather than their indices
# Last updated Feb 28th 2011

source("unifyVD.R") # calling for the unifyVD function; it can be in any folder in which case you'd specify the whole path. e.g., "/Users/username/Desktop/my_functions/unifyVD.R"
source("AMP.R")
source("OSymm.R")

input <- "ex input.txt" # specify any name (or path) you want for the input file
m <- 3 # specify number of dimensions

cor <- c(); if (m==3) {cor <- c("x","y","z")}; if (m==2) {cor <- c("x","y")}
X <- read.table(input, row.names=1, skip=1) 
l1 <- scan(input, sep="\t", what="character", nlines=1) #reads first line of file which contains specimen ID No.'s
id <- l1[which(l1!="")] #yields the vector of specimens id's
k1 <- nrow(X) #number of landmarks including left and right and two of each of the dorsal-ventral common ones
p <- k1*m
N <- ncol(X)/m # number of specimens; also length(id)
XX <- array(as.matrix(X), dim=c(k1, m, N), dimnames=list(LM=rownames(X), cor=cor, spec_id=id)) #each matrix in the array is now a specimen

# sort by specimen number for easy linking with other specimen properties:
#i <- order(id, index.return = TRUE)
# id <- id[i] 
# XX <- XX[,,i]

#### specifying the landmarks; the landmark names in the 'left' vector is assumed to be given in the same order as the equivalent landmarks in the 'right' vector, and the same for the comV and comD vectors.
dorsal <- c("tPM", "IOF", "P3", "abM1", "adJS", "vZA", "PSO", "PS(as)", "JFO", "LFO", "JLO", "aLJ", "pMN", "SNa", "aN", "FFN", "FFP", "SOP", "tPM_L", "IOF_L", "P3_L", "abM1_L", "adJS_L", "vZA_L", "PSO_L", "PS(as)_L", "JFO_L", "LFO_L", "JLO_L", "aLJ_L", "pMN_L", "SNa_L", "aN_L")
ventral <- c("V1", "V3", "BMF", "OMF", "SOP_V","BOC", "F8", "PS(as)_V", "JFO_V", "JLO_V", "pJM", "pbM3", "P3_V", "P2dia", "BOC_L", "F8_L", "PS(as)_LV", "JFO_LV", "JLO_LV", "pJM_L", "pbM3_L", "P3_LV", "P2dia_L")
comD <- c("SOP", "PS(as)", "JFO", "JLO", "P3", "PS(as)_L", "JFO_L", "JLO_L", "P3_L") # the dorsal common LM's
comV <- c("SOP_V", "PS(as)_V", "JFO_V", "JLO_V", "P3_V", "PS(as)_LV", "JFO_LV", "JLO_LV", "P3_LV") # the ventral common LM's
midline <- c("FFN", "FFP", "SOP", "V1", "V3", "BMF", "OMF", "SOP_V")
right <- c("tPM", "IOF", "P3", "abM1", "adJS", "vZA", "PSO", "PS(as)", "JFO", "LFO", "JLO", "aLJ", "pMN", "SNa", "aN", "BOC", "F8", "PS(as)_V", "JFO_V", "JLO_V", "pJM", "pbM3", "P3_V", "P2dia")
left <- c("tPM_L", "IOF_L", "P3_L", "abM1_L", "adJS_L", "vZA_L", "PSO_L", "PS(as)_L", "JFO_L", "LFO_L", "JLO_L", "aLJ_L", "pMN_L", "SNa_L", "aN_L", "BOC_L", "F8_L", "PS(as)_LV", "JFO_LV", "JLO_LV", "pJM_L", "pbM3_L", "P3_LV", "P2dia_L")

##### Attaching ventral and dorsal
VD <- array(NA, dim=c(k1, m, N), dimnames=list(LM=c(dorsal, ventral), cor=cor, spec_id=id)) # array into which the unified data is entered
Er_vd <- matrix(nr=length(comD), nc=N, dimnames=list(comD, id)) # matrix of the unification errors
for (h in 1:N) {
		vd <- unifyVD(XX[,,h], dorsal, ventral, comD, comV)
		VD[,,h] <- vd$unified
		Er_vd[,h] <- vd$errors}


centsize <- function(X, ncl=m){
	M <- matrix(colMeans(X, na.rm=TRUE), byrow=TRUE, nr=nrow(X), nc=ncl)
	sqrt(sum((X-M)^2, na.rm=TRUE))}


CS <- matrix(apply(VD, 3, centsize), byrow=TRUE, nr=length(comD), nc=N, dimnames=list(comD, id))
sdErvd <- 100*Er_vd/CS # unification errors standardized by centroid size

write.table(sdErvd, file="sdErvd.tab", sep="\t")
write.table(VD, file="VD.tab", sep="\t")
# Unless you need to work with these files in other software, it is actually better to save objects as R objects

########## aligning along the midline plane and fill in NA's by reflecting from other side
AA <- array(NA, dim=c(k1, m, N), dimnames=list(LM=c(midline, right, left), cor=cor, spec_id=id))
for (h in 1:N) { 
	AA[,,h] <- AMP(VD[,,h], midline, right, left, fillNA=TRUE)$aligned}

write.table(AA, file="aligned_noNA.tab", sep="\t")

which(is.na(AA), arr.ind=TRUE)

######## alternatively you can use object symmetry (following Klingenberg et al. 2002 Evolution 56:1909-1920) to fill in missing landmarks and then AMP to align along the plane of symmetry

AA <- array(NA, dim=c(k1, m, N), dimnames=list(LM=c(midline, right, left), cor=cor, spec_id=id))
for (h in 1:N) { 
	AA[,,h] <- AMP(OSymm(VD[,,h],midline, right, left)$rec.orig, midline, right, left, fillNA=FALSE)$aligned}

plot(AA[,1,1],AA[,m,1], pch=NA)
text(AA[,1,1],AA[,m,1], labels=c(midline, right, left), cex=0.5)