******************************************************************
*THESIS - SALTMARSH; CLUSTER									 *
*ANALYSIS, ANOVA, TUKEY'S TEST ETC. USING "VEGAN" AND "PERMUTE"  *
* PACKAGES.												   		 *
******************************************************************

# FIRST CHANGE WORKING DIRECTORY (UNDER MISC) TO FOLDER CONTAINING THE DATA FILE

# READ THE DATA FILE INTO R BY DRAGGING THE FILE NAME INTO THE FOLLOWING COMMAND
    # load the two data files
    # - one contains the vegetation coverage
    # - one contains the environmental variables
    # Both files are row labelled by Quadrat number
#  data file should be in the format: col.1 = quadID, col for each species and abundance in rows. keep the env data in a diff file with diff quadIDs so that r doesnt think they are the same thing


vegCover<-read.csv("cover_quad_data.csv")
vegEnv<-read.csv("environmental_quad_data.csv")

# look at the structure of the two data tables
str(vegCover)
str(vegEnv)

####################
# Cluster analysis #
####################

# turn the quadrat species data into a distance matrix, leaving out the quadrat numbers in column one

vegDist <- dist(vegCover[,-1])

# and now the clustering - using Ward's minimum variance

vegCluster <- hclust(vegDist, method="ward")


x<-as.dendrogram(z) 
plot(x, xlab="mydata complete-LINKAGE", ylim=c(0,4)) #visualization of the dendrogram 
clusters<-cutree(z, h=1.6) #obtain clusters at cutoff height=1.6 
ord<-cmdscale(d, k=2) #Multidimensional scaling of the data down to 2 dimensions 
clusplot(ord,clusters, color=TRUE, shade=TRUE,labels=4, lines=0) #visualization of the clusters in 2D map 
var1<-var(clusters==1) #variance of cluster 1 

#extract cluster memberships: 
clids = as.data.frame(clusters) 
names(clids) = c("id") 
clids$cdr = row.names(clids) 
row.names(clids) = c(1:dim(clids)[1]) 
clstructure = lapply(unique(clids$id), function(x){clids[clids$id == x,'cdr']}) 

clstructure[[1]] #get memberships of cluster 1 


# now, plot out the dendrogram, and decide on a height to divide the tree into clusters.

plot(vegCluster, hang=-1, xlab="Quadrat number")
cutHeight<-4000
rect.hclust(vegCluster, h=cutHeight)

# Now assign the membership of those clusters
vegEnv$Cluster <- cutree(vegCluster, h=cutHeight)

clids=as.data.frame(vegEnv$Cluster)
names(clids)=c("id")
clids$cdr=row.names(clids)
row.names(clids)=c(1:dim(clids)[1])
clstructure=lapply(unique(clid$id),function(x){clids[clids$id==x,'cdr']})

#  turn that from integers (1,2,3, etc.) into a categorical variable
vegEnv$Cluster <- factor(vegEnv$Cluster)

# and look at how many members there are of each cluster
table(vegEnv$Cluster)
    
# The order of the cluster numbers does _not_ relate to the order on the plot,
# so plot the clusters out again with the cluster numbers added to the labels
vegEnv$dendroLabels <- with(vegEnv, paste('Q ', Code, ': Cl ', Cluster, sep=''))
plot(vegCluster, labels=vegEnv$dendroLabels, hang=-1)
rect.hclust(vegCluster, h=cutHeight)
    
# put those clusters into the community data to find the abundances of each species
vegCover$Cluster <- vegEnv$Cluster
justVeg <- subset(vegCover, select=2:23)
# drop the quadrat labels and clusters
vegClusters <- aggregate(justVeg, by=list(Cluster=vegCover$Cluster), mean)
    
rownames(vegClusters) <- vegClusters$Cluster
print(t(vegClusters[, -1]), digits=2)

*********** LARGE DATA SETS - CUTTING DENDROGRAM AND VIEWING NODES ************

_________________

# extract relevant code from this: 
hc <- hclust(dist(USArrests), "ave")
(dend1 <- as.dendrogram(hc)) # "print()" method
str(dend1)          # "str()" method
str(dend1, max = 2, last.str= "'") # only the first two sub-levels
oo <- options(str.dendrogram.last = "\\") # yet another possibility
str(dend1, max = 2) # only the first two sub-levels
options(oo)# .. resetting them

# make plot window into two cols, two rows for graphs to be plotted as a matrix of four: 
op <- par(mfrow= c(2,2), mar = c(5,2,1,4))
plot(dend1)

## "triangle" type and show inner nodes:
# this code says: label the nodes with small circles [pch=c(1,NA) tells r not to label the inner nodes and leaves differently], make the shape triangle and centre the nodes
plot(dend1, nodePar=list(pch = c(1,NA), cex=0.8, lab.cex = 0.8),type = "t", center=TRUE)
# this code changes the look of the branches by making them dotted lines instead of solid:
plot(dend1, edgePar=list(col = 1:2, lty = 2:3),dLeaf=1, edge.root = TRUE)
# this code plot the dendrogram horizontally 
plot(dend1, nodePar=list(pch = 2:1,cex=.4*2:1, col = 2:3),horiz=TRUE)

## simple test for as.hclust() as the inverse of as.dendrogram():
stopifnot(identical(as.hclust(dend1)[1:4], hc[1:4]))

dend2 <- cut(dend1, h=4000)
plot(dend2$upper)
## leaves are wrong horizontally:
plot(dend2$upper, nodePar=list(pch = c(1,7), col = 2:1))
##  dend2$lower is *NOT* a dendrogram, but a list of .. :
plot(dend2$lower[[3]], nodePar=list(col=4))
## "inner" and "leaf" edges in different type & color :
plot(dend2$lower[[2]], nodePar=list(col=1),edgePar = list(lty=1:2, col=2:1), edge.root=TRUE)
par(op)
d3 <- dend2$lower[[2]][[2]][[1]]
plot(d3)
stopifnot(identical(d3, dend2$lower[[2]][[c(2,1)]]))
str(d3, last.str="'")

## merge() to join dendrograms:
(d13 <- merge(dend2$lower[[1]], dend2$lower[[3]]))
## merge() all parts back (using default 'height' instead of original one):
den.1 <- Reduce(merge, dend2$lower)
## or merge() all four parts at same height --> 4 branches (!)
d. <- merge(dend2$lower[[1]], dend2$lower[[2]], dend2$lower[[3]],
            dend2$lower[[4]])
## (with a warning) or the same using  do.call :
stopifnot(identical(d., do.call(merge, dend2$lower)))
plot(d., main="merge(d1, d2, d3, d4)  |->  dendrogram with a 4-split")

## "Zoom" in to the first dendrogram :
plot(dend1, xlim = c(1,20), ylim = c(1,50))

nP <- list(col=3:2, cex=c(2.0, 0.75), pch= 21:22,
           bg= c("light blue", "pink"),
           lab.cex = 0.75, lab.col = "tomato")
plot(d3, nodePar= nP, edgePar = list(col="gray", lwd=2), horiz = TRUE)

addE <- function(n) {
      if(!is.leaf(n)) {
        attr(n, "edgePar") <- list(p.col="plum")
        attr(n, "edgetext") <- paste(attr(n,"members"),"members")
      }
      n
}
d3e <- dendrapply(d3, addE)
plot(d3e, nodePar= nP)
plot(d3e, nodePar= nP, leaflab = "textlike")
--------------------------

# as.dendrogram converts class "hclust" objects in to class "dendrogram" 
dend1<-as.dendrogram(vegCluster)
plot(dend1)
dend2<-cut(dend1, h=4000)
# plots the 3rd cluster from beneath the cut
plot(dend2$upper) 





********************************************************************************

# code to extract the information from the cluster data 

str(unclass(vegCluster))

List of 7
 $ merge      : int [1:1777, 1:2] -3 -163 -170 -213 -6 -7 -131 -233 -295 -365 ...
 $ height     : num [1:1777] 0 0 0 0 0 0 0 0 0 0 ...
 $ order      : int [1:1778] 166 409 215 164 165 294 45 399 404 81 ...
 $ labels     : NULL
 $ method     : chr "ward"
 $ call       : language hclust(d = vegDist, method = "ward")
 $ dist.method: chr "euclidean"

order<-(vegCluster$order)
order

# code to try and select a node on the dendrogram and then extract information for that cluster (almost works) 

clus1<-identify(vegCluster)
str(clus1)


########################
# Analysis of Variance #
########################

# use a function to extract the cluster membership of each quadrat
vegEnv$Cluster <- cutree(vegCluster, h=4000)

#  turn that from integers (1,2,3, etc.) into a categorical variable
vegEnv$Cluster <- factor(vegEnv$Cluster)

# and look at how many members there are of each cluster
table(vegEnv$Cluster)

# now we can calculate analysis of variance using elevation
# to see if the elevation differs between clusters

elevMod <- aov(Elevation ~ Cluster, data=vegEnv)
summary(elevMod)

# same thing using distance from the estuary confluence

disMod <- aov(Distance ~ Cluster, data=vegEnv)
summary(disMod)

# Plotting a barplot with error bars in R
#  - this is a little bit awkward - but it is in Excel too!

# The first thing to do is to get the means 
# - the command tapply() breaks up a variable into groups for each cluster
#   and calculates a value within each group 

#MEANS
elevMean <- with(vegEnv, tapply(Elevation, Cluster, mean))
print(elevMean)

disMean <- with(vegEnv, tapply(Distance, Cluster, mean))
print(disMean)

# now create a function to calculate standard errors
# and do the tapply thing again 

se <- function(x) {sd(x)/sqrt(length(x))}
elevStErr <- with(vegEnv, tapply(Elevation, Cluster, se))

# now plot the barplot - increase the y axis length to fit in the error bars (ylim)

elevMids <- barplot(elevMean, ylab="Elevation (m)", xlab="Vegetation cluster", ylim=c(0,5))

# For mysterious reasons, there is no built in error bar command, 
# so we will use the arrows() command to draw them. 

arrows(x0=elevMids, y0=elevMean - elevStErr, x1=elevMids, 
       y1=elevMean + elevStErr, ang=90, code=3)

#BARPLOT FOR DISTANCE

# Plotting a barplot with error bars in R
#  - this is a little bit awkward - but it is in Excel too!

# The first thing to do is to get the means 
# - the command tapply() breaks up a variable into groups for each cluster
#   and calculates a value within each group 

disMean <- with(vegEnv, tapply(Distance, Cluster, mean))
print(disMean)

# now create a function to calculate standard errors
# and do the tapply thing again 

se <- function(x) {sd(x)/sqrt(length(x))}
disStErr <- with(vegEnv, tapply(Distance, Cluster, se))

# now plot the barplot - increase the y axis length to fit in the error bars (ylim)

disMids <- barplot(disMean, ylab="Distance from origin (m)", xlab="Vegetation cluster", ylim=c(0,7000))

# For mysterious reasons, there is no built in error bar command, 
# so we will use the arrows() command to draw them. 

arrows(x0=disMids, y0=disMean - disStErr, x1=disMids, 
       y1=disMean + disStErr, ang=90, code=3)

# Now - ANOVA tells us only that the clusters explain 
# significantly more variance than the null hypothesis.
# We don't know which clusters are different. We have to use
# Tukey's Honest Signficant Difference test to see which pairs
# of clusters are significantly different - the plot shows the
# observed difference between each pair and a confidence limit

# So first for depth
elevHSD <- TukeyHSD(elevMod)
plot(elevHSD)
abline(v=0, col='red')

# now distance from origin

disHSD <- TukeyHSD(disMod)
plot(disHSD)
abline(v=0, col='red')

#LINEAR REGRESSION
seaframe<-read.csv("Seacliff.csv")
str(seaframe)
Soil.depthmod<-lm(Soil.depth~Distance.from.origin, data=seaframe)
summary(Soil.depthmod)
plot(Soil.depth~Distance.from.origin, data=seaframe)
abline(Soil.depthmod)

seaframe<-read.csv("Seacliff.csv")
str(seaframe)
conductivitymod<-lm(log(conductivity)~Distance.from.origin, data=seaframe)
summary(conductivitymod)
plot(log(conductivity)~Distance.from.origin, data=seaframe)
abline(conductivitymod)

seaframe<-read.csv("Seacliff.csv")
str(seaframe)
pHmod<-lm(pH~Distance.from.origin, data=seaframe)
summary(conductivitymod)
plot(pH~Distance.from.origin, data=seaframe)
abline(conductivitymod)

Bare.ground  

seaframe<-read.csv("Seacliff.csv")
str(seaframe)
baregroundmod<-lm(Bare.ground~Distance.from.origin, data=seaframe)
summary(baregroundmod)
plot(Bare.ground~Distance.from.origin, data=seaframe)
abline(baregroundmod)

----------------------- 
# INFERENCES FROM COMMUNITY ANALYSES RELY ON WHETHER OR NOT THE COMMUNITY HAS BEEN ADEQUATELY SAMPLED. TO ASSESS THIS CAN CREATE SPECIES ACCUMULATION CURVES USING FUNCTION AVAILABLE IN THE "vegan" PACKAGE

# FIRST LOAD THE VEGAN PACKAGE (THIS REQUIRES THE permute PACKAGE TO BE LOADED FIRST):

library(permute)
library(vegan)

# create a matrix of only the species from the data file i.e.

# vegMat<-(vegCover[rows,cols])

# NOW USE THE specaccum FUNCTION TO GENERATE THE SPECIES CCUMULATION CURVE INFORMATION:

spp.curve<-specaccum(comm=vegMat,method="random",permutations=1000)
plot(spp.curve)

# HEATMAPS
# HEATMAPS TYPICALLY USED TO PRESENT MOLECULAR DATA (E.G. MICROARRAY) HOWEVER, COULD USE IT TO LOOK AT THE ENTIRETY OF THE DATASET BEFORE PROCEEDING TO ORDINATION TECHNIQUES

image(com)

----------------------
# HERE I WAS TRYING TO ADJUST THE AXIS SO THAT THEY WERE LABELLED WITH THE IND QUADRAT CODE AND THE SPECIES BUT I COULDN'T GET THEM TO FIT / CHANGE THE INTERVAL THAT THE TICK MARKS WERE LABELLED
quad=seq(1,1809,by=1)
species=seq(1,22,by=1)
image(x=quad,y=species,z=com,xaxt="n",yaxt="n")
axis(side=1,tick=TRUE,at=quad,labels=code)
axis(side=2,tick=TRUE,at=species,labels=colnames(com))
par(las=2,cex.axis=0.5)
-----------------------
