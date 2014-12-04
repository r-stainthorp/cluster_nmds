# # # # # # # # # # # # # # # # # # # # # # # 
# NON METRIC MULTIDIMENSIONAL SCALING IN R  #
# # # # # # # # # # # # # # # # # # # # # # # 

# DEFINE WORKING DIRECTORY AND LOAD TABLE CONTAINING SPECIES ABUNDANCES

marsh<-read.csv("marsh_types.csv")

# LOAD PERMUTE AND VEGAN PACKAGES - FOR COMMUNITY ANALYSIS

library(permute)
library(vegan)

# EXTRACT MARSH TYPE COLUMN AND CLASSIFY AS A CATEGORICAL VARIABLE

type<-factor(marsh$Type)

# EXTRACT SPECIES ABUNDANCES BY REMOVING THE MARSH TYPE COLUMN

veg<-marsh[-1]

# PERFORM THE ORDINATION USING THE WRAPPER FUNCTION metaMDS:

vegmds<-metaMDS(veg,trymax=100)

# ONCE THE ORDINATION HAS COMPLETED (WHICH MAY TAKE SOME TIME) ASK R TO SHOW THE OUPPUT:

vegmds

# DISPLAY ORDINATION RESULTS, INITIALLY PLOTTING A BLANK FIGURE SO THAT THE SITES AND SPECIES POINTS CAN BE CUSTOMISED:

plot(vegmds, type = "n")

points(vegmds,display="sites", cex=0.8,pch=21,col="red",bg="yellow")
text(vegmds,display="spec",cex=0.7,col="blue")

# ALTERNATIVE CUSTOMISATION ARGUMENTS FOR THE POINTS PLOTTED ON THE GRAPH:
 
ordilabel(vegmds,dis="sites",cex=1.2,font=3,fill="hotpink",col="blue")
ordilabel(vegmds,dis="spec",cex=0.6,font=3,priority=colSums(veg),border=NA,col="blue")

# IF THE PLOT IS TOO CLUTTERED CAN 'ZOOM IN' ON PARTS OF THE PLOT:
plot(vegmds,xlim=c(-0.4,0.5),ylim=c(-0.4,0.5),type="n")
text(vegmds,display="sites",cex=0.7,col="blue")
points(vegmds,display="sites", cex=0.8,pch=20,col="red")
points(vegmds,display="spec", cex=0.8,pch="+",col="black")

cnam<-make.cepnames(names(veg))
stems<-colSums(veg)
orditorp(vegmds,dis="spec",labels=cnam,priority=stems,font=3,col="blue",pch="+",pcol="grey",air=0.6)

tnam<-make.cepnames(names(type))
orditorp(vegmds,dis="sites",labels=type,priority=stems,pch=20,pcol="yellow")

ordihull(vegmds,type,col="blue")
ordiellipse(vegmds,type,col=3,lwd=2)
ordispider(vegmds,type,col="red",label=TRUE)

marshenv<-read.csv("marsh_env.csv")
summary(marshenv)
ordfit<-envfit(vegmds~Elevation+Type,data=marshenv,perm=1000)
ordfit
plot(vegmds,xlim=c(-0.4,0.5),ylim=c(-0.4,0.5),type="n")
points(vegmds,display="sites", cex=0.5,pch=16,col="grey")
points(vegmds,display="spec", cex=0.8,pch="+",col="black")
plot(ordfit)

disfit<-envfit(vegmds~Distance+Type,data=marshenv,perm=1000)
disfit
plot(vegmds,xlim=c(-0.4,0.5),ylim=c(-0.4,0.5),type="n")
points(vegmds,display="sites", cex=0.5,pch=16,col="grey")
points(vegmds,display="spec", cex=0.8,pch="+",col="black")
plot(disfit)

help(ordicluster)

The vegan package has a special function to display the cluster fusions in ordination. The ordicluster function combines sites and cluster centroids similarly as the average linkage method:
R> ordiplot(ord, dis = "si") 
R> ordicluster("name of ordination e.g. vegMDS", "name of cluster analysis e.g. vegCluster")
