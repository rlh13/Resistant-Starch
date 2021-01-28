library(qiime2R)
library(phyloseq)
library(microbiome)
library(vegan)
#Read in your distance matrix from QIIME
DM <-read_qza("/Users/Riley/Box Sync/Riley's Documents/Keim Lab/Resistant Starch/RS Data/Microbiome/DADA2_350_unfiltered/diversity_DADA2_nomock/weighted_unifrac_distance_matrix.qza")
DM <- DM$data
attr(DM,"Labels")
map <- read.csv("/Users/Riley/Box Sync/Riley's Documents/Keim Lab/Resistant Starch/RS Data/Microbiome/FL103_Mapfile.csv",sep=",",header=T,stringsAsFactors=F)
map.sub <- map[map$X.SampleID %in%  attr(DM, "Labels"), ] #remove samples not in distance matrix

all(attr(DM, "Labels") == map.sub$X.SampleID) #check to see if map and dist match
ord = match(attr(DM, "Labels"), map.sub$X.SampleID) # make a vector of the order of the samples in your distance matrix
map1 = map.sub[ord,] # re-order the mapping file to match the distance matrix
all(map1$X.SampleID == attr(DM, "Labels")) # check to make sure it worked

# #NA values cause adonis to error
# map2 <- map1[!is.na(map1$diff_time_hrs),]

DM <- as.data.frame(as.matrix(DM))
#DM <- DM[!is.na(map1$diff_time_hrs),!is.na(map1$diff_time_hrs)]

all(map1$X.SampleID == rownames(DM)) # check to make sure it worked

#Graph the results
#Perform Principal Coordinate analysis
PCo <- cmdscale(DM, k=16, eig = T)
all(map1$X.SampleID == row.names(PCo$points)) #make sure everything stayed in order

p_Trt <- ordiplot(PCo, choices = c(1,2), type = "none", display = "sites", 
                  cex = 0.5)

points(p_Trt, what = "sites", select = c(map1$Trt == "RS"),
       col = "deepskyblue3", cex = 0.5, pch = 1)
points(p_Trt, what = "sites", select = c(map1$Trt == "C"),
       col = "forestgreen", cex = 0.5, pch = 1)
points(p_Trt, what = "sites", select = c(map1$Trt == "base"),
       col = "firebrick1", cex = 0.5, pch = 1)
points(p_Trt, what = "sites", select = c(map1$Trt == "none"),
       col = "darkorchid4", cex = 0.5, pch = 1)

ordispider(p_Trt, groups = map1$Trt , display="sites", 
           spiders = c("median"),
           show.groups = c("RS"),
           col = c("deepskyblue3"))
ordispider(p_Trt, groups = map1$Trt , display="sites", 
           spiders = c("median"),
           show.groups = c("C"),
           col = c("forestgreen"))
ordispider(p_Trt, groups = map1$Trt , display="sites", 
           spiders = c("median"),
           show.groups = c("base"),
           col = c("firebrick1"))
ordispider(p_Trt, groups = map1$Trt , display="sites", 
           spiders = c("median"),
           show.groups = c("none"),
           col = c("darkorchid4"))

legend("topright", legend = c("RS","Control","Baseline", "Washout"), 
       col = c("deepskyblue3", "forestgreen", "firebrick1", "darkorchid4"),
       lty = 1, lwd = 10, cex = 0.5)

#Repeat for Trt2
p_Trt2 <- ordiplot(PCo, choices = c(1,2), type = "none", display = "sites", 
                   cex = 0.5)

points(p_Trt2, what = "sites", select = c(map1$Trt2 == "RS"),
       col = "deepskyblue3", cex = 0.5, pch = 1)
points(p_Trt2, what = "sites", select = c(map1$Trt2 == "C"),
       col = "forestgreen", cex = 0.5, pch = 1)
points(p_Trt2, what = "sites", select = c(map1$Trt2 == "before"),
       col = "firebrick1", cex = 0.5, pch = 1)

ordispider(p_Trt2, groups = map1$Trt2 , display="sites", 
           spiders = c("median"),
           show.groups = c("RS"),
           col = c("deepskyblue3"))
ordispider(p_Trt2, groups = map1$Trt2 , display="sites", 
           spiders = c("median"),
           show.groups = c("C"),
           col = c("forestgreen"))
ordispider(p_Trt2, groups = map1$Trt2 , display="sites", 
           spiders = c("median"),
           show.groups = c("before"),
           col = c("firebrick1"))

legend("topright", legend = c("RS","Control","Before"), 
       col = c("deepskyblue3", "forestgreen", "firebrick1"),
       lty = 1, lwd = 10, cex = 0.5)

#Repeat for Trt3
pdf("/Users/Riley/Box Sync/Riley's Documents/Keim Lab/Resistant Starch/RS Data/Microbiome/DADA2_350_unfiltered/Figures/Spider_Trt3.pdf", width=10, height=10, pointsize = 19)
p_Trt3 <- ordiplot(PCo, choices = c(1,2), type = "none", display = "sites", 
                  cex.axis=1.5,cex.lab=1.5,cex=0.5)

points(p_Trt3, what = "sites", select = c(map1$Trt3 == "RS"),
       col = "deepskyblue3", cex = 0.5, pch = 1)
points(p_Trt3, what = "sites", select = c(map1$Trt3 == "C"),
       col = "forestgreen", cex = 0.5, pch = 1)
points(p_Trt3, what = "sites", select = c(map1$Trt3 == "Pre_C"),
       col = "firebrick1", cex = 0.5, pch = 1)
points(p_Trt3, what = "sites", select = c(map1$Trt3 == "Pre_RS"),
       col = "darkorchid4", cex = 0.5, pch = 1)

ordispider(p_Trt3, groups = map1$Trt3 , display="sites", 
           spiders = c("median"),
           show.groups = c("RS"),
           col = c("deepskyblue3"))
ordispider(p_Trt3, groups = map1$Trt3 , display="sites", 
           spiders = c("median"),
           show.groups = c("C"),
           col = c("forestgreen"))
ordispider(p_Trt3, groups = map1$Trt3 , display="sites", 
           spiders = c("median"),
           show.groups = c("Pre_C"),
           col = c("firebrick1"))
ordispider(p_Trt3, groups = map1$Trt3 , display="sites", 
           spiders = c("median"),
           show.groups = c("Pre_RS"),
           col = c("darkorchid4"))

legend("topright", legend = c("RS","Control","Pre Control", "Pre RS"), 
       col = c("deepskyblue3", "forestgreen", "firebrick1", "darkorchid4"),
       lty = 1, lwd = 10, cex = 1.25)
dev.off()

#Arrow plot
points <- data.frame(PCo$points[,1:2]) #get points for first 2 principal coordinates
#How to get variance explained by each PC in PCo?
eig <- PCo$eig
eig[eig < 0] = 0 #QIIME2 makes all negative eigenvalues zero
eig[1]/sum(eig) #0.3595839 = 35.96%
eig[2]/sum(eig) #0.1002939 = 10.03% 

plot(points, type="n", asp=1, xlab="PC1 (35.96%)", ylab="PC2 (10.03%)", cex.axis=1.5, cex.lab=2)
#Connect Pre_C and Pre_RS (FL.103.XXX.1 and FL.103.XXX.3) by dotted lines
id <- unique(map1$SubID)
for(i in 1:length(id)){
  x1 <- points[paste0(id[i],".1"),"X1"] #get PC1 coordinate of FL.103.XXX.1
  x2 <- points[paste0(id[i],".3"),"X1"]  #get PC1 coordinate of FL.103.XXX.3
  y1 <- points[paste0(id[i],".1"),"X2"]  #get PC2 coordinate of FL.103.XXX.1
  y2 <- points[paste0(id[i],".3"),"X2"] #get PC2 coordinate of FL.103.XXX.3
  segments(x1,y1,x2,y2,lty=3,lwd=2) #draw segment between FL.103.XXX.1 and FL.103.XXX.3
}
#Connect Pre_C and C by dotted lines
preC_id <- map1[map$Trt3 %in% c("Pre_C","C"),]$X.SampleID
for(i in seq(1,(length(preC_id)-1),2)){ #go through preC_ID by increments of 2
  x1 <- points[preC_id[i],"X1"] #get PC1 coordinate of Pre_C
  x2 <- points[preC_id[i+1],"X1"]  #get PC1 coordinate of C
  y1 <- points[preC_id[i],"X2"]  #get PC2 coordinate of Pre_C
  y2 <- points[preC_id[i+1],"X2"] #get PC2 coordinate of C
  arrows(x1,y1,x2,y2,length="0.15",angle="20",code="2", col="forestgreen", lty=1, lwd=2.5) #draw segment between Pre_C and C
}

preRS_id <- map1[map$Trt3 %in% c("Pre_RS","RS"),]$X.SampleID
preRS_id <- head(preRS_id,-1) #remove last sample (FL.103.058.4) that has no Pre_RS correlate
for(i in seq(1,(length(preRS_id)-1),2)){ #go through preC_ID by increments of 2
  x1 <- points[preRS_id[i],"X1"] #get PC1 coordinate of Pre_RS
  x2 <- points[preRS_id[i+1],"X1"]  #get PC1 coordinate of RS
  y1 <- points[preRS_id[i],"X2"]  #get PC2 coordinate of Pre_RS
  y2 <- points[preRS_id[i+1],"X2"] #get PC2 coordinate of RS
  arrows(x1,y1,x2,y2,length="0.15",angle="20",code="2", col="deepskyblue3", lty=1, lwd=2.5) #draw segment between Pre_C and C
}
