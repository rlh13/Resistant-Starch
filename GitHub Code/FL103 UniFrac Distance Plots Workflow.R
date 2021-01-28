library(ggplot2)
#Violin plots of UniFrac Distances
#read in treatment distance matrices
trt_DM <- read.csv("/Users/Riley/Box Sync/Riley's Documents/Keim Lab/Resistant Starch/RS Data/Microbiome/DADA2_350_unfiltered/diversity_DADA2_nomock/weighted-unifrac-treatment-matrix.csv",sep=",",header=TRUE)
trt2_DM <- read.csv("/Users/Riley/Box Sync/Riley's Documents/Keim Lab/Resistant Starch/RS Data/Microbiome/DADA2_350_unfiltered/diversity_DADA2_nomock/weighted-unifrac-treatment2-matrix.csv",sep=",",header=TRUE)
trt3_DM <- read.csv("/Users/Riley/Box Sync/Riley's Documents/Keim Lab/Resistant Starch/RS Data/Microbiome/DADA2_350_unfiltered/diversity_DADA2_nomock/weighted-unifrac-treatment3-matrix.csv",sep=",",header=TRUE)
#Group 1 column is the "distance from" variable
#Get intra-individual distances only
#First 10 characters of Subject ID are FL.103.XXX
intra_trt_DM <- trt_DM[substr(trt_DM$SubjectID1,1,10)==substr(trt_DM$SubjectID2,1,10),]
intra_trt2_DM <- trt2_DM[substr(trt2_DM$SubjectID1,1,10)==substr(trt2_DM$SubjectID2,1,10),]
intra_trt3_DM <- trt3_DM[substr(trt3_DM$SubjectID1,1,10)==substr(trt3_DM$SubjectID2,1,10),]

#Get inter-individual distances only
inter_trt_DM <- trt_DM[substr(trt_DM$SubjectID1,1,10)!=substr(trt_DM$SubjectID2,1,10),]
inter_trt2_DM <- trt2_DM[substr(trt2_DM$SubjectID1,1,10)!=substr(trt2_DM$SubjectID2,1,10),]
inter_trt3_DM <- trt3_DM[substr(trt3_DM$SubjectID1,1,10)!=substr(trt3_DM$SubjectID2,1,10),]

#Plot intra-individual distances for Trt2 (before, Control, RS)
intra_trt2_plot <- ggplot(intra_trt2_DM,aes(x=factor(Group1,level=c("before","C","RS")),y=Distance,fill=Group1))+
  geom_violin(postition=position_dodge(1)) + 
  stat_summary(fun.y=mean, geom="point",position=position_dodge(1)) +
  scale_fill_manual(values=c("firebrick1", "forestgreen","deepskyblue3")) +
  labs(title="Within Subject",x="Treatment Group",y="Weighted UniFrac Distance") +
  scale_x_discrete(labels=c("Before","Control","RS")) +
  theme(plot.title = element_text(size=25, hjust = 0.5),
        axis.text=element_text(size=15), axis.title=element_text(size=20),
        legend.position="none",panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(colour="black"),axis.text.y = element_text(colour="black"))

#Plot inter-individual distances for Trt (base, C, none, RS) and Trt3 (Pre_C, C, Pre_RS, RS)
inter_trt3_plot <- ggplot(intra_trt3_DM,aes(x=factor(Group1,level=c("Pre_C","C","Pre_RS","RS")),y=Distance,fill=Group1))+
  geom_violin(postition=position_dodge(1)) + 
  stat_summary(fun.y=mean, geom="point",position=position_dodge(1)) +
  scale_fill_manual(values=c("forestgreen","firebrick1","darkorchid4","deepskyblue3")) +
  labs(title="Between Subjects",x="Treatment Group",y="Weighted UniFrac Distance") +
  scale_x_discrete(labels=c("Pre-Control","Control","Pre-RS","RS")) +
  theme(plot.title = element_text(size=25, hjust = 0.5),
        axis.text=element_text(size=15), axis.title=element_text(size=20),
        legend.position="none",panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(colour="black"),axis.text.y = element_text(colour="black"))
ggarrange(intra_trt2_plot,inter_trt3_plot,ncol=2)

#Trt doesn't look much different than Trt3 (use Trt3)
inter_trt_plot <- ggplot(intra_trt_DM,aes(x=factor(Group1,level=c("base","C","none","RS")),y=Distance,fill=Group1))+
  geom_violin(postition=position_dodge(1)) + 
  stat_summary(fun.y=mean, geom="point",position=position_dodge(1)) +
  scale_fill_manual(values=c("firebrick1", "forestgreen","darkorchid4","deepskyblue3")) +
  labs(title="Between Subjects",x="Treatment Group",y="Weighted UniFrac Distance") +
  scale_x_discrete(labels=c("Baseline","Control","Washout","RS")) +
  theme(plot.title = element_text(size=25, hjust = 0.5),
        axis.text=element_text(size=15), axis.title=element_text(size=20),
        legend.position="none",panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(colour="black"),axis.text.y = element_text(colour="black"))

#Signficant differences?
#Perform Wilcox post-test
intra_trt2_before_C <- intra_trt2_DM[intra_trt2_DM$Group1=="before" | intra_trt2_DM$Group1=="C",]
intra_trt2_before_RS <- intra_trt2_DM[intra_trt2_DM$Group1=="before" | intra_trt2_DM$Group1=="RS",]
intra_trt2_C_RS <- intra_trt2_DM[intra_trt2_DM$Group1=="C" | intra_trt2_DM$Group1=="RS",]
wilcox.test(Distance ~ Group1, data=intra_trt2_before_C, paired=FALSE, alt="two.sided", conf.int=TRUE, correct=FALSE)
#p-value = 0.9594
wilcox.test(Distance ~ Group1, data=intra_trt2_before_RS, paired=FALSE, alt="two.sided", conf.int=TRUE, correct=FALSE)
#p-value = 0.09745
wilcox.test(Distance ~ Group1, data=intra_trt2_C_RS, paired=FALSE, alt="two.sided", conf.int=TRUE, correct=FALSE)
#p-value = 0.1406
p<-c(0.9594, 0.09745, 0.1406)
p.adjust(p, method="fdr")
#0.9594 0.2109 0.2109

inter_trt3_pre_pre <- inter_trt3_DM[inter_trt3_DM$Group1=="Pre_C" | inter_trt3_DM$Group1=="Pre_RS",]
inter_trt3_pre_C <- inter_trt3_DM[inter_trt3_DM$Group1=="Pre_C" | inter_trt3_DM$Group1=="C",]
inter_trt3_pre_RS <- inter_trt3_DM[inter_trt3_DM$Group1=="Pre_RS" | inter_trt3_DM$Group1=="RS",]
inter_trt3_C_RS <- inter_trt3_DM[inter_trt3_DM$Group1=="C" | inter_trt3_DM$Group1=="RS",]
wilcox.test(Distance ~ Group1, data=inter_trt3_pre_pre, paired=FALSE, alt="two.sided", conf.int=TRUE, correct=FALSE)
#p-value = 0.006066
wilcox.test(Distance ~ Group1, data=inter_trt3_pre_C, paired=FALSE, alt="two.sided", conf.int=TRUE, correct=FALSE)
#p-value = 0.6442
wilcox.test(Distance ~ Group1, data=inter_trt3_pre_RS, paired=FALSE, alt="two.sided", conf.int=TRUE, correct=FALSE)
#p-value = 0.00000001413
wilcox.test(Distance ~ Group1, data=inter_trt3_C_RS, paired=FALSE, alt="two.sided", conf.int=TRUE, correct=FALSE)
#p-value = 0.0006174
p<-c(0.006066, 0.6442, 0.00000001413,0.0006174)
p.adjust(p, method="fdr")
#0.00808800000 0.64420000000 0.00000005652 0.00123480000
