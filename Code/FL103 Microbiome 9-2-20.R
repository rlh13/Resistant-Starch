#Load metadata and otu/tax tables from the above at genus, family, and phylum levels (relative abundance)
metadata <- read.csv("/Users/Riley/Box Sync/Riley's Documents/Keim Lab/Resistant Starch/RS Data/Microbiome/FL103_Mapfile.csv",sep=",",header=T,stringsAsFactors=F)
otu_w_tax_genusrel <- read.csv("/Users/Riley/Box Sync/Riley's Documents/Keim Lab/Resistant Starch/RS Data/Microbiome/DADA2_350_unfiltered/otu_w_tax_genusrel.csv",sep=",",header = T) 
otu_w_tax_genusrel$X <- NULL #remove X column
otu_w_tax_famrel <- read.csv("/Users/Riley/Box Sync/Riley's Documents/Keim Lab/Resistant Starch/RS Data/Microbiome/DADA2_350_unfiltered/otu_w_tax_famrel.csv",sep=",",header = T)
otu_w_tax_famrel$X <- NULL
otu_w_tax_phyrel <- read.csv("/Users/Riley/Box Sync/Riley's Documents/Keim Lab/Resistant Starch/RS Data/Microbiome/DADA2_350_unfiltered/otu_w_tax_phyrel.csv",sep=",",header = T)
otu_w_tax_phyrel$X <- NULL

#Define taxa of interest at genus, family, and phylum level
gen_taxa <- c("g__Bifidobacterium","g__Bacteroides","g__Prevotella","g__Faecalibacterium","g__Ruminococcus","g__[Eubacterium]","g__Gemmiger","g__Roseburia")
fam_taxa <- c("f__Bifidobacteriaceae", "f__Bacteroidaceae","f__Prevotellaceae", "f__Ruminococcaceae", "f__Eubacteriaceae")
phy_taxa <- c("p__Firmicutes","p__Bacteroidetes","p__Actinobacteria")
tax_names <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")

#Get taxa of interest and add to metadata
genus_table <- otu_w_tax_genusrel[which(otu_w_tax_genusrel$Genus %in% gen_taxa), ]
#Get rid of g__Ruminococcus belonging to f__Lachnospiraceae
genus_table <- genus_table[-which(genus_table$Genus %in% c("g__Ruminococcus") & genus_table$Family %in% c("f__Lachnospiraceae")), ]
fam_table <- otu_w_tax_famrel[which(otu_w_tax_famrel$Family %in% fam_taxa), ]
phy_table <- otu_w_tax_phyrel[which(otu_w_tax_phyrel$Phylum %in% phy_taxa), ]

rownames(genus_table) <- genus_table$Genus
genus_table <- genus_table[ ,-which(names(genus_table) %in% tax_names)] #remove tax columns
genus_table <- t(genus_table)
rownames(fam_table) <- fam_table$Family
fam_table <- fam_table[ ,-which(names(fam_table) %in% tax_names)] #remove tax columns
fam_table <- t(fam_table)
rownames(phy_table) <- phy_table$Phylum
phy_table <- phy_table[ ,-which(names(phy_table) %in% tax_names)] #remove tax columns
phy_table <- t(phy_table)

metadata_w_tax <- cbind(metadata,genus_table,fam_table,phy_table)
#Remove the mock community
metadata_w_tax <- metadata_w_tax[-which(metadata_w_tax$X.SampleID %in% "MockComm"), ]

#Calculate Shannon diversity (use relative abundance) and Chao1 richness (use counts) and add to metadata
#use non-glommed data to calculate
library(vegan)
library(fossil)
otu_w_noglomtax_rel <- read.csv("/Users/Riley/Box Sync/Riley's Documents/Keim Lab/Resistant Starch/RS Data/Microbiome/DADA2_350_unfiltered/otu_w_noglomtax_rel.csv",sep=",",header = T)
otu_w_noglomtax_rel$X <- NULL
otu_w_noglomtax <- read.csv("/Users/Riley/Box Sync/Riley's Documents/Keim Lab/Resistant Starch/RS Data/Microbiome/DADA2_350_unfiltered/otu_w_noglomtax.csv",sep=",",header = T)
otu_w_noglomtax$X <- NULL
shannon <- diversity(t(otu_w_noglomtax_rel[,1:121]), index = "shannon")
metadata_w_tax$Shannon_div <- head(shannon,-1) #removes mock community sample

chao <- as.matrix(apply(otu_w_noglomtax[,1:121],2,chao1))#exclude columns with taxonomy
colnames(chao) <- "Chao1" 
metadata_w_tax$Chao1 <- head(chao,-1) #removes mock community sample

#Violin plots
library(ggplot2)
library(reshape2)
meta_gen_melt <- metadata_w_tax[,c("X.SampleID","Trt3",gen_taxa)]
meta_gen_melt <- melt(meta_gen_melt)
ggplot(meta_gen_melt,aes(x=factor(Trt3,level=c("Pre_C","C","Pre_RS","RS")),y=value,fill=variable))+
  geom_violin(postition=position_dodge(1)) + 
  stat_summary(fun.y=mean, geom="point",position=position_dodge(1)) +
  labs(x="Treatment",y="Relative Abundance") +
  theme(plot.title = element_text(size=25, hjust = 0.5),
        axis.text=element_text(size=15), axis.title=element_text(size=20),
        legend.position="none",panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(colour="black"),axis.text.y = element_text(colour="black")) +
  facet_wrap(~variable)
#only genera that showed up in DESeq2
deseq_genus <- c("g__Bifidobacterium","g__Bacteroides","g__Faecalibacterium","g__Ruminococcus","g__Gemmiger","g__Roseburia")
meta_dsqgen_melt <- metadata_w_tax[,c("X.SampleID","Trt3",deseq_genus)]
meta_dsqgen_melt <- melt(meta_dsqgen_melt)
ggplot(meta_dsqgen_melt,aes(x=factor(Trt3,level=c("Pre_C","C","Pre_RS","RS")),y=value,fill=variable))+
  geom_violin(postition=position_dodge(1)) + 
  stat_summary(fun.y=mean, geom="point",position=position_dodge(1)) +
  labs(x="Treatment",y="Relative Abundance") +
  theme(plot.title = element_text(size=25, hjust = 0.5),
        axis.text=element_text(size=15), axis.title=element_text(size=20),
        legend.position="none",panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(colour="black"),axis.text.y = element_text(colour="black")) +
  facet_wrap(~variable)

meta_fam_melt <- metadata_w_tax[,c("X.SampleID","Trt3",fam_taxa)]
meta_fam_melt <- melt(meta_fam_melt)
ggplot(meta_fam_melt,aes(x=factor(Trt2,level=c("Pre_C","C","Pre_RS","RS")),y=value,fill=variable))+
  geom_violin(postition=position_dodge(1)) + 
  stat_summary(fun.y=mean, geom="point",position=position_dodge(1)) +
  labs(x="Treatment",y="Relative Abundance") +
  theme(plot.title = element_text(size=25, hjust = 0.5),
        axis.text=element_text(size=15), axis.title=element_text(size=20),
        legend.position="none",panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(colour="black"),axis.text.y = element_text(colour="black")) +
  facet_wrap(~variable)

meta_phy_melt <- metadata_w_tax[,c("X.SampleID","Trt3",phy_taxa)]
meta_phy_melt <- melt(meta_phy_melt)
ggplot(meta_phy_melt,aes(x=factor(Trt3,level=c("Pre_C","C","Pre_RS","RS")),y=value,fill=variable))+
  geom_violin(postition=position_dodge(1)) + 
  stat_summary(fun.y=mean, geom="point",position=position_dodge(1)) +
  labs(x="Treatment",y="Relative Abundance") +
  theme(plot.title = element_text(size=25, hjust = 0.5),
        axis.text=element_text(size=15), axis.title=element_text(size=20),
        legend.position="none",panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(colour="black"),axis.text.y = element_text(colour="black")) +
  facet_wrap(~variable)

meta_div_melt <- metadata_w_tax[,c("X.SampleID","Trt3",c("Shannon_div","Chao1"))]
meta_div_melt <- melt(meta_div_melt)
div.labs <- c("Shannon Diversity", "Chao1 Richness")
names(div.labs) <- c("Shannon_div", "Chao1")

pdf("/Users/Riley/Box Sync/Riley's Documents/Keim Lab/Resistant Starch/RS Data/Microbiome/DADA2_350_unfiltered/Figures/Violin_shannon_chao.pdf", width=24, height=10)
ggplot(meta_div_melt,aes(x=factor(Trt3,level=c("Pre_C","C","Pre_RS","RS")),y=value,fill=variable))+
  geom_violin(postition=position_dodge(1)) + 
  stat_summary(fun.y=mean, geom="point",position=position_dodge(1)) +
  labs(x="Treatment",y="") +
  theme(plot.title = element_text(size=70, hjust = 0.5),
        axis.text=element_text(size=50), axis.title=element_text(size=60),
        legend.position="none",panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(colour="black"),axis.text.y = element_text(colour="black"),
        strip.text = element_text(size=50)) +
  facet_wrap(~variable,scales = "free_y",labeller = labeller(variable=div.labs))
dev.off()

#P/B ratio calculation
#using log2 transformed values according to Roager et al (2014)
PB_ratio <- log(metadata_w_tax$g__Prevotella/metadata_w_tax$g__Bacteroides,2)
PB_ratio <- data.frame(PB_ratio=PB_ratio[!is.infinite(PB_ratio)&!is.na(PB_ratio)])
mean(PB_ratio$PB_ratio) #-2.620479
ggplot(PB_ratio, aes(x=PB_ratio)) + 
  geom_histogram(colour="black", fill="white",binwidth = 1.75) + 
  #geom_density(alpha=.2, fill="#FF6666") +
  geom_vline(aes(xintercept=mean(PB_ratio)), color="blue", linetype="dashed", size=1) +
  labs(title="P/B ratio histogram",x="P/B ratio", y = "Count") +
  theme(plot.title = element_text(size=25, hjust = 0.5),
        axis.text=element_text(size=15), axis.title=element_text(size=20),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(colour="black"),axis.text.y = element_text(colour="black"))
PB_table <- cbind(X.SampleID=metadata_w_tax$X.SampleID,PB_ratio)
samps <- data.frame(PB_table[!is.infinite(PB_ratio) & !is.na(PB_ratio),])$X.SampleID
#Too many values are NA or infinite (only 16 numeric values), try at family level
metadata_w_tax$PB_ratio_fam <- log(metadata_w_tax$f__Prevotellaceae/metadata_w_tax$f__Bacteroidaceae,2)
#Values are the same, double checked otu tables and found that all Bacteroidaceae and Prevotellaceae were classified at genus level (no additional family level only classification)
ggplot(metadata_w_tax[metadata_w_tax$X.SampleID %in% samps,],aes(x=log(g__Bacteroides,2),y=log(g__Prevotella,2),color=SubID))+
  geom_point(size=5,aes(shape=Trt)) + 
  labs(x="log(Bacteroides)",y="log(Prevotella)") +
  theme(plot.title = element_text(size=25, hjust = 0.5),
        axis.text=element_text(size=15), axis.title=element_text(size=20),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(colour="black"),axis.text.y = element_text(colour="black"))
#F/B ratio
metadata_w_tax$FB_ratio <- metadata_w_tax$p__Firmicutes/metadata_w_tax$p__Bacteroidetes
ggplot(metadata_w_tax,aes(x=p__Bacteroidetes,y=p__Firmicutes))+
  geom_point() + 
  labs(x="Bacteroidetes (%)",y="Firmicutes (%)") +
  theme(plot.title = element_text(size=25, hjust = 0.5),
        axis.text=element_text(size=15), axis.title=element_text(size=20),
        legend.position="none",panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(colour="black"),axis.text.y = element_text(colour="black"))
ggplot(metadata_w_tax,aes(x=Trt3,y=FB_ratio))+
  geom_violin(postition=position_dodge(1),fill="#F1C40F") + 
  stat_summary(fun.y=mean, geom="point",position=position_dodge(1)) +
  labs(x="Treatment",y="Firmicutes/Bacteroidetes\nRatio") +
  theme(plot.title = element_text(size=25, hjust = 0.5),
        axis.text=element_text(size=15), axis.title=element_text(size=20),
        legend.position="none",panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(colour="black"),axis.text.y = element_text(colour="black"))


#Do abundances of taxa correlate with relevant health outcomes?
library(Hmisc)
library(dplyr)
library(corrplot)
library(Deducer)
library(ggpubr)
library(grid)
library(gridExtra)
library(purrr)
#FB ratio correlation with BMI at baseline?
metadata_base <- metadata_w_tax[(metadata_w_tax$Trt=="base"),]
rcorr(metadata_base$FB_ratio,metadata_base$BMI,type="pearson")$P #0.8737361, no correlation

#Get all taxa (not just those from DESeq)
#Get last filled tax level and make row names
otu_w_tax_id <- read.csv("/Users/Riley/Box Sync/Riley's Documents/Keim Lab/Resistant Starch/RS Data/Microbiome/DADA2_350_unfiltered/otu_w_tax_genusrel_OTUID.csv",sep=",",header = T) 
otu_w_tax_id$X <- NULL #remove X column
tax_names <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
tax_level <- apply(otu_w_tax_id[,tax_names], 1, function(x) tail(na.omit(x), 1))
otu_w_tax_id <- cbind(otu_w_tax_id, LastTax = tax_level)
rownames(otu_w_tax_id) <- paste(otu_w_tax_id$LastTax,otu_w_tax_id$Row.names,sep="_")
otu_w_tax_id[,c(tax_names,"Row.names","LastTax","MockComm")] <- NULL #remove non-numeric columns and Mock Community sample
all_otu_w_tax_id <- t(otu_w_tax_id)
#instances of [ or - mess up the function so have to remove these
colnames(all_otu_w_tax_id) <- gsub("\\[|\\]",replacement="",colnames(all_otu_w_tax_id))
colnames(all_otu_w_tax_id) <- gsub("-",replacement="",colnames(all_otu_w_tax_id))
colnames(all_otu_w_tax_id) <- paste(colnames(all_otu_w_tax_id),"per",sep="_")
metadata_w_alltax <- cbind(head(metadata,-1),all_otu_w_tax_id)
#with absolute abundance
otu_w_abstax_id <- read.csv("/Users/Riley/Box Sync/Riley's Documents/Keim Lab/Resistant Starch/RS Data/Microbiome/DADA2_350_unfiltered/otu_w_tax_genus_OTUID.csv",sep=",",header = T) 
otu_w_abstax_id$X <- NULL #remove X column
tax_names <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
tax_level <- apply(otu_w_abstax_id[,tax_names], 1, function(x) tail(na.omit(x), 1))
otu_w_abstax_id <- cbind(otu_w_abstax_id, LastTax = tax_level)
rownames(otu_w_abstax_id) <- paste(otu_w_abstax_id$LastTax,otu_w_abstax_id$Row.names,sep="_")
otu_w_abstax_id[,c(tax_names,"Row.names","LastTax","MockComm")] <- NULL #remove non-numeric columns and Mock Community sample
all_otu_w_abstax_id <- t(otu_w_abstax_id)
#instances of [ or - mess up the function so have to remove these
colnames(all_otu_w_abstax_id) <- gsub("\\[|\\]",replacement="",colnames(all_otu_w_abstax_id))
colnames(all_otu_w_abstax_id) <- gsub("-",replacement="",colnames(all_otu_w_abstax_id))
metadata_w_alltax <- cbind(metadata_w_alltax,all_otu_w_abstax_id)
write.csv(metadata_w_alltax,"/Users/Riley/Box Sync/Riley's Documents/Keim Lab/Resistant Starch/RS Data/Microbiome/DADA2_350_unfiltered/metadata_w_alltax_absrel.csv")
#get top 20 most abundant OTUs
means <- colMeans(all_otu_w_tax_id,na.rm=T)
top20_otu_w_tax_id <- all_otu_w_tax_id[,order(means)]
top20_otu_w_tax_id <- subset(top20_otu_w_tax_id,select = -c(1:(ncol(top20_otu_w_tax_id)-20)))
top20_tax <- colnames(top20_otu_w_tax_id)
metadata_w_top20tax <- cbind(head(metadata,-1),top20_otu_w_tax_id)
write.csv(metadata_w_top20tax,"/Users/Riley/Box Sync/Riley's Documents/Keim Lab/Resistant Starch/RS Data/Microbiome/DADA2_350_unfiltered/metadata_w_top20tax.csv")

#subset to only baseline taxa and replicate each value
base_samps <- metadata[which(metadata_w_alltax$Trt=="base"),]$X.SampleID
otu_w_tax_base <- all_otu_w_tax_id[which(rownames(all_otu_w_tax_id) %in% base_samps),]
#Subset to only top 20 taxa 
base_means <- colMeans(otu_w_tax_base,na.rm=T)
top20_otu_w_tax_base <- otu_w_tax_base[,order(base_means)]
top20_otu_w_tax_base <- subset(top20_otu_w_tax_base,select = -c(1:(ncol(top20_otu_w_tax_base)-20)))
#instances of [ or - mess up the function so have to remove these
colnames(top20_otu_w_tax_base) <- gsub("\\[|\\]",replacement="",colnames(top20_otu_w_tax_base))
colnames(top20_otu_w_tax_base) <- gsub("-",replacement="",colnames(top20_otu_w_tax_base))
top20_tax_base <- colnames(top20_otu_w_tax_base)
top20_otu_w_tax_repbase <- as.data.frame(apply(top20_otu_w_tax_base,2,rep,each=4))
metadata_w_top20tax_base <- cbind(head(metadata,-1),top20_otu_w_tax_repbase)
write.csv(metadata_w_top20tax_base,"/Users/Riley/Box Sync/Riley's Documents/Keim Lab/Resistant Starch/RS Data/Microbiome/DADA2_350_unfiltered/metadata_w_top20tax_base.csv")

#subset to Pre-RS taxa and replicate each value
preRS_samps <- metadata[which(metadata_w_alltax$Trt3=="Pre_RS"),]$X.SampleID
otu_w_tax_preRS <- all_otu_w_tax_id[which(rownames(all_otu_w_tax_id) %in% preRS_samps),]
#Subset to only top 20 taxa 
preRS_means <- colMeans(otu_w_tax_preRS,na.rm=T)
top20_otu_w_tax_preRS <- otu_w_tax_preRS[,order(preRS_means)]
top20_otu_w_tax_preRS <- subset(top20_otu_w_tax_preRS,select = -c(1:(ncol(top20_otu_w_tax_preRS)-20)))
#instances of [ or - mess up the function so have to remove these
colnames(top20_otu_w_tax_preRS) <- gsub("\\[|\\]",replacement="",colnames(top20_otu_w_tax_preRS))
colnames(top20_otu_w_tax_preRS) <- gsub("-",replacement="",colnames(top20_otu_w_tax_preRS))
top20_tax_preRS <- colnames(top20_otu_w_tax_preRS)
top20_otu_w_tax_rep_preRS <- as.data.frame(apply(top20_otu_w_tax_preRS,2,rep,each=4))
metadata_w_top20tax_preRS <- cbind(head(metadata,-1),top20_otu_w_tax_rep_preRS)
write.csv(metadata_w_top20tax_preRS,"/Users/Riley/Box Sync/Riley's Documents/Keim Lab/Resistant Starch/RS Data/Microbiome/DADA2_350_unfiltered/metadata_w_top20tax_preRS.csv")

#Health outcomes: glucose and insulin (fasting and iAUC), HOMA-IR, stool pH, hydrogen and methane (AUC), SCFA, BA
out <- c("GluciAUC","InsiAUC","Gluc0","Ins0","HOMAIR","stoolPH","CH4AUC","H2AUC",
         "total_SCFA", "acetate", "propionate", "butyrate", "acetate_per", "propionate_per", "butyrate_per",
         "total_BA", "total_primary_BA", "total_secondary_BA", "total_conj_BA", "CA", "CDCA", "TCDCA", "GCA",
          "GCDCA", "DCA", "LCA", "TDCA", "GDCA", "UDCA")
glyc <- c("GluciAUC","InsiAUC","Gluc0","Ins0","HOMAIR")
ferm <- c("stoolPH","CH4AUC","H2AUC")
scfa <- c("total_SCFA", "acetate", "propionate", "butyrate", "acetate_per", "propionate_per", "butyrate_per")
bile <- c("total_BA", "total_primary_BA", "total_secondary_BA", "total_conj_BA", "CA", "CDCA", "TCDCA", "GCA",
          "GCDCA", "DCA", "LCA", "TDCA", "GDCA", "UDCA")
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}
outlierreplacement <- function(dataframe){
  dataframe %>%          
    map_if(is.numeric, ~ replace(.x, .x %in% boxplot.stats(.x)$out, NA)) %>%
    bind_cols 
}
correltest <- function(data,vars,taxon,method){
  dat <- data[,c(vars,taxon)]
  #dat <- data.matrix(dat)
  #browser()
  #dat <- outlierreplacement(dat) #replace outliers with NA
  corrtable <- rcorr(as.matrix(dat),type=method)
  rtable <- round(corrtable$r,2)
  ptable <- round(corrtable$P,2)
  corrtable <- flattenCorrMatrix(corrtable$r,corrtable$P)
  corrtable <- corrtable %>% mutate_if(is.numeric, ~round(., 2))
  output <- list("corrtable"=corrtable,"rtable"=rtable,"ptable"=ptable)
  return(output)
}
corrscat <- function(data,vars,taxon,subset=NULL,method){
  dat <- correltest(data,vars,taxon,method)
  corr <- dat$corrtable
  sig <- corr[corr$p<=0.05,]
  sig <- sig[!is.na(sig$row),] #remove rows with NA
  sig <- sig[(sig$row %in% vars) & (sig$column %in% taxon),] #get rows with outcomes and taxa
  if(!is.null(subset)){
    sig <- sig[(sig$row %in% subset),]
  }
  sig <- sig %>% mutate_if(is.factor, as.character)
  sig <- sig[order(sig$row),]
  x_labels <- sub("_[^_]+$", "", sig[,"column"])
  data <- outlierreplacement(data)
  plots <- list()
  for(i in 1:nrow(sig)){
    yvar <- sig[i,1]
    xvar <- sig[i,2]
    plots[[i]] <- ggscatter(data, x = xvar, y = yvar, 
                            add = "reg.line", conf.int = TRUE, 
                            cor.coef = TRUE, cor.method = method,
                            xlab = x_labels[i], ylab = yvar)
  }
  return(plots)
}

alltax_out <- correltest(metadata_w_alltax,out,all_tax,method="pearson")
alltax_out_corr <- alltax_out$corrtable
alltax_out_sig <- alltax_out_corr[alltax_out_corr$p<=0.05,]
alltax_out_sig <- alltax_out_sig[!is.na(alltax_out_sig$row),] #remove rows with NA
alltax_out_sig <- alltax_out_sig[(alltax_out_sig$row %in% out) & (alltax_out_sig$column %in% all_tax),] #get rows with outcomes and taxa
#Plots for glycemic response
alltax_glyc_plots <- corrscat(metadata_w_alltax,out,all_tax,subset=glyc,method="pearson")
grid.arrange(grobs=alltax_glyc_plots,ncol=4) #1000x700
#Plots for fermentation response
alltax_ferm_plots <- corrscat(metadata_w_alltax,out,all_tax,subset=ferm,method="pearson")
grid.arrange(grobs=alltax_ferm_plots,ncol=3) #1015x600
#Plots for SCFA
alltax_scfa_plots <- corrscat(metadata_w_alltax,out,all_tax,subset=scfa,method="pearson")
grid.arrange(grobs=alltax_scfa_plots,ncol=5)
#Plots for BA
alltax_bile_plots <- corrscat(metadata_w_alltax,out,all_tax,subset=bile,method="pearson")
grid.arrange(grobs=alltax_bile_plots,ncol=6)
#Separate into categories: primary, secondary, conjugated, and totals
p_BA <- c("CA","CDCA")
s_BA <- c("DCA","LCA")
c_BA <- c("TCDCA", "GCA", "GCDCA", "TDCA", "GDCA", "UDCA")
tot_BA <- c("total_BA","total_conj_BA","total_primary_BA","total_secondary_BA")
alltax_pBA_plots <- corrscat(metadata_w_alltax,out,all_tax,subset=p_BA,method="pearson")
grid.arrange(grobs=alltax_pBA_plots,ncol=3)
alltax_sBA_plots <- corrscat(metadata_w_alltax,out,all_tax,subset=s_BA,method="pearson")
grid.arrange(grobs=alltax_sBA_plots,ncol=2)
alltax_cBA_plots <- corrscat(metadata_w_alltax,out,all_tax,subset=c_BA,method="pearson")
grid.arrange(grobs=alltax_cBA_plots,ncol=6)
alltax_totBA_plots <- corrscat(metadata_w_alltax,out,all_tax,subset=tot_BA,method="pearson")
grid.arrange(grobs=alltax_totBA_plots,ncol=6)

#repeat for all baseline taxa and difference in outcomes
for(i in 1:length(out)){
  diff_var <- paste("diff",out[i],sep="_")
  metadata_C <-metadata_w_alltax_base[metadata_w_alltax_base$Trt=="C",] 
  metadata_RS <-metadata_w_alltax_base[metadata_w_alltax_base$Trt=="RS",] 
  diff_val <- rep(metadata_RS[,out[i]] - metadata_C[,out[i]],each=4)
  metadata_w_alltax_base[,diff_var] <- diff_val
}
diff_var <- paste("diff",out,sep="_")
metadata_RS <- metadata_w_alltax_base[(metadata_w_alltax_base$Trt3=="RS"),] #Subset to only RS samples
#correlation
allbasetax_diffout <- correltest(metadata_RS,diff_var,all_tax_base,method="pearson")
allbasetax_diffout_corr <- allbasetax_diffout$corrtable
allbasetax_diffout_sig <- allbasetax_diffout_corr[allbasetax_diffout_corr$p<=0.05,]
allbasetax_diffout_sig <- allbasetax_diffout_sig[!is.na(allbasetax_diffout_sig$row),] #remove rows with NA
allbasetax_diffout_sig <- allbasetax_diffout_sig[(allbasetax_diffout_sig$row %in% diff_var) & (allbasetax_diffout_sig$column %in% all_tax_base),] #get rows with outcomes and taxa
allbasetax_diffout_plots <- corrscat(metadata_RS,diff_var,all_tax_base,method="pearson")
grid.arrange(grobs=allbasetax_diffout_plots,ncol=8)
#Get subsets of plots
diff_glyc <- unique (grep(paste(glyc,collapse="|"), 
                          diff_var, value=TRUE))
diff_ferm <- unique (grep(paste(ferm,collapse="|"), 
                          diff_var, value=TRUE))
diff_scfa <- unique (grep(paste(scfa,collapse="|"), 
                          diff_var, value=TRUE))
diff_ba <- unique (grep(paste(bile,collapse="|"), 
                          diff_var, value=TRUE))
#Glycemic outcomes
allbasetax_diffglyc_plots <- corrscat(metadata_RS,diff_var,all_tax_base,subset=diff_glyc,method="pearson")
grid.arrange(grobs=allbasetax_diffglyc_plots,ncol=2)
#Fermentation outcomes
allbasetax_diffferm_plots <- corrscat(metadata_RS,diff_var,all_tax_base,subset=diff_ferm,method="pearson")
grid.arrange(grobs=allbasetax_diffferm_plots,ncol=2)
#SCFA outcomes
allbasetax_diffscfa_plots <- corrscat(metadata_RS,diff_var,all_tax_base,subset=diff_scfa,method="pearson")
grid.arrange(grobs=allbasetax_diffscfa_plots,ncol=2)
#Bile acid outcomes
allbasetax_diffba_plots <- corrscat(metadata_RS,diff_var,all_tax_base,subset=diff_ba,method="pearson")
grid.arrange(grobs=allbasetax_diffba_plots,ncol=6)
#Subset by category
diff_pBA <- paste("diff",p_BA,sep="_")
diff_sBA <- paste("diff",s_BA,sep="_")
diff_cBA <- paste("diff",c_BA,sep="_")
diff_totBA <- paste("diff",tot_BA,sep="_")
allbasetax_diffpBA_plots <- corrscat(metadata_RS,diff_var,all_tax_base,subset=diff_pBA,method="pearson")
grid.arrange(grobs=allbasetax_diffpBA_plots,ncol=2)
allbasetax_diffsBA_plots <- corrscat(metadata_RS,diff_var,all_tax_base,subset=diff_sBA,method="pearson")
grid.arrange(grobs=allbasetax_diffsBA_plots,ncol=3)
allbasetax_diffcBA_plots <- corrscat(metadata_RS,diff_var,all_tax_base,subset=diff_cBA,method="pearson")
grid.arrange(grobs=allbasetax_diffcBA_plots,ncol=4)
allbasetax_difftotBA_plots <- corrscat(metadata_RS,diff_var,all_tax_base,subset=diff_totBA,method="pearson")
grid.arrange(grobs=allbasetax_difftotBA_plots,ncol=3)

#Repeat with top 20 and top 20 baseline correlation with fermentation and SCFAs
metadata_w_top20tax <- read.csv("/Users/Riley/Box Sync/Riley's Documents/Keim Lab/Resistant Starch/RS Data/Microbiome/DADA2_350_unfiltered/metadata_w_top20tax.csv",sep = ",",row.names = 1,stringsAsFactors = FALSE)
metadata_w_top20tax_base <- read.csv("/Users/Riley/Box Sync/Riley's Documents/Keim Lab/Resistant Starch/RS Data/Microbiome/DADA2_350_unfiltered/metadata_w_top20tax_base.csv",sep = ",",row.names = 1,stringsAsFactors = FALSE)
top20tax_fermscfa <- correltest(metadata_w_top20tax,c(ferm,scfa),top20_tax,method="pearson")
top20tax_fermscfa_corr <- top20tax_fermscfa$corrtable
top20tax_fermscfa_corr <- top20tax_fermscfa_corr[!is.na(top20tax_fermscfa_corr$row),] #remove rows with NA
top20tax_fermscfa_corr <- top20tax_fermscfa_corr[(top20tax_fermscfa_corr$row %in% c(ferm,scfa)) & (top20tax_fermscfa_corr$column %in% top20_tax),] #get rows with outcomes and taxa
write.csv(top20tax_fermscfa_corr,"/Users/Riley/Box Sync/Riley's Documents/Keim Lab/Resistant Starch/RS Data/Microbiome/DADA2_350_unfiltered/top20tax_fermscfa_corr.csv")
top20tax_fermscfa_sig <- top20tax_fermscfa_corr[top20tax_fermscfa_corr$p<=0.05,]

top20tax_base_fermscfa <- correltest(metadata_w_top20tax_base,c(ferm,scfa),top20_tax_base,method="pearson")
top20tax_base_fermscfa_corr <- top20tax_base_fermscfa$corrtable
top20tax_base_fermscfa_corr <- top20tax_base_fermscfa_corr[!is.na(top20tax_base_fermscfa_corr$row),] #remove rows with NA
top20tax_base_fermscfa_corr <- top20tax_base_fermscfa_corr[(top20tax_base_fermscfa_corr$row %in% c(ferm,scfa)) & (top20tax_base_fermscfa_corr$column %in% top20_tax_base),] #get rows with outcomes and taxa
write.csv(top20tax_base_fermscfa_corr,"/Users/Riley/Box Sync/Riley's Documents/Keim Lab/Resistant Starch/RS Data/Microbiome/DADA2_350_unfiltered/top20tax_base_fermscfa_corr.csv")
top20tax_base_fermscfa_sig <- top20tax_base_fermscfa_corr[top20tax_base_fermscfa_corr$p<=0.05,]

#Repeat with top 20 baseline and top 20 pre-RS correlation with fermentation and SCFAs after RS treatment only
metadata_w_top20tax_base_RS <- subset(metadata_w_top20tax_base,metadata_w_top20tax_base$Trt3=="RS")
top20tax_base_fermscfa_RS <- correltest(metadata_w_top20tax_base_RS,c(ferm,scfa),top20_tax_base,method="pearson")
top20tax_base_fermscfa_RS_corr <- top20tax_base_fermscfa_RS$corrtable
top20tax_base_fermscfa_RS_corr <- top20tax_base_fermscfa_RS_corr[!is.na(top20tax_base_fermscfa_RS_corr$row),] #remove rows with NA
top20tax_base_fermscfa_RS_corr <- top20tax_base_fermscfa_RS_corr[(top20tax_base_fermscfa_RS_corr$row %in% c(ferm,scfa)) & (top20tax_base_fermscfa_RS_corr$column %in% top20_tax_base),] #get rows with outcomes and taxa
write.csv(top20tax_base_fermscfa_RS_corr,"/Users/Riley/Box Sync/Riley's Documents/Keim Lab/Resistant Starch/RS Data/Microbiome/DADA2_350_unfiltered/top20tax_base_fermscfa_RS_corr.csv")
top20tax_base_fermscfa_RS_sig <- top20tax_base_fermscfa_RS_corr[top20tax_base_fermscfa_RS_corr$p<=0.05,]

metadata_w_top20tax_preRS_RS <- subset(metadata_w_top20tax_preRS,metadata_w_top20tax_preRS$Trt3=="RS")
top20tax_preRS_fermscfa_RS <- correltest(metadata_w_top20tax_preRS_RS,c(ferm,scfa),top20_tax_preRS,method="pearson")
top20tax_preRS_fermscfa_RS_corr <- top20tax_preRS_fermscfa_RS$corrtable
top20tax_preRS_fermscfa_RS_corr <- top20tax_preRS_fermscfa_RS_corr[!is.na(top20tax_preRS_fermscfa_RS_corr$row),] #remove rows with NA
top20tax_preRS_fermscfa_RS_corr <- top20tax_preRS_fermscfa_RS_corr[(top20tax_preRS_fermscfa_RS_corr$row %in% c(ferm,scfa)) & (top20tax_preRS_fermscfa_RS_corr$column %in% top20_tax_preRS),] #get rows with outcomes and taxa
write.csv(top20tax_preRS_fermscfa_RS_corr,"/Users/Riley/Box Sync/Riley's Documents/Keim Lab/Resistant Starch/RS Data/Microbiome/DADA2_350_unfiltered/top20tax_preRS_fermscfa_RS_corr.csv")
top20tax_preRS_fermscfa_RS_sig <- top20tax_preRS_fermscfa_RS_corr[top20tax_preRS_fermscfa_RS_corr$p<=0.05,]

#Create heatmaps for top 20 taxa correlations
top20tax_fermscfa_corr <- read.csv("/Users/Riley/Box Sync/Riley's Documents/Keim Lab/Resistant Starch/RS Data/Microbiome/DADA2_350_unfiltered/top20tax_fermscfa_corr.csv",sep = ",",row.names = 1,stringsAsFactors = FALSE)
top20_matrix <- matrix(data = top20tax_fermscfa_corr$cor, nrow = length(unique(top20tax_fermscfa_corr$row)), ncol = length(unique(top20tax_fermscfa_corr$column)))
rownames(top20_matrix) <- unique(top20tax_fermscfa_corr$row)
colnames(top20_matrix) <- unique(top20tax_fermscfa_corr$column)
colnames(top20_matrix) <- gsub("(_[^_]+)_.*", "\\1", colnames(top20_matrix))

top20tax_base_fermscfa_corr <- read.csv("/Users/Riley/Box Sync/Riley's Documents/Keim Lab/Resistant Starch/RS Data/Microbiome/DADA2_350_unfiltered/top20tax_base_fermscfa_corr.csv",sep = ",",row.names = 1,stringsAsFactors = FALSE)
top20_base_matrix <- matrix(data = top20tax_base_fermscfa_corr$cor, nrow = length(unique(top20tax_base_fermscfa_corr$row)), ncol = length(unique(top20tax_base_fermscfa_corr$column)))
rownames(top20_base_matrix) <- unique(top20tax_base_fermscfa_corr$row)
colnames(top20_base_matrix) <- unique(top20tax_base_fermscfa_corr$column)
colnames(top20_base_matrix) <- gsub("(_[^_]+)_.*", "\\1", colnames(top20_base_matrix))

top20tax_base_RS_fermscfa_corr <- read.csv("/Users/Riley/Box Sync/Riley's Documents/Keim Lab/Resistant Starch/RS Data/Microbiome/DADA2_350_unfiltered/top20tax_base_fermscfa_RS_corr.csv",sep = ",",row.names = 1,stringsAsFactors = FALSE)
top20_base_RS_matrix <- matrix(data = top20tax_base_RS_fermscfa_corr$cor, nrow = length(unique(top20tax_base_RS_fermscfa_corr$row)), ncol = length(unique(top20tax_base_RS_fermscfa_corr$column)))
rownames(top20_base_RS_matrix) <- unique(top20tax_base_RS_fermscfa_corr$row)
colnames(top20_base_RS_matrix) <- unique(top20tax_base_RS_fermscfa_corr$column)
colnames(top20_base_RS_matrix) <- gsub("(_[^_]+)_.*", "\\1", colnames(top20_base_RS_matrix))

top20tax_preRS_RS_fermscfa_corr <- read.csv("/Users/Riley/Box Sync/Riley's Documents/Keim Lab/Resistant Starch/RS Data/Microbiome/DADA2_350_unfiltered/top20tax_preRS_fermscfa_RS_corr.csv",sep = ",",row.names = 1,stringsAsFactors = FALSE)
top20_preRS_RS_matrix <- matrix(data = top20tax_preRS_RS_fermscfa_corr$cor, nrow = length(unique(top20tax_preRS_RS_fermscfa_corr$row)), ncol = length(unique(top20tax_preRS_RS_fermscfa_corr$column)))
rownames(top20_preRS_RS_matrix) <- unique(top20tax_preRS_RS_fermscfa_corr$row)
colnames(top20_preRS_RS_matrix) <- unique(top20tax_preRS_RS_fermscfa_corr$column)
colnames(top20_preRS_RS_matrix) <- gsub("(_[^_]+)_.*", "\\1", colnames(top20_preRS_RS_matrix))

library("RColorBrewer")
library("pheatmap")
col <- colorRampPalette(brewer.pal(9, "RdBu"))(256)
newnames_20 <- lapply(colnames(top20_matrix),function(x) bquote(italic(.(x)))) #make taxa names italic
top20_heatmap <- pheatmap(t(top20_matrix), cluster_rows = FALSE, cluster_cols = FALSE, 
                           col = col, scale="none", margins=c(5,10),angle_col = "45",
                           labels_row = as.expression(newnames_20))
newnames_20base <- lapply(colnames(top20_base_matrix),function(x) bquote(italic(.(x))))
top20_base_heatmap <- pheatmap(t(top20_base_matrix), cluster_rows = FALSE, cluster_cols = FALSE, 
                          col = col, scale="none", margins=c(5,10),angle_col = "45",
                          labels_row = as.expression(newnames_20base))
newnames_20baseRS <- lapply(colnames(top20_base_RS_matrix),function(x) bquote(italic(.(x))))
top20_base_RS_heatmap <- pheatmap(t(top20_base_RS_matrix), cluster_rows = FALSE, cluster_cols = FALSE, 
                               col = col, scale="none", margins=c(5,10),angle_col = "45",
                               labels_row = as.expression(newnames_20baseRS))
newnames_20preRSRS <- lapply(colnames(top20_preRS_RS_matrix),function(x) bquote(italic(.(x))))
top20_preRS_RS_heatmap <- pheatmap(t(top20_preRS_RS_matrix), cluster_rows = FALSE, cluster_cols = FALSE, 
                                  col = col, scale="none", margins=c(5,10),angle_col = "45",
                                  labels_row = as.expression(newnames_20preRSRS))
