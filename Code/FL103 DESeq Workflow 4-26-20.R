#DESeq2
#BiocManager::install("DESeq2")
library(DESeq2)
# Note 1: in the code for the "Waste Not, Want Not" paper, they used the Wald test with a 
# parametric fit (and resorted to a local fit, only if the parametric fit yielded an error) 
# Note 2: The multiple-inference correction is Benjamini-Hochberg occurs within the DESeq function
make.report <- function(res.df, taxtable.df, otutable.df){
  sig.df = res.df[which(res.df$pvalue < 0.05),]	
  sig.df.withNames = merge(as.data.frame(sig.df), taxtable.df, by.x="row.names", by.y="row.names", all.x=TRUE, all.y=FALSE)
  names(sig.df.withNames)[names(sig.df.withNames)=="Row.names"] <- "OTU"	 
  df <- merge(sig.df.withNames, otutable.df, by.x="OTU", by.y="row.names", all=FALSE, all.x=TRUE) 
  #move otu abundance column to the second position
  col_idx <- grep("otu.percent", names(df))
  vec <- c(1, col_idx, ((2:ncol(df))[-col_idx+1]))
  df <- df[, vec]
  #sort by otu.percent in decreasing order
  #browser()
  sorted.df <- df[order(df$otu.percent, decreasing=TRUE), ]
  return(sorted.df)
}
phy <- readRDS("/Users/Riley/Box Sync/Riley's Documents/Keim Lab/Resistant Starch/RS Data/Microbiome/DADA2_350_unfiltered/phy.RDS")
#Remove mock community sample
phy <- subset_samples(phy, Trt!="MC")
#Glom by genus but keep higher level classification
phy_genus <- tax_glom(phy,"Genus",NArm = FALSE)

#Is there a difference between intervention groups, not accounting for any other variation
#Remove low abundance sample from phy_genus
lowsamps <- colnames(otu_table(phy_genus))[grep("FL.103.058", colnames(otu_table(phy_genus)))]
samps_to_keep <- !c(colnames(otu_table(phy_genus)) %in% lowsamps)
phy_genus <- prune_samples(samps_to_keep, phy_genus)

taxtable.df <- as.data.frame(phy_genus@tax_table@.Data)
otutable.df <- as.data.frame(otu_table(phy_genus))
options(scipen=999)
otu.means <- rowMeans(otutable.df)
total <- sum(otu.means)
otu.percent <- (otu.means/total)*100
otutable.wperc <- cbind(otutable.df, otu.percent)
#Make sure your comparison group is a factor
sample_data(phy_genus)$Trt3 <- as.factor(sample_data(phy_genus)$Trt3)
sample_data(phy_genus)$SubID <- as.factor(sample_data(phy_genus)$SubID)
#Run DESeq analysis on all groups and extract results
dsq = phyloseq_to_deseq2(phy_genus, ~SubID + Trt3)
dsq_obj <- DESeq(dsq)
resultsNames(dsq_obj)
#Individually create results files for each comparison of interest
res_RSvC <- results(dsq_obj, contrast = c("Trt3", "RS", "C"))
res_PrevC <- results(dsq_obj, contrast = c("Trt3", "C", "Pre_C"))
res_PrevRS <- results(dsq_obj, contrast = c("Trt3", "RS", "Pre_RS"))
res_PreRSvPreC <- results(dsq_obj, contrast = c("Trt3", "Pre_RS", "Pre_C"))
table(res_RSvC$pvalue <0.05) #5 significant results
res_RSvC_df = make.report(as.data.frame(res_RSvC), taxtable.df, otutable.wperc)
write.csv(res_RSvC_df, "/Users/Riley/Box Sync/Riley's Documents/Keim Lab/Resistant Starch/RS Data/Microbiome/DADA2_350_unfiltered/DESeq_FL103_RSvC_allreport.csv")

table(res_PrevC$pvalue <0.05) #2 significant results
res_PrevC_df = make.report(as.data.frame(res_PrevC), taxtable.df, otutable.wperc)
write.csv(res_PrevC_df, "/Users/Riley/Box Sync/Riley's Documents/Keim Lab/Resistant Starch/RS Data/Microbiome/DADA2_350_unfiltered/DESeq_FL103_PrevC_allreport.csv")

table(res_PrevRS$pvalue <0.05) #6 significant results
res_PrevRS_df = make.report(as.data.frame(res_PrevRS), taxtable.df, otutable.wperc)
write.csv(res_PrevRS_df, "/Users/Riley/Box Sync/Riley's Documents/Keim Lab/Resistant Starch/RS Data/Microbiome/DADA2_350_unfiltered/DESeq_FL103_PrevRS_allreport.csv")

table(res_PreRSvPreC$pvalue <0.05) #0 significant results
res_PreRSvPreC_df = make.report(as.data.frame(res_PreRSvPreC), taxtable.df, otutable.wperc)
write.csv(res_PreRSvPreC_df, "/Users/Riley/Box Sync/Riley's Documents/Keim Lab/Resistant Starch/RS Data/Microbiome/DADA2_350_unfiltered/DESeq_FL103_PreRSvPreC_allreport.csv")

#Run DESeq analysis on all groups and extract results
#Account for Treatment sequence
sample_data(phy_genus)$TrtSeq <- as.factor(sample_data(phy_genus)$TrtSeq)
# This design models the difference of the difference of Subject, Treatment, 
# and the Sequence-Specific effect of Treatment
# (Sequence main effect excluded because redundant with SubID - only varies between subjects)
dsq_seq = phyloseq_to_deseq2(phy_genus, ~ TrtSeq + Trt3  + TrtSeq:Trt3)
# This performs a likelihood ratio test, which determines if the increased likelihood of the data
# in the full model, where we include all terms, is more than expected if the terms excluded in 
# the reduced model are really zero.  
# In the reduced model, we're removing the interaction term
# Transcripts with small p values from this test are those which 
# showed a response treatment-specific effect
dsq_seq_obj <- DESeq(dsq_seq, test="LRT", reduced = ~  TrtSeq + Trt3 )
resultsNames(dsq_seq_obj)
res_seq <- as.data.frame(results(dsq_seq_obj))
table(res_seq$pvalue <0.05) #none
