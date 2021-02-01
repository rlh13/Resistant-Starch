#Load metadata
metadata <- read.csv("/Users/Riley/Box Sync/Riley's Documents/Keim Lab/Resistant Starch/RS Data/Microbiome/FL103_Mapfile.csv",sep=",",header=T,stringsAsFactors=F)
rownames(metadata) <- metadata$X.SampleID
metadata$Trt <- factor(metadata$Trt)
metadata$TrtSeq <- factor(metadata$TrtSeq)
#Load QIIME2 data
library(devtools)
library(qiime2R)
library(phyloseq)
#DADA2_unfilt
phy = qza_to_phyloseq("/Users/Riley/Box Sync/Riley's Documents/Keim Lab/Resistant Starch/RS Data/Microbiome/DADA2_350_unfiltered/table.qza", 
                      "/Users/Riley/Box Sync/Riley's Documents/Keim Lab/Resistant Starch/RS Data/Microbiome/DADA2_350_unfiltered/unrooted-tree.qza",
                      "/Users/Riley/Box Sync/Riley's Documents/Keim Lab/Resistant Starch/RS Data/Microbiome/DADA2_350_unfiltered/trained_tax/taxonomy.qza")
phy <- merge_phyloseq(phy, sample_data(metadata))
phy = prune_taxa(taxa_sums(phy) > 1, phy)
phy_rel <- transform_sample_counts(phy, function(x) x / sum(x) )
saveRDS(phy,"/Users/Riley/Box Sync/Riley's Documents/Keim Lab/Resistant Starch/RS Data/Microbiome/DADA2_350_unfiltered/phy.RDS")
saveRDS(phy_rel,"/Users/Riley/Box Sync/Riley's Documents/Keim Lab/Resistant Starch/RS Data/Microbiome/DADA2_350_unfiltered/phy_rel.RDS")
tax=as.data.frame(phy@tax_table@.Data)
otu <- as.data.frame(otu_table(phy))
otu_w_tax <- merge(otu,tax,by="row.names")
otu_w_tax <- otu_w_tax[ , -which(names(otu_w_tax) %in% "Row.names")]
write.csv(otu_w_tax,"/Users/Riley/Box Sync/Riley's Documents/Keim Lab/Resistant Starch/RS Data/Microbiome/DADA2_350_unfiltered/otu_w_noglomtax.csv")
tax_rel=as.data.frame(phy_rel@tax_table@.Data)
otu_rel <- as.data.frame(otu_table(phy_rel))
otu_w_tax_rel <- merge(otu_rel,tax_rel,by="row.names")
otu_w_tax_rel <- otu_w_tax_rel[ , -which(names(otu_w_tax_rel) %in% "Row.names")]
write.csv(otu_w_tax_rel,"/Users/Riley/Box Sync/Riley's Documents/Keim Lab/Resistant Starch/RS Data/Microbiome/DADA2_350_unfiltered/otu_w_noglomtax_rel.csv")

phy_species <- tax_glom(phy,"Species",NArm = FALSE)
phy_speciesrel <- transform_sample_counts(phy_species, function(x) x / sum(x) )
tax_species=as.data.frame(phy_speciesrel@tax_table@.Data)
otu_species <- as.data.frame(otu_table(phy_speciesrel))
otu_w_speciestax <- merge(otu_species,tax_species,by="row.names")
write.csv(otu_w_speciestax,"/Users/Riley/Box Sync/Riley's Documents/Keim Lab/Resistant Starch/RS Data/Microbiome/DADA2_350_unfiltered/otu_w_speciestax_rel.csv")

phy_genus <- tax_glom(phy,"Genus",NArm = FALSE)
phy_genusrel <- transform_sample_counts(phy_genus, function(x) x / sum(x) )
saveRDS(phy_genus,"/Users/Riley/Box Sync/Riley's Documents/Keim Lab/Resistant Starch/RS Data/Microbiome/DADA2_350_unfiltered/phy_genus.RDS")
saveRDS(phy_genusrel,"/Users/Riley/Box Sync/Riley's Documents/Keim Lab/Resistant Starch/RS Data/Microbiome/DADA2_350_unfiltered/phy_genusrel.RDS")
phy_fam <- tax_glom(phy,"Family",NArm = FALSE)
phy_famrel <- transform_sample_counts(phy_fam, function(x) x / sum(x) )
saveRDS(phy_fam,"/Users/Riley/Box Sync/Riley's Documents/Keim Lab/Resistant Starch/RS Data/Microbiome/DADA2_350_unfiltered/phy_fam.RDS")
saveRDS(phy_famrel,"/Users/Riley/Box Sync/Riley's Documents/Keim Lab/Resistant Starch/RS Data/Microbiome/DADA2_350_unfiltered/phy_famrel.RDS")
phy_phy <- tax_glom(phy,"Phylum",NArm = FALSE)
phy_phyrel <- transform_sample_counts(phy_phy, function(x) x / sum(x) )
saveRDS(phy_phy,"/Users/Riley/Box Sync/Riley's Documents/Keim Lab/Resistant Starch/RS Data/Microbiome/DADA2_350_unfiltered/phy_phy.RDS")
saveRDS(phy_phyrel,"/Users/Riley/Box Sync/Riley's Documents/Keim Lab/Resistant Starch/RS Data/Microbiome/DADA2_350_unfiltered/phy_phyrel.RDS")
#For phy_genus
tax_genus=as.data.frame(phy_genus@tax_table@.Data)
otu_genus <- as.data.frame(otu_table(phy_genus))
otu_w_genustax <- merge(otu_genus,tax_genus,by="row.names")
write.csv(otu_w_genustax,"/Users/Riley/Box Sync/Riley's Documents/Keim Lab/Resistant Starch/RS Data/Microbiome/DADA2_350_unfiltered/otu_w_tax_genus_OTUID.csv")
otu_w_genustax <- otu_w_genustax[ , -which(names(otu_w_genustax) %in% "Row.names")]
write.csv(otu_w_genustax,"/Users/Riley/Box Sync/Riley's Documents/Keim Lab/Resistant Starch/RS Data/Microbiome/DADA2_350_unfiltered/otu_w_genustax.csv")
tax_genusrel=as.data.frame(phy_genusrel@tax_table@.Data)
otu_genusrel <- as.data.frame(otu_table(phy_genusrel))
otu_w_tax_genusrel <- merge(otu_genusrel,tax_genusrel,by="row.names")
write.csv(otu_w_tax_genusrel,"/Users/Riley/Box Sync/Riley's Documents/Keim Lab/Resistant Starch/RS Data/Microbiome/DADA2_350_unfiltered/otu_w_tax_genusrel_OTUID.csv")
otu_w_tax_genusrel <- otu_w_tax_genusrel[ , -which(names(otu_w_tax_genusrel) %in% "Row.names")]
write.csv(otu_w_tax_genusrel,"/Users/Riley/Box Sync/Riley's Documents/Keim Lab/Resistant Starch/RS Data/Microbiome/DADA2_350_unfiltered/otu_w_tax_genusrel.csv")
#For phy_fam
tax_fam=as.data.frame(phy_fam@tax_table@.Data)
otu_fam <- as.data.frame(otu_table(phy_fam))
otu_w_famtax <- merge(otu_fam,tax_fam,by="row.names")
otu_w_famtax <- otu_w_famtax[ , -which(names(otu_w_famtax) %in% "Row.names")]
write.csv(otu_w_famtax,"/Users/Riley/Box Sync/Riley's Documents/Keim Lab/Resistant Starch/RS Data/Microbiome/DADA2_350_unfiltered/otu_w_famtax.csv")
tax_famrel=as.data.frame(phy_famrel@tax_table@.Data)
otu_famrel <- as.data.frame(otu_table(phy_famrel))
otu_w_tax_famrel <- merge(otu_famrel,tax_famrel,by="row.names")
otu_w_tax_famrel <- otu_w_tax_famrel[ , -which(names(otu_w_tax_famrel) %in% "Row.names")]
write.csv(otu_w_tax_famrel,"/Users/Riley/Box Sync/Riley's Documents/Keim Lab/Resistant Starch/RS Data/Microbiome/DADA2_350_unfiltered/otu_w_tax_famrel.csv")
#For phy_phy
tax_phy=as.data.frame(phy_phy@tax_table@.Data)
otu_phy <- as.data.frame(otu_table(phy_phy))
otu_w_phytax <- merge(otu_phy,tax_phy,by="row.names")
otu_w_phytax <- otu_w_phytax[ , -which(names(otu_w_phytax) %in% "Row.names")]
write.csv(otu_w_phytax,"/Users/Riley/Box Sync/Riley's Documents/Keim Lab/Resistant Starch/RS Data/Microbiome/DADA2_350_unfiltered/otu_w_phytax.csv")
tax_phyrel=as.data.frame(phy_phyrel@tax_table@.Data)
otu_phyrel <- as.data.frame(otu_table(phy_phyrel))
otu_w_tax_phyrel <- merge(otu_phyrel,tax_phyrel,by="row.names")
otu_w_tax_phyrel <- otu_w_tax_phyrel[ , -which(names(otu_w_tax_phyrel) %in% "Row.names")]
write.csv(otu_w_tax_phyrel,"/Users/Riley/Box Sync/Riley's Documents/Keim Lab/Resistant Starch/RS Data/Microbiome/DADA2_350_unfiltered/otu_w_tax_phyrel.csv")
