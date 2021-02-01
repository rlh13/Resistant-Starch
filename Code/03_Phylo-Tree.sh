#building phylogenetic tree 
#Conduct alignment for tree building 
qiime alignment mafft \
--i-sequences ~/FL103/DADA2_350_unfiltered/representative_sequences.qza \
--o-alignment ~/FL103/DADA2_350_unfiltered/representative_sequences_aligned.qza

#Mask (i.e., filter) unconserved and highly gapped columns from an alignment.
qiime alignment mask \
--i-alignment ~/FL103/DADA2_350_unfiltered/representative_sequences_aligned.qza \
--o-masked-alignment ~/FL103/DADA2_350_unfiltered/representative_sequences_aligned_masked.qza

#Generate a tree for phylogenetic diversity analyses
qiime phylogeny fasttree \
--i-alignment ~/FL103/DADA2_350_unfiltered/representative_sequences_aligned_masked.qza \
--o-tree ~/FL103/DADA2_350_unfiltered/unrooted-tree.qza

#place the root at the midpoint of the longest tip-to-tip distance on the tree
qiime phylogeny midpoint-root \
--i-tree ~/FL103/DADA2_350_unfiltered/unrooted-tree.qza \
--o-rooted-tree ~/FL103/DADA2_350_unfiltered/rooted-tree.qza

