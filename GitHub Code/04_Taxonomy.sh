#Classifier need to be re-trained if diff lib prep method used 
#Training classifier
#download gg_13_8_otus.tar.gz file
#extract files
tar xvf ~/FL103/gg_13_8_otus.tar.gz -C ~/FL103/gg_files/

qiime tools import \
--type 'FeatureData[Sequence]' \
--input-path ~/FL103/gg_files/99_otus.fasta \
--output-path ~/FL103/gg_files/99_otus.qza

qiime tools import \
--type 'FeatureData[Taxonomy]' \
--input-format HeaderlessTSVTaxonomyFormat \
--input-path ~/FL103/gg_files/99_otu_taxonomy.txt \
--output-path ~/FL103/gg_files/ref-taxonomy.qza

qiime feature-classifier extract-reads \
  --i-sequences ~/FL103/gg_files/99_otus.qza \
  --p-f-primer ACTCCTACGGGAGGCAGCAG \
  --p-r-primer GGACTACHVGGGTWTCTAAT \
  --p-trunc-len 487 \
  --p-min-length 300 \
  --p-max-length 500 \
  --o-reads ~/FL103/gg_files/ref-seqs.qza

qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads ~/FL103/gg_files/ref-seqs.qza \
  --i-reference-taxonomy ~/FL103/gg_files/ref-taxonomy.qza \
  --o-classifier ~/FL103/gg_files/classifier.qza

#Using trained classifier
qiime feature-classifier classify-sklearn \
--i-classifier ~/FL103/gg_files/classifier.qza \
--i-reads ~/FL103/DADA2_350_unfiltered/representative_sequences.qza \
--o-classification ~/FL103/DADA2_350_unfiltered/trained_tax/taxonomy.qza

#visualize/adding metadata to the taxonomy file
qiime metadata tabulate \
--m-input-file ~/FL103/DADA2_350_unfiltered/taxonomy.qza \
--o-visualization ~/FL103/DADA2_350_unfiltered/taxonomy.qzv

#filter out chloroplasts and mitochondria
qiime taxa filter-table \
--i-table ~/FL103/DADA2_350_unfiltered/table.qza \
--i-taxonomy ~/FL103/DADA2_350_unfiltered/taxonomy.qza \
--p-exclude mitochondria,chloroplast \
--o-filtered-table ~/FL103/DADA2_350_unfiltered/table-no-mitochondria-no-chloroplast.qza

srun --mem=64000 --partition=high --time=07:00:00   \
qiime feature-table filter-samples \
--i-table ~/FL103/DADA2_350_unfiltered/table-no-mitochondria-no-chloroplast.qza \
--m-metadata-file ~/FL103/FL103_MapFile.txt \
--p-exclude-ids \
--p-where "SubID ='MockComm'" \
--o-filtered-table ~/FL103/DADA2_350_unfiltered/table_filtered.qza

qiime feature-table summarize \
  --i-table ~/FL103/DADA2_350_unfiltered/table-no-mitochondria-no-chloroplast.qza \
  --o-visualization ~/FL103/DADA2_350_unfiltered/table-no-mitochondria-no-chloroplast.qzv \
  --m-sample-metadata-file ~/FL103/FL103_MapFile.txt

qiime feature-table summarize \
  --i-table ~/FL103/DADA2_350_unfiltered/table_filtered.qza \
  --o-visualization ~/FL103/DADA2_350_unfiltered/table_filtered.qzv \
  --m-sample-metadata-file ~/FL103/FL103_MapFile.txt

#generate interactive barplot of taxonomy 
qiime taxa barplot \
--i-table ~/FL103/DADA2_350_unfiltered/diversity_DADA2_nomock/rarefied_table.qza \
--i-taxonomy ~/FL103/DADA2_350_unfiltered/taxonomy.qza \
--m-metadata-file ~/FL103/FL103_MapFile.txt \
--o-visualization ~/FL103/DADA2_350_unfiltered/taxa_summary.qzv

