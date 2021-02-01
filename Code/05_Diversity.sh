#core diversity analyses (rarefaction, beta diversity, alpha diversity)
#sampling depth of 503 chosen from table.qzv interactive sample detail. Excludes 1 sample #(FL103.058.2; sample depth = 57)
qiime diversity core-metrics-phylogenetic \
  --i-phylogeny ~/FL103/DADA2_350_unfiltered/rooted-tree.qza \
  --i-table ~/FL103/DADA2_350_unfiltered/table_filtered.qza \
  --p-sampling-depth 533 \
  --m-metadata-file ~/FL103/FL103_MapFile.txt \
  --output-dir ~/FL103/DADA2_350_unfiltered/diversity_DADA2_nomock/

#Plot alpha div
#used previously created data with table_filtered.qza (excludes mock community sample)
qiime diversity alpha-rarefaction \
  --i-table ~/FL103/DADA2_350_unfiltered/diversity_DADA2_nomock/rarefied_table.qza \
  --i-phylogeny ~/FL103/DADA2_350_unfiltered/rooted-tree.qza \
  --m-metadata-file ~/FL103/FL103_MapFile.txt \
  --p-max-depth 533 \
  --o-visualization ~/FL103/DADA2_350_unfiltered/diversity_DADA2_nomock/alpha-rarefaction.qzv

#UniFrac
qiime diversity beta-group-significance \
--i-distance-matrix ~/FL103/DADA2_350_unfiltered/diversity_DADA2_nomock/unweighted_unifrac_distance_matrix.qza \
--m-metadata-file ~/FL103/FL103_MapFile.txt \
--m-metadata-column Trt3 \
--o-visualization ~/FL103/DADA2_350_unfiltered/diversity_DADA2_nomock/unweighted-unifrac-treatment-significance.qzv \
--p-pairwise

qiime diversity beta-group-significance \
--i-distance-matrix ~/FL103/DADA2_350_unfiltered/diversity_DADA2_nomock/weighted_unifrac_distance_matrix.qza \
--m-metadata-file ~/FL103/FL103_MapFile.txt \
--m-metadata-column Trt \
--o-visualization ~/FL103/DADA2_350_unfiltered/diversity_DADA2_nomock/weighted-unifrac-treatment-significance.qzv \
--p-pairwise

qiime diversity beta-group-significance \
--i-distance-matrix ~/FL103/DADA2_350_unfiltered/diversity_DADA2_nomock/weighted_unifrac_distance_matrix.qza \
--m-metadata-file ~/FL103/FL103_MapFile.txt \
--m-metadata-column Trt2 \
--o-visualization ~/FL103/DADA2_350_unfiltered/diversity_DADA2_nomock/weighted-unifrac-treatment2-significance.qzv \
--p-pairwise

qiime diversity beta-group-significance \
--i-distance-matrix ~/FL103/DADA2_350_unfiltered/diversity_DADA2_nomock/weighted_unifrac_distance_matrix.qza \
--m-metadata-file ~/FL103/FL103_MapFile.txt \
--m-metadata-column Trt3 \
--o-visualization ~/FL103/DADA2_350_unfiltered/diversity_DADA2_nomock/weighted-unifrac-treatment3-significance.qzv \
--p-pairwise

