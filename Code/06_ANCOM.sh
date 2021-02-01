#Differential abundance testing (ANCOM)
qiime composition add-pseudocount \
  --i-table ~/FL103/DADA2_350_unfiltered/table_filtered.qza \
  --o-composition-table ~/FL103/DADA2_350_unfiltered/table_filt_pseudocount.qza

qiime composition ancom \
  --i-table ~/FL103/DADA2_350_unfiltered/table_filt_pseudocount.qza \
  --m-metadata-file ~/FL103/FL103_MapFile.txt \
  --m-metadata-column Trt \
  --o-visualization ~/FL103/DADA2_350_unfiltered/diversity_DADA2_nomock/ANCOM_table_treatment.qzv

qiime composition ancom \
  --i-table ~/FL103/DADA2_350_unfiltered/table_filt_pseudocount.qza \
  --m-metadata-file ~/FL103/FL103_MapFile.txt \
  --m-metadata-column Trt2 \
  --o-visualization ~/FL103/DADA2_350_unfiltered/diversity_DADA2_nomock/ANCOM_table_treatment2.qzv

qiime composition ancom \
  --i-table ~/FL103/DADA2_350_unfiltered/table_filt_pseudocount.qza \
  --m-metadata-file ~/FL103/FL103_MapFile.txt \
  --m-metadata-column Trt3 \
  --o-visualization ~/FL103/DADA2_350_unfiltered/diversity_DADA2_nomock/ANCOM_table_treatment3.qzv

#Collapse at genus level and re-run
qiime taxa collapse \
  --i-table ~/FL103/DADA2_350_unfiltered/table_filtered.qza \
  --i-taxonomy ~/FL103/DADA2_350_unfiltered/taxonomy.qza \
  --p-level 6 \
  --o-collapsed-table ~/FL103/DADA2_350_unfiltered/table_filtered-L6.qza

qiime composition add-pseudocount \
  --i-table ~/FL103/DADA2_350_unfiltered/table_filtered-L6.qza \
  --o-composition-table ~/FL103/DADA2_350_unfiltered/table_filtered-L6_pseudocount.qza

qiime composition ancom \
  --i-table ~/FL103/DADA2_350_unfiltered/table_filtered-L6_pseudocount.qza \
  --m-metadata-file ~/FL103/FL103_MapFile.txt \
  --m-metadata-column Trt \
  --o-visualization ~/FL103/DADA2_350_unfiltered/diversity_DADA2_nomock/ANCOM_L6-table_treatment.qzv

qiime composition ancom \
  --i-table ~/FL103/DADA2_350_unfiltered/table_filtered-L6_pseudocount.qza \
  --m-metadata-file ~/FL103/FL103_MapFile.txt \
  --m-metadata-column Trt2 \
  --o-visualization ~/FL103/DADA2_350_unfiltered/diversity_DADA2_nomock/ANCOM_L6-table_treatment2.qzv

qiime composition ancom \
  --i-table ~/FL103/DADA2_350_unfiltered/table_filtered-L6_pseudocount.qza \
  --m-metadata-file ~/FL103/FL103_MapFile.txt \
  --m-metadata-column Trt3 \
  --o-visualization ~/FL103/DADA2_350_unfiltered/diversity_DADA2_nomock/ANCOM_L6-table_treatment3.qzv

#Collapse at family level and re-run
qiime taxa collapse \
  --i-table ~/FL103/DADA2_350_unfiltered/table_filtered.qza \
  --i-taxonomy ~/FL103/DADA2_350_unfiltered/taxonomy.qza \
  --p-level 5 \
  --o-collapsed-table ~/FL103/DADA2_350_unfiltered/table_filtered-L5.qza

qiime composition add-pseudocount \
  --i-table ~/FL103/DADA2_350_unfiltered/table_filtered-L5.qza \
  --o-composition-table ~/FL103/DADA2_350_unfiltered/table_filtered-L5_pseudocount.qza

qiime composition ancom \
  --i-table ~/FL103/DADA2_350_unfiltered/table_filtered-L5_pseudocount.qza \
  --m-metadata-file ~/FL103/FL103_MapFile.txt \
  --m-metadata-column Trt3 \
  --o-visualization ~/FL103/DADA2_350_unfiltered/diversity_DADA2_nomock/ANCOM_L5-table_treatment3.qzv


