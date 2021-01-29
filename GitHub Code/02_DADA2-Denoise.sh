#run DADA2 denoising
qiime dada2 denoise-single \
--i-demultiplexed-seqs ~/FL103/demux_unfiltered.qza \
--p-trim-left 50 \
--p-trunc-len 350 \
--output-dir ~/FL103/DADA2_350_unfiltered/ 

qiime metadata tabulate \
--m-input-file ~/FL103/DADA2_350_unfiltered/denoising_stats.qza \
--o-visualization ~/FL103/DADA2_350_unfiltered/denoising_stats.qzv

qiime feature-table summarize \
  --i-table ~/FL103/DADA2_350_unfiltered/table.qza \
  --o-visualization ~/FL103/DADA2_350_unfiltered/table.qzv \
  --m-sample-metadata-file ~/FL103/FL103_MapFile.txt

qiime feature-table tabulate-seqs \
--i-data ~/FL103/DADA2_350_unfiltered/representative_sequences.qza \
--o-visualization ~/FL103/DADA2_350_unfiltered/representative_sequences.qzv


