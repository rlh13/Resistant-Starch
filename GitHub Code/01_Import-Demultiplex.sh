module load bio3
module load qiime/2-2020.2
source activate qiime2-2020.2

#Import data (sequences + barcodes)
#Files in FL103 folder are “sequences.fastq.gz” and “barcodes.fastq.gz”
qiime tools import \
--type EMPSingleEndSequences \
--input-path ~/FL103/barcodes_unfiltered/ \
--output-path ~/FL103/FL103seqs_unfiltered.qza

#demultiplexing 
qiime demux emp-single \
--i-seqs ~/FL103/FL103seqs_unfiltered.qza \
--m-barcodes-file ~/FL103/FL103_MapFile.txt \
--m-barcodes-column BarcodeSequence \
--p-no-golay-error-correction \
--o-per-sample-sequences ~/FL103/demux_unfiltered.qza \
--o-error-correction-details ~/FL103/demux-details_unfiltered.qza

#summarize demux results
qiime demux summarize \
--i-data ~/FL103/demux_unfiltered.qza \
--o-visualization ~/FL103/demux_unfiltered.qzv

