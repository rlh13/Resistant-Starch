#Activate QIIME1 conda
conda create -n qiime1 python=2.7 qiime matplotlib=1.4.3 mock nose -c bioconda
source activate qiime1

#separate reads from barcodes 
extract_barcodes.py --fastq1 ~/FL103/AllReads_unfiltered.fastq --bc1_len 8 -m ~/FL103/FL103_MapFile.txt -o ~/FL103/barcodes_unfiltered

conda deactivate

#Zip barcode file
gzip ~/FL103/barcodes_unfiltered/barcodes.fastq -r ~/FL103/barcodes_unfiltered/barcodes.fastq.gz 

#Rename “reads” to “sequences” for QIIME2 formatting
mv ~/FL103/barcodes_unfiltered/reads.fastq ~/FL103/barcodes_unfiltered/sequences.fastq

#Zip sequences file
gzip ~/FL103/barcodes_unfiltered/sequences.fastq -r ~/FL103/barcodes_unfiltered/sequences.fastq.gz