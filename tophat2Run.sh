#!/bin/bash 
#
batch1=./Seq_data/141119_SN674_0296_AC5P9CACXX/Unaligned/Project_C5P9CACXX
batch2=./Seq_data/150108_SN674_0300_AC5RNUACXX/Unaligned/Project_C5RNUACXX
outputRoot=./tophat2_out
refRoot=/home/gmsmiao/nextgen1/fibrosis/annotation/Homo_sapiens/Ensembl/CRCh38.78.primary
transcriptome=/home/gmsmiao/nextgen1/fibrosis/annotation/Homo_sapiens/Ensembl/CRCh38.78.primary/transcriptome

for sample in $batch1/*;do 
	sampleName=`basename $sample`
	b1_read1=`ls $batch1/$sampleName/*_R1_001.fastq.gz`
	b2_read1=`ls $batch2/$sampleName/*_R1_001.fastq.gz`
	b1_read2=`ls $batch1/$sampleName/*_R2_001.fastq.gz`
	b2_read2=`ls $batch2/$sampleName/*_R2_001.fastq.gz`
	outputDir=$outputRoot/$sampleName
	qsub -N $sampleName -pe smp 12 "tophat2 --num-threads 12 --prefilter-multihits --read-mismatches 2 --read-edit-dist 2 --mate-inner-dist 0 --read-realign-edit-dist 0 -o $outputDir --transcriptome-index=$transcriptome/genes $refRoot/genome $b1_read1,$b2_read1 $b1_read2,$b2_read2"
done 
	
