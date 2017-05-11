#!/bin/bash 
#
annot=/home/gmsmiao/nextgen1/fibrosis/annotation/Homo_sapiens/Ensembl/CRCh38.78.primary/genes.gtf
for sampleName in ./tophat2_out/tophat2_sorted/*.bam ; do
	name=`echo $sampleName | cut -d. -f2 | cut -d/ -f4`  
	qsub -N $name -pe smp 1 "htseq-count --order=name --format=bam --stranded=no --type=exon --idattr=gene_id $sampleName $annot > ./count_out/$name.count" 
done 
