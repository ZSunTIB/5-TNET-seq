#!/bin/sh
##########################################################################################################
fq_R2=Raw_fastq_R2
##########################################################################################################
#Extract all unmapped reads
bowtie /Users/sz/Data_analysis/Genome/Bowtie_genome/MG1655_v2 -q \
$fq_R2\_R12_cutadapt_remD_6Ntrim_remEM.fastq \
-S $fq_R2\_R12_cutadapt_remD_6Ntrim_remEM_multiple.sam

samtools sort $fq_R2\_R12_cutadapt_remD_6Ntrim_remEM_multiple.sam -o \
$fq_R2\_R12_multiple_sort.sam

mkdir Processed_5TNET_R21_unmapped

#Remove all mapped reads
samtools view -f 4 $fq_R2\_R12_multiple_sort.sam -b -o ./Processed_5TNET_R21_unmapped/$fq_R2\_R12_unmapped.bam

#convert .bam file to .fastq file.
bedtools bamtofastq -i ./Processed_5TNET_R21_unmapped/$fq_R2\_R12_unmapped.bam -fq ./Processed_5TNET_R21_unmapped/$fq_R2\_R12_unmapped.fastq

awk 'BEGIN {OFS="\n"} 
     !/^@/ {next} 
     {header = $0; getline seq; getline plus; getline qual} 
     !seen[header]++ {print header, seq, plus, qual}
' ./Processed_5TNET_R21_unmapped/$fq_R2\_R12_unmapped.fastq > ./Processed_5TNET_R21_unmapped/$fq_R2\_R12_unmapped_rm.fastq

mv ./Processed_5TNET_R21_unmapped/$fq_R2\_R12_unmapped_rm.fastq ./Processed_5TNET_R21_unmapped/$fq_R2\_R12_unmapped.fastq

rm $fq_R2\_R12_cutadapt_remD_6Ntrim_remEM_multiple.sam $fq_R2\_R12_multiple_sort.sam
#rm $fq_R2\_R12_cutadapt_remD_6Ntrim_remEM.fastq
