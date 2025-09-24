#!/bin/sh
##########################################################################################################
fq_R1=Raw_fastq_R1
fq_R2=Raw_fastq_R2
##########################################################################################################
cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -o $fq_R1\_cutadapt.fastq \
$fq_R1\.fastq

cutadapt -a CTGTAGGCACCATCAATGATCGTCGGAGTGT -o $fq_R2\_cutadapt.fastq \
$fq_R2\.fastq

#Merge R2 and R1 reads with 4-nt overlap
flash -m 4 $fq_R2\_cutadapt.fastq $fq_R1\_cutadapt.fastq

#Rename the merged fastq file
mv out.extendedFrags.fastq $fq_R2\_R12_cutadapt.fastq
rm out.histogram out.hist out.notCombined_1.fastq out.notCombined_2.fastq

#Remove PCR duplicates
clumpify.sh in=$fq_R2\_R12_cutadapt.fastq \
out=$fq_R2\_R12_cutadapt_remD.fastq subs=0 dedupe=t

#Remove the 6-nt barcode.
cutadapt --cut -6 -o $fq_R2\_R12_cutadapt_remD_6Ntrim.fastq \
$fq_R2\_R12_cutadapt_remD.fastq

bioawk -c fastx 'length($seq) > 3 {print "@"$name" "$comment"\n"$seq"\n+\n"$qual}' \
$fq_R2\_R12_cutadapt_remD_6Ntrim.fastq \
> $fq_R2\_R12_cutadapt_remD_6Ntrim_remEM.fastq

bowtie /Users/sz/Data_analysis/Genome/Bowtie_genome/MG1655_v2 -q \
$fq_R2\_R12_cutadapt_remD_6Ntrim_remEM.fastq \
--best --strata -m 1 -S $fq_R2\_R12_cutadapt_remD_6Ntrim_remEM.sam

samtools sort $fq_R2\_R12_cutadapt_remD_6Ntrim_remEM.sam -o $fq_R2\_R12_sort.sam


mkdir Processed_5TNET_R21

#Remove all unmapped and multiple mapped reads
samtools view -F 4 $fq_R2\_R12_sort.sam -b -o ./Processed_5TNET_R21/$fq_R2\_R12_uniq.bam
samtools index ./Processed_5TNET_R21/$fq_R2\_R12_uniq.bam

#Remove all uniquely mapped reads
samtools view -f 4 $fq_R2\_R12_sort.sam -b -o ./Processed_5TNET_R21/$fq_R2\_R12_unmapped_multiple.bam

rm $fq_R2\_R12_cutadapt_remD_6Ntrim_remEM.sam $fq_R2\_R12_sort.sam
##########
#Generate the 5' end reads count for each coordinate
bedtools genomecov -d -5 -strand + -ibam ./Processed_5TNET_R21/$fq_R2\_R12_uniq.bam \
> ./Processed_5TNET_R21/$fq_R2\_R12_uniq_Positive.wig
bedtools genomecov -d -5 -strand - -ibam ./Processed_5TNET_R21/$fq_R2\_R12_uniq.bam \
> ./Processed_5TNET_R21/$fq_R2\_R12_uniq_Negative.wig
##########
#Generate the reads coverage for each coordinate.
bedtools genomecov -d -strand + -ibam ./Processed_5TNET_R21/$fq_R2\_R12_uniq.bam \
> ./Processed_5TNET_R21/$fq_R2\_R12_uniq_Positive_Allread_depth.wig
bedtools genomecov -d -strand - -ibam ./Processed_5TNET_R21/$fq_R2\_R12_uniq.bam \
> ./Processed_5TNET_R21/$fq_R2\_R12_uniq_Negative_Allread_depth.wig
##########
#Generate the read depth (3' END) of each coordiante used for PS selection.
bedtools genomecov -d -3 -strand + -ibam ./Processed_5TNET_R21/$fq_R2\_R12_uniq.bam \
> ./Processed_5TNET_R21/$fq_R2\_R12_uniq_3end_Positive.wig
bedtools genomecov -d -3 -strand - -ibam ./Processed_5TNET_R21/$fq_R2\_R12_uniq.bam \
> ./Processed_5TNET_R21/$fq_R2\_R12_uniq_3end_Negative.wig
##########
rm $fq_R1\_cutadapt.fastq $fq_R2\_cutadapt.fastq $fq_R2\_R12_cutadapt.fastq
rm $fq_R2\_R12_cutadapt_remD.fastq $fq_R2\_R12_cutadapt_remD_6Ntrim.fastq
