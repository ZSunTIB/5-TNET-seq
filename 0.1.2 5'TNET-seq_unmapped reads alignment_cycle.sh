#!/bin/sh
##########################################################################################################
fq_R2=Raw_fastq_R2
##########################################################################################################
#This program is used to cut the nucleotide of unmapped reads from 5' end and 
#align the residual reads to the genome.

cd ./Processed_5TNET_R21_unmapped
##########################################################################################################
#Remove reads that are shorter than 12-nt
seqkit seq -m 12 $fq_R2\_R12_unmapped.fastq -o $fq_R2\_R12_unmapped_m12.fastq
mv $fq_R2\_R12_unmapped_m12.fastq $fq_R2\_R12_unmapped.fastq

cp $fq_R2\_R12_unmapped.fastq $fq_R2\_R12_unmapped_c0_unmapped.fastq

for((cut=1; cut<=40; cut++))
do
	T=$((cut-1))
	#Collect the read name, 5' end nucleotide and whole sequence.
    awk 'NR%4==1 || NR%4==2 {print $0}' $fq_R2\_R12_unmapped_c$T\_unmapped.fastq \
    > 1.1_$fq_R2\_R12_unmapped_c$cut\_unmapped_2row.fastq
    
    awk 'NR%2==1 {print $0} NR%2==0 {print substr($0,1,1)}' 1.1_$fq_R2\_R12_unmapped_c$cut\_unmapped_2row.fastq \
    > 1.2_$fq_R2\_R12_unmapped_c$cut\_unmapped_5end1.fastq
    
    awk -v OFS="\t" 'NR%2 {a=$0;next} {print a,$0}' 1.2_$fq_R2\_R12_unmapped_c$cut\_unmapped_5end1.fastq \
    > 1.3_$fq_R2\_R12_unmapped_c$cut\_unmapped_5end1_merge.fastq
    ####################
	#Remove the 5' end nucleotide.
    cutadapt --cut 1 -o \
	$fq_R2\_R12_unmapped_c$cut.fastq \
	$fq_R2\_R12_unmapped_c$T\_unmapped.fastq

	#Align the residual read sequence to the genome and collect uniquely mapped reads.
	bowtie /Users/sz/Data_analysis/Genome/Bowtie_genome/MG1655_v2 -q \
	$fq_R2\_R12_unmapped_c$cut.fastq \
	--best --strata -m 1 -S $fq_R2\_R12_unmapped_c$cut\_uniq.sam
	
	samtools view -F 4 $fq_R2\_R12_unmapped_c$cut\_uniq.sam \
	-b -o $fq_R2\_R12_unmapped_c$cut\_uniq.bam
	
	bedtools bamtofastq -i $fq_R2\_R12_unmapped_c$cut\_uniq.bam \
	-fq $fq_R2\_R12_unmapped_c$cut\_uniq.fastq
	
	clumpify.sh in=$fq_R2\_R12_unmapped_c$cut\_uniq.fastq \
	out=$fq_R2\_R12_unmapped_c$cut\_uniq_rem.fastq subs=0 dedupe=t

	mv $fq_R2\_R12_unmapped_c$cut\_uniq_rem.fastq $fq_R2\_R12_unmapped_c$cut\_uniq.fastq

	###
	#Collect the read name and sequence.
	awk 'NR%4==1 || NR%4==2 {print $0}' $fq_R2\_R12_unmapped_c$cut\_uniq.fastq \
	> 1.1_$fq_R2\_R12_unmapped_c$cut\_uniq_2row.fastq
	
	awk -v OFS="\t" 'NR%2 {a=$0;next} {print a,$0}' 1.1_$fq_R2\_R12_unmapped_c$cut\_uniq_2row.fastq \
	> 1.3_$fq_R2\_R12_unmapped_c$cut\_uniq_5end1_merge.fastq

	awk -v OFS="\t" 'NR==FNR {reads[$1]=$2;next} {print $0,reads[$1]}' \
	1.3_$fq_R2\_R12_unmapped_c$cut\_uniq_5end1_merge.fastq \
	1.3_$fq_R2\_R12_unmapped_c$cut\_unmapped_5end1_merge.fastq \
	> 2.1_$fq_R2\_R12_unmapped_c$cut\_uniq_all.txt

	#Pick up the uniquely mapped reads after cutting one nucleotide from 5' end.
	awk -v OFS="\t" 'NF==3 {print $0}' 2.1_$fq_R2\_R12_unmapped_c$cut\_uniq_all.txt \
	> 2.0_$fq_R2\_R12_unmapped_c$cut\_uniq.txt
	##########
	#Align the residual read sequence to the genome, collect multiple mapped and unmapped reads.
	bowtie /Users/sz/Data_analysis/Genome/Bowtie_genome/MG1655_v2 -q \
	$fq_R2\_R12_unmapped_c$cut.fastq \
	-S $fq_R2\_R12_unmapped_c$cut\_multiple.sam
	
	samtools view -F 4 $fq_R2\_R12_unmapped_c$cut\_multiple.sam \
	-b -o $fq_R2\_R12_unmapped_c$cut\_multiple.bam

	bedtools bamtofastq -i $fq_R2\_R12_unmapped_c$cut\_multiple.bam \
	-fq $fq_R2\_R12_unmapped_c$cut\_multiple.fastq

	clumpify.sh in=$fq_R2\_R12_unmapped_c$cut\_multiple.fastq \
	out=$fq_R2\_R12_unmapped_c$cut\_multiple_rem.fastq subs=0 dedupe=t

	mv $fq_R2\_R12_unmapped_c$cut\_multiple_rem.fastq $fq_R2\_R12_unmapped_c$cut\_multiple.fastq

	###
	#Collect the read name and sequence.
	awk 'NR%4==1 || NR%4==2 {print $0}' $fq_R2\_R12_unmapped_c$cut\_multiple.fastq \
	> 1.1_$fq_R2\_R12_unmapped_c$cut\_multiple_2row.fastq

	awk -v OFS="\t" 'NR%2 {a=$0;next} {print a,$0}' 1.1_$fq_R2\_R12_unmapped_c$cut\_multiple_2row.fastq \
	> 1.3_$fq_R2\_R12_unmapped_c$cut\_multiple_5end1_merge.fastq

	awk -v OFS="\t" 'NR==FNR {reads[$1]=$2;next} {print $0,reads[$1]}' \
	1.3_$fq_R2\_R12_unmapped_c$cut\_multiple_5end1_merge.fastq \
	1.3_$fq_R2\_R12_unmapped_c$cut\_unmapped_5end1_merge.fastq \
	> 2.1_$fq_R2\_R12_unmapped_c$cut\_multiple_all.txt

	#Pick up unmapped reads after cutting one nucleotide from 5' end .
	awk -v OFS="\t" 'NF==2 {print $0}' 2.1_$fq_R2\_R12_unmapped_c$cut\_multiple_all.txt \
	> 2.0_$fq_R2\_R12_unmapped_c$cut\_unmapped.txt
	##########
	#Remove all mapped reads including multiple mapping.
	samtools view -f 4 $fq_R2\_R12_unmapped_c$cut\_multiple.sam \
	-b -o $fq_R2\_R12_unmapped_c$cut\_unmapped.bam

	bedtools bamtofastq -i $fq_R2\_R12_unmapped_c$cut\_unmapped.bam \
	-fq $fq_R2\_R12_unmapped_c$cut\_unmapped.fastq

	clumpify.sh in=$fq_R2\_R12_unmapped_c$cut\_unmapped.fastq \
	out=$fq_R2\_R12_unmapped_c$cut\_unmapped_rem.fastq subs=0 dedupe=t

	mv $fq_R2\_R12_unmapped_c$cut\_unmapped_rem.fastq $fq_R2\_R12_unmapped_c$cut\_unmapped.fastq
	##########
	rm 1.1_$fq_R2\_R12_unmapped_c$cut\_unmapped_2row.fastq 1.1_$fq_R2\_R12_unmapped_c$cut\_multiple_2row.fastq
	rm 1.1_$fq_R2\_R12_unmapped_c$cut\_uniq_2row.fastq 1.2_$fq_R2\_R12_unmapped_c$cut\_unmapped_5end1.fastq
	rm 1.3_$fq_R2\_R12_unmapped_c$cut\_multiple_5end1_merge.fastq
	rm 1.3_$fq_R2\_R12_unmapped_c$cut\_uniq_5end1_merge.fastq
	rm 2.1_$fq_R2\_R12_unmapped_c$cut\_uniq_all.txt 2.1_$fq_R2\_R12_unmapped_c$cut\_multiple_all.txt

	rm $fq_R2\_R12_unmapped_c$cut\_uniq.sam $fq_R2\_R12_unmapped_c$cut\_uniq.fastq
	rm $fq_R2\_R12_unmapped_c$cut\_multiple.sam $fq_R2\_R12_unmapped_c$cut\_multiple.bam
	rm $fq_R2\_R12_unmapped_c$cut\_multiple.fastq $fq_R2\_R12_unmapped_c$T\_unmapped.fastq
	rm $fq_R2\_R12_unmapped_c$cut\_unmapped.bam $fq_R2\_R12_unmapped_c$cut.fastq
	##########
	if [[ $cut == 1 ]];then
		continue
	else
		#Combine the nucleotides cut from different cycles of 5' end trimming
		awk -v OFS="" 'NR==FNR {reads[$1]=$2;next} {print $0,reads[$1]}' \
		1.3_$fq_R2\_R12_unmapped_c$cut\_unmapped_5end1_merge.fastq \
		1.3_$fq_R2\_R12_unmapped_c$T\_unmapped_5end1_merge.fastq \
		> 1.3_$fq_R2\_R12_unmapped_c$cut\_unmapped_5end1_merge_int.fastq

		mv 1.3_$fq_R2\_R12_unmapped_c$cut\_unmapped_5end1_merge_int.fastq 1.3_$fq_R2\_R12_unmapped_c$cut\_unmapped_5end1_merge.fastq
		
		awk -v OFS="\t" 'NR==FNR {reads[$1]=$2;next} {print $1,reads[$1],$3}' \
		1.3_$fq_R2\_R12_unmapped_c$cut\_unmapped_5end1_merge.fastq \
		2.0_$fq_R2\_R12_unmapped_c$cut\_uniq.txt \
		> 2.0_$fq_R2\_R12_unmapped_c$cut\_uniq_int.txt

		mv 2.0_$fq_R2\_R12_unmapped_c$cut\_uniq_int.txt 2.0_$fq_R2\_R12_unmapped_c$cut\_uniq.txt
		
		cat 2.0_$fq_R2\_R12_unmapped_c$T\_uniq.txt 2.0_$fq_R2\_R12_unmapped_c$cut\_uniq.txt \
		> 2.0_$fq_R2\_R12_unmapped_c$cut\_uniq_int.txt
		mv 2.0_$fq_R2\_R12_unmapped_c$cut\_uniq_int.txt 2.0_$fq_R2\_R12_unmapped_c$cut\_uniq.txt

		awk -v OFS="\t" 'NR==FNR {reads[$1]=$2;next} {print $1,reads[$1],$3}' \
		1.3_$fq_R2\_R12_unmapped_c$cut\_unmapped_5end1_merge.fastq \
		2.0_$fq_R2\_R12_unmapped_c$cut\_unmapped.txt \
		> 2.0_$fq_R2\_R12_unmapped_c$cut\_unmapped_int.txt

		mv 2.0_$fq_R2\_R12_unmapped_c$cut\_unmapped_int.txt 2.0_$fq_R2\_R12_unmapped_c$cut\_unmapped.txt
		
		#Merge the bam files containing uniquely mapped reads generated from difference cycles.
		samtools merge $fq_R2\_R12_unmapped_c$cut\_uniq_int.bam \
		$fq_R2\_R12_unmapped_c$T\_uniq.bam $fq_R2\_R12_unmapped_c$cut\_uniq.bam

		mv $fq_R2\_R12_unmapped_c$cut\_uniq_int.bam \
		$fq_R2\_R12_unmapped_c$cut\_uniq.bam

		rm 1.3_$fq_R2\_R12_unmapped_c$T\_unmapped_5end1_merge.fastq
		rm 2.0_$fq_R2\_R12_unmapped_c$T\_unmapped.txt
		rm 2.0_$fq_R2\_R12_unmapped_c$T\_uniq.txt
		rm $fq_R2\_R12_unmapped_c$T\_uniq.bam

	fi
done

T=$((cut-1))

samtools sort $fq_R2\_R12_unmapped_c$T\_uniq.bam \
-o $fq_R2\_R12_unmapped_c$T\_uniq_sort.bam

samtools index $fq_R2\_R12_unmapped_c$T\_uniq_sort.bam

rm $fq_R2\_R12_unmapped_c$T\_uniq.bam $fq_R2\_R12_unmapped.bam
rm $fq_R2\_R12_unmapped_c$T\_unmapped.fastq $fq_R2\_R12_unmapped.fastq
rm 2.0_$fq_R2\_R12_unmapped_c$T\_unmapped.txt

mv 2.0_$fq_R2\_R12_unmapped_c$T\_uniq.txt $fq_R2\_R12_unmapped_c$T\_uniq.txt
mv $fq_R2\_R12_unmapped_c$T\_uniq_sort.bam $fq_R2\_R12_unmapped_c$T\_uniq.bam
mv $fq_R2\_R12_unmapped_c$T\_uniq_sort.bam.bai $fq_R2\_R12_unmapped_c$T\_uniq.bam.bai
###############################################################################################
#Generate the 3' end reads count of each coordiante for 5'end trimmed reads.
bedtools genomecov -d -3 -strand + -ibam $fq_R2\_R12_unmapped_c$T\_uniq.bam \
> $fq_R2\_unmapped_c$T\_uniq_3End_Positive.wig

bedtools genomecov -d -3 -strand - -ibam $fq_R2\_R12_unmapped_c$T\_uniq.bam \
> $fq_R2\_unmapped_c$T\_uniq_3End_Negative.wig
