#!/bin/sh
###############################################################################################
fq_R2=Raw_fastq_R2
###############################################################################################
#Extract reads mapped to both the sense and antisense strands, containing mismatches
#in the 1st, 2nd, or 3rd coordinates of the sense strand
#OR in the lAST 1st, 2nd, or 3rd coordinate of the antisense strand.

mkdir Bam_unmapped_align

cd ./Bam_unmapped_align

#Classify all reads aligned to the genome after cutting 5' nucleotides based on mismatch types.
cp ../Processed_5TNET_R21_unmapped/$fq_R2\_R12_unmapped_c40_uniq.bam \
$fq_R2\_R12_unmapped_c40_uniq.bam

cp ../Processed_5TNET_R21_unmapped/1.3_$fq_R2\_R12_unmapped_c40_unmapped_5end1_merge.fastq \
1.3_$fq_R2\_R12_unmapped_c40_unmapped_5end1_merge.fastq

#Classify the 5' end trimmed reads with mismatches
samtools view $fq_R2\_R12_unmapped_c40_uniq.bam | awk -v OFS="\t" 'BEGIN {print "Mismatch_type","Strand","5_coordinate","Read_length","Read_sequence","Position","Qname"} \
{if(($2 == "0") && ($13 ~/MD:Z:0[ATCG]/) && ($13 !~/MD:Z:0[ATCG][0]/) && ($13 !~/MD:Z:0[ATCG][1][ATCG]/)) {print "1M","+",$4,length($10),$10,$4,$1} \
if (($2 == "0") && ($13 ~/MD:Z:1[ATCG]/) && ($13 !~/MD:Z:0[ATCG][0]/) && ($13 !~/MD:Z:1[ATCG][0][ATCG]/)) {print "2M","+",$4,length($10),$10,$4,$1} \
if (($2 == "0") && ($13 ~/MD:Z:2[ATCG]/)) {print "3M","+",$4,length($10),$10,$4,$1} \
if (($2 == "0") && ($13 ~/MD:Z:0[ATCG]0[ATCG]/)) {print "1M2M","+",$4,length($10),$10,$4,$1} \
if (($2 == "0") && ($13 ~/MD:Z:0[ATCG]1[ATCG]/)) {print "1M3M","+",$4,length($10),$10,$4,$1} \
if (($2 == "0") && ($13 ~/MD:Z:1[ATCG]0[ATCG]/)) {print "2M3M","+",$4,length($10),$10,$4,$1} \
if (($2 == "16") && ($13 ~/[ATCG]0$/) && ($13 !~/[ATCG][0][ATCG]0$/) && ($13 !~/[ATCG][1][ATCG]0$/)) {print "-1M","-",$4+length($10)-1,length($10),$10,$4,$1} \
if (($2 == "16") && ($13 ~/[ATCG]1$/) && ($13 !~/[ATCG][0][ATCG]0$/) && ($13 !~/[ATCG][0][ATCG]1$/)) {print "-2M","-",$4+length($10)-1,length($10),$10,$4,$1} \
if (($2 == "16") && ($13 ~/[ATCG]2$/)) {print "-3M","-",$4+length($10)-1,length($10),$10,$4,$1} \
if (($2 == "16") && ($13 ~/[ATCG]0[ATCG]0$/)) {print "-1M2M","-",$4+length($10)-1,length($10),$10,$4,$1} \
if (($2 == "16") && ($13 ~/[ATCG]1[ATCG]0$/)) {print "-1M3M","-",$4+length($10)-1,length($10),$10,$4,$1} \
if (($2 == "16") && ($13 ~/[ATCG]0[ATCG]1$/)) {print "-2M3M","-",$4+length($10)-1,length($10),$10,$4,$1}}' \
> 3.1.1_$fq_R2\_R12_uniq_1to3_mismatch_reads.txt

#Remove the reads that don't have mutation in the 1st, 2nd, or 3rd coordinates.
samtools view -h $fq_R2\_R12_unmapped_c40_uniq.bam | \
awk '(($2 == "0") && ($13 ~/MD:Z:0[ATCG]/) && ($13 !~/MD:Z:0[ATCG][0]/) && ($13 !~/MD:Z:0[ATCG][1][ATCG]/)) || \
(($2 == "0") && ($13 ~/MD:Z:1[ATCG]/) && ($13 !~/MD:Z:0[ATCG][0]/) && ($13 !~/MD:Z:1[ATCG][0][ATCG]/)) || \
(($2 == "0") && ($13 ~/MD:Z:2[ATCG]/)) || \
(($2 == "0") && ($13 ~/MD:Z:0[ATCG]0[ATCG]/)) || \
(($2 == "0") && ($13 ~/MD:Z:0[ATCG]1[ATCG]/)) || \
(($2 == "0") && ($13 ~/MD:Z:1[ATCG]0[ATCG]/)) || \
(($2 == "16") && ($13 ~/[ATCG]0$/) && ($13 !~/[ATCG][0][ATCG]0$/) && ($13 !~/[ATCG][1][ATCG]0$/)) || \
(($2 == "16") && ($13 ~/[ATCG]1$/) && ($13 !~/[ATCG][0][ATCG]0$/) && ($13 !~/[ATCG][0][ATCG]1$/)) || \
(($2 == "16") && ($13 ~/[ATCG]2$/)) || \
(($2 == "16") && ($13 ~/[ATCG]0[ATCG]0$/)) || \
(($2 == "16") && ($13 ~/[ATCG]1[ATCG]0$/)) || \
(($2 == "16") && ($13 ~/[ATCG]0[ATCG]1$/)) || \
$1 ~ /^@/' | 
samtools view -bS > $fq_R2\_R12_unmapped_c40_uniq_sort_rem.bam

#Transform the read sequence mapped to the "-" strand to the reverse-complement sequence.
awk 'BEGIN {getline; print "Rec_com_sequence"} {if($2=="-") {print $5} else {print "+"}}' \
3.1.1_$fq_R2\_R12_uniq_1to3_mismatch_reads.txt \
> 3.1.3_negative_reads_seq.txt

perl -pe 'chomp;tr/ACGTacgt/TGCAtgca/;$_=reverse."\n"' 3.1.3_negative_reads_seq.txt \
> 3.1.4_negative_reads_rev_com_seq.txt

paste -d "\t" 3.1.1_$fq_R2\_R12_uniq_1to3_mismatch_reads.txt 3.1.4_negative_reads_rev_com_seq.txt \
> 3.1.5_$fq_R2\_R12_uniq_1to3_mismatch_reads.txt

awk -v OFS="\t" 'BEGIN {getline; print "Mismatch_type","Strand","5end_coordinate","Read_length","5to3_sequence","Qname"} \
{if($2=="+") {print $1,$2,$3,$4,$5,$7} \
else {print $1,$2,$3,$4,$8,$7}}' \
3.1.5_$fq_R2\_R12_uniq_1to3_mismatch_reads.txt \
> 3.2.1_$fq_R2\_R12_uniq_1to3_mismatch_reads.txt

#Remove the reads with "N" in the mismatched position
awk 'BEGIN {getline; print $0} \
{if ($5 ~/^N/ || $5 ~/^[ATCG]N/ || $5 ~/^[ATCG][ATCG]N/) {} \
else {print $0}}' 3.2.1_$fq_R2\_R12_uniq_1to3_mismatch_reads.txt \
> 3.2.2_$fq_R2\_R12_uniq_1to3_mismatch_reads_remN.txt

samtools view -h $fq_R2\_R12_unmapped_c40_uniq_sort_rem.bam | 
awk 'BEGIN {getline; print $0} \
{if ((($2 == "0") && ($10 ~/^N/ || $10 ~/^[ATCG]N/ || $10 ~/^[ATCG][ATCG]N/)) \
|| (($2 == "16") && ($10 ~/N$/ || $10 ~/N[ATCG]$/ || $10 ~/N[ATCG][ATCG]$/))) {} \
else {print $0}}' \
> $fq_R2\_R12_unmapped_c40_uniq_sort_deleteN.sam

samtools view $fq_R2\_R12_unmapped_c40_uniq_sort_deleteN.sam -b -o $fq_R2\_R12_unmapped_c40_uniq_sort_deleteN.bam

mv $fq_R2\_R12_unmapped_c40_uniq_sort_deleteN.sam $fq_R2\_R12_unmapped_c40_uniq_sort.sam
mv $fq_R2\_R12_unmapped_c40_uniq_sort_deleteN.bam $fq_R2\_R12_unmapped_c40_uniq.bam
rm $fq_R2\_R12_unmapped_c40_uniq_sort.sam

samtools index $fq_R2\_R12_unmapped_c40_uniq.bam

rm $fq_R2\_R12_unmapped_c40_uniq_sort_rem.bam
rm 3.1.1_$fq_R2\_R12_uniq_1to3_mismatch_reads.txt \
3.1.3_negative_reads_seq.txt 3.1.4_negative_reads_rev_com_seq.txt \
3.1.5_$fq_R2\_R12_uniq_1to3_mismatch_reads.txt \
3.2.1_$fq_R2\_R12_uniq_1to3_mismatch_reads.txt


awk -v OFS="\t" 'NR==FNR {wig[$1]=$2;next} {print $0"\t"wig["@"$6]}' \
1.3_$fq_R2\_R12_unmapped_c40_unmapped_5end1_merge.fastq 3.2.2_$fq_R2\_R12_uniq_1to3_mismatch_reads_remN.txt \
> 3.3.1_$fq_R2\_R12_uniq_1to3_mismatch_reads_remN_posRD.txt

awk -v OFS="\t" 'NR==1 {print $0,"5Cut_seq";next} {print $0}' 3.3.1_$fq_R2\_R12_uniq_1to3_mismatch_reads_remN_posRD.txt \
> 3.3.2_$fq_R2\_R12_uniq_1to3_mismatch_reads_remN_cut.txt

rm 3.3.1_$fq_R2\_R12_uniq_1to3_mismatch_reads_remN_posRD.txt 3.2.2_$fq_R2\_R12_uniq_1to3_mismatch_reads_remN.txt
##########
#Find the real 5' end of the extracted reads with mismatches after removing the mismatched nucleotides.
#The first nucleotide downstream of the last mismatched nucleotide, which is identical with the first nucleotide of the read, is treated as the real 5' end.
dn=0
up=15

#Extract DNA sequence from the E. coli genome NC_000913.2
awk -v OFS="\t" 'BEGIN{getline;} \
{if($2=="+") {print "NC_000913.2",$3-"'$dn'"-1,$3+"'$up'",$3,1,"+"} \
else {print "NC_000913.2",$3-"'$up'"-1,$3+"'$dn'",$3,1,"-"}}' \
3.3.2_$fq_R2\_R12_uniq_1to3_mismatch_reads_remN_cut.txt \
> 4.0_TSS+4_coor_tem.bed

bedtools getfasta -fi /Users/sz/Data_analysis/Genome/NC_000913.2.fasta \
-bed 4.0_TSS+4_coor_tem.bed \
-s -name -fo 4.0_TSS+4_coor_tem.fasta

awk "NR%2==0" 4.0_TSS+4_coor_tem.fasta \
> 4.0_TSS+4_coor_ext_tem.fasta

mv 4.0_TSS+4_coor_ext_tem.fasta 4.0_TSS+4_coor_tem.fasta

awk -v OFS="\t" 'BEGIN{print "5end_+30"} {print $0}' \
4.0_TSS+4_coor_tem.fasta > 4.1_TSS+4_coor_tem.fasta

paste -d "\t" 3.3.2_$fq_R2\_R12_uniq_1to3_mismatch_reads_remN_cut.txt 4.1_TSS+4_coor_tem.fasta \
> 3.3.2.1_$fq_R2\_R12_uniq_1to3_mismatch_reads_remN_cut_5end+30.txt

mv 3.3.2.1_$fq_R2\_R12_uniq_1to3_mismatch_reads_remN_cut_5end+30.txt 3.3.2_$fq_R2\_R12_uniq_1to3_mismatch_reads_remN_cut.txt
rm 4.0_TSS+4_coor_tem.bed 4.0_TSS+4_coor_tem.fasta 4.1_TSS+4_coor_tem.fasta

#Extract the first nucleotide of each read sequence, and the seuqences after the last mismatched nucleotide.
awk -v OFS="\t" 'NR==1 {print $1,$2,$3,$4,$5,$6,$7,"1st_string","Dn_mm_strings";next} \
{if($1 == "1M") {print $1,$2,$3,$4,$5,$6,$7,substr($7,1,1),substr($8,2)} \
if($1 == "2M") {print $1,$2,$3,$4,$5,$6,$7,substr($7,1,1),substr($8,3)} \
if($1 == "3M") {print $1,$2,$3,$4,$5,$6,$7,substr($7,1,1),substr($8,4)} \
if($1 == "1M2M") {print $1,$2,$3,$4,$5,$6,$7,substr($7,1,1),substr($8,3)} \
if($1 == "1M3M") {print $1,$2,$3,$4,$5,$6,$7,substr($7,1,1),substr($8,4)} \
if($1 == "2M3M") {print $1,$2,$3,$4,$5,$6,$7,substr($7,1,1),substr($8,4)} \
if($1 == "-1M") {print $1,$2,$3,$4,$5,$6,$7,substr($7,1,1),substr($8,2)} \
if($1 == "-2M") {print $1,$2,$3,$4,$5,$6,$7,substr($7,1,1),substr($8,3)} \
if($1 == "-3M") {print $1,$2,$3,$4,$5,$6,$7,substr($7,1,1),substr($8,4)} \
if($1 == "-1M2M") {print $1,$2,$3,$4,$5,$6,$7,substr($7,1,1),substr($8,3)} \
if($1 == "-1M3M") {print $1,$2,$3,$4,$5,$6,$7,substr($7,1,1),substr($8,4)} \
if($1 == "-2M3M") {print $1,$2,$3,$4,$5,$6,$7,substr($7,1,1),substr($8,4)} \
}' \
3.3.2_$fq_R2\_R12_uniq_1to3_mismatch_reads_remN_cut.txt \
> 3.4.4_$fq_R2\_R12_uniq_1to3_mismatch_reads_RD_substr.txt

#The first nucleotide downstream of the last mismatched nucleotide, which is identical with the first nucleotide of the read, is treated as the real 5' end.
awk -v OFS="\t" 'NR==1 {print $0,"TSS_index";next} {print $0,index($9,$8)}' \
3.4.4_$fq_R2\_R12_uniq_1to3_mismatch_reads_RD_substr.txt \
> 3.4.5_$fq_R2\_R12_uniq_1to3_mismatch_reads_RD_index.txt

awk -v OFS="\t" 'NR==1 {print $0,"TSS_coor","TSS_index";next} \
{if($1 == "1M") {print $0,$3+1-1+$10,1+$10} \
if($1 == "2M") {print $0,$3+2-1+$10,2+$10} \
if($1 == "3M") {print $0,$3+3-1+$10,3+$10} \
if($1 == "1M2M") {print $0,$3+2-1+$10,2+$10} \
if($1 == "1M3M") {print $0,$3+3-1+$10,3+$10} \
if($1 == "2M3M") {print $0,$3+3-1+$10,3+$10} \
if($1 == "-1M") {print $0,-($3-1+1-$10),-(-1-$10)} \
if($1 == "-2M") {print $0,-($3-2+1-$10),-(-2-$10)} \
if($1 == "-3M") {print $0,-($3-3+1-$10),-(-3-$10)} \
if($1 == "-1M2M") {print $0,-($3-2+1-$10),-(-2-$10)} \
if($1 == "-1M3M") {print $0,-($3-3+1-$10),-(-3-$10)} \
if($1 == "-2M3M") {print $0,-($3-3+1-$10),-(-3-$10)} \
}' \
3.4.5_$fq_R2\_R12_uniq_1to3_mismatch_reads_RD_index.txt \
> 3.4.6_$fq_R2\_R12_uniq_1to3_mismatch_reads_RD_TSS.txt

awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$7,$11,$12}' \
3.4.6_$fq_R2\_R12_uniq_1to3_mismatch_reads_RD_TSS.txt \
> $fq_R2\_R12_uniq_1to3_mismatch_reads_RD_TSS.txt

rm 3.3.2_$fq_R2\_R12_uniq_1to3_mismatch_reads_remN_cut.txt
rm 3.4.4_$fq_R2\_R12_uniq_1to3_mismatch_reads_RD_substr.txt
rm 3.4.5_$fq_R2\_R12_uniq_1to3_mismatch_reads_RD_index.txt 3.4.6_$fq_R2\_R12_uniq_1to3_mismatch_reads_RD_TSS.txt
mv 1.3_$fq_R2\_R12_unmapped_c40_unmapped_5end1_merge.fastq $fq_R2\_R12_unmapped_c40_unmapped_5end1_merge.fastq
