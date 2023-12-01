#!/bin/sh
###############################################################################################
fq_R2=Raw_fastq_R2
###############################################################################################
#Extract reads mapped to both the sense and antisense strands, containing mismatches
#in the 1st, 2nd, or 3rd coordinates of the sense strand
#OR in the lAST 1st, 2nd, or 3rd coordinate of the antisense strand.

samtools view -h ../Processed_5TNET_R21/$fq_R2\_R12_uniq.bam | awk '{if ((($2 == "0") && \
($13 ~/MD:Z:0[ATCG]/ || $13 ~/MD:Z:1[ATCG]/ || $13 ~/MD:Z:2[ATCG]/)) || ($1 ~ /^@/)) {print $0} \
else if ((($2 == "16") && ($13 ~/[ATCG]0$/ || $13 ~/[ATCG]1$/ || $13 ~/[ATCG]2$/)) || ($1 ~ /^@/)) {print $0}}' \
> 0.0_$fq_R2\_R12_uniq_1to3_mismatch.sam


samtools view 0.0_$fq_R2\_R12_uniq_1to3_mismatch.sam -b -o 0.0_$fq_R2\_R12_uniq_1to3_mismatch.bam
samtools index 0.0_$fq_R2\_R12_uniq_1to3_mismatch.bam

#Classify the extracted reads with mismatches
samtools view 0.0_$fq_R2\_R12_uniq_1to3_mismatch.bam | awk -v OFS="\t" 'BEGIN {print "Mismatch_type","Strand","5_coordinate","Read_length","Read_sequence","Position","Qname"} \
{if(($2 == "0") && ($13 ~/MD:Z:0[ATCG]/) && ($13 !~/MD:Z:0[ATCG][0]/) && ($13 !~/MD:Z:0[ATCG][1][ATCG]/)) {print "1M","+",$4,length($10),$10,$4,$1} \
if (($2 == "0") && ($13 ~/MD:Z:1[ATCG]/) && ($13 !~/MD:Z:0[ATCG][0]/) && ($13 !~/MD:Z:1[ATCG][0][ATCG]/)) {print "2M","+",$4,length($10),$10,$4,$1} \
if (($2 == "0") && ($13 ~/MD:Z:2[ATCG]/)) {print "3M","+",$4,length($10),$10,$4,$1} \
if (($2 == "0") && ($13 ~/MD:Z:0[ATCG]0[ATCG]/)) {print "1M2M","+",$4,length($10),$10,$4,$1} \
if (($2 == "0") && ($13 ~/MD:Z:0[ATCG]1[ATCG]/)) {print "1M3M","+",$4,length($10),$10,$4,$1} \
if (($2 == "0") && ($13 ~/MD:Z:1[ATCG]0[ATCG]/)) {print "2M3M","+",$4,length($10),$10,$4,$1} \
if(($2 == "16") && ($13 ~/[ATCG]0$/) && ($13 !~/[ATCG][0][ATCG]0$/) && ($13 !~/[ATCG][1][ATCG]0$/)) {print "-1M","-",$4+length($10)-1,length($10),$10,$4,$1} \
if (($2 == "16") && ($13 ~/[ATCG]1$/) && ($13 !~/[ATCG][0][ATCG]0$/) && ($13 !~/[ATCG][0][ATCG]1$/)) {print "-2M","-",$4+length($10)-1,length($10),$10,$4,$1} \
if (($2 == "16") && ($13 ~/[ATCG]2$/)) {print "-3M","-",$4+length($10)-1,length($10),$10,$4,$1} \
if (($2 == "16") && ($13 ~/[ATCG]0[ATCG]0$/)) {print "-1M2M","-",$4+length($10)-1,length($10),$10,$4,$1} \
if (($2 == "16") && ($13 ~/[ATCG]1[ATCG]0$/)) {print "-1M3M","-",$4+length($10)-1,length($10),$10,$4,$1} \
if (($2 == "16") && ($13 ~/[ATCG]0[ATCG]1$/)) {print "-2M3M","-",$4+length($10)-1,length($10),$10,$4,$1}}' \
> 3.1.1_$fq_R2\_R12_uniq_1to3_mismatch_reads.txt

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

samtools view -h 0.0_$fq_R2\_R12_uniq_1to3_mismatch.bam | 
awk 'BEGIN {getline; print $0} \
{if ((($2 == "0") && ($10 ~/^N/ || $10 ~/^[ATCG]N/ || $10 ~/^[ATCG][ATCG]N/)) \
|| (($2 == "16") && ($10 ~/N$/ || $10 ~/N[ATCG]$/ || $10 ~/N[ATCG][ATCG]$/))) {} \
else {print $0}}' \
> 0.0_$fq_R2\_R12_uniq_1to3_mismatch_deleteN.sam

samtools view 0.0_$fq_R2\_R12_uniq_1to3_mismatch_deleteN.sam -b -o 0.0_$fq_R2\_R12_uniq_1to3_mismatch_deleteN.bam

mv 0.0_$fq_R2\_R12_uniq_1to3_mismatch_deleteN.sam 0.0_$fq_R2\_R12_uniq_1to3_mismatch.sam
mv 0.0_$fq_R2\_R12_uniq_1to3_mismatch_deleteN.bam 0.0_$fq_R2\_R12_uniq_1to3_mismatch.bam

samtools index 0.0_$fq_R2\_R12_uniq_1to3_mismatch.bam

rm 3.1.1_$fq_R2\_R12_uniq_1to3_mismatch_reads.txt \
3.1.3_negative_reads_seq.txt 3.1.4_negative_reads_rev_com_seq.txt \
3.1.5_$fq_R2\_R12_uniq_1to3_mismatch_reads.txt \
3.2.1_$fq_R2\_R12_uniq_1to3_mismatch_reads.txt

#Generate the 5' end reads count for the extracted reads with mismatches.
awk -v OFS="\t" 'NR==FNR {wig[$2]=$2"\t"$3;next} {if (($3 in wig) && ($2=="+")) {print $2"\t"wig[$3]}}' \
../Processed_5TNET_R21/$fq_R2\_R12_uniq_Positive.wig 3.2.2_$fq_R2\_R12_uniq_1to3_mismatch_reads_remN.txt \
> 3.3.1_$fq_R2\_R12_uniq_1to3_mismatch_reads_remN_posRD.txt

awk -v OFS="\t" 'NR==FNR {file1[$1"\t"$2]=$3;next} \
{if($2=="+") {print $0"\t"file1[$2"\t"$3]} else {print $0}}' \
3.3.1_$fq_R2\_R12_uniq_1to3_mismatch_reads_remN_posRD.txt 3.2.2_$fq_R2\_R12_uniq_1to3_mismatch_reads_remN.txt \
> 3.4.1_$fq_R2\_R12_uniq_1to3_mismatch_reads_posRD.txt

awk -v OFS="\t" 'NR==FNR {wig[$2]=$2"\t"$3;next} {if (($3 in wig) && ($2=="-")) {print $2"\t"wig[$3]}}' \
../Processed_5TNET_R21/$fq_R2\_R12_uniq_Negative.wig 3.2.2_$fq_R2\_R12_uniq_1to3_mismatch_reads_remN.txt \
> 3.3.2_$fq_R2\_R12_uniq_1to3_mismatch_reads_remN_negRD.txt

awk -v OFS="\t" 'NR==FNR {file1[$1"\t"$2]=$3;next} \
{if($2=="-") {print $0"\t"file1[$2"\t"$3]} else {print $0}}' \
3.3.2_$fq_R2\_R12_uniq_1to3_mismatch_reads_remN_negRD.txt 3.4.1_$fq_R2\_R12_uniq_1to3_mismatch_reads_posRD.txt \
> 3.4.2_$fq_R2\_R12_uniq_1to3_mismatch_reads_RD.txt

awk -v OFS="\t" 'NR==1 {print $0,"Read_count";next} {print $0}' 3.4.2_$fq_R2\_R12_uniq_1to3_mismatch_reads_RD.txt \
> 3.4.3_$fq_R2\_R12_uniq_1to3_mismatch_reads_RD.txt

rm 3.2.2_$fq_R2\_R12_uniq_1to3_mismatch_reads_remN.txt \
3.3.1_$fq_R2\_R12_uniq_1to3_mismatch_reads_remN_posRD.txt \
3.3.2_$fq_R2\_R12_uniq_1to3_mismatch_reads_remN_negRD.txt \
3.4.1_$fq_R2\_R12_uniq_1to3_mismatch_reads_posRD.txt \
3.4.2_$fq_R2\_R12_uniq_1to3_mismatch_reads_RD.txt

#Generate the 5' end reads count of each coordiante for the extracted reads with mismatches.
bedtools genomecov -d -5 -strand + -ibam 0.0_$fq_R2\_R12_uniq_1to3_mismatch.bam \
> $fq_R2\_R12_uniq_m3to0_Positive.wig
bedtools genomecov -d -5 -strand - -ibam 0.0_$fq_R2\_R12_uniq_1to3_mismatch.bam \
> $fq_R2\_R12_uniq_m3to0_Negative.wig
##########
#Find the real 5' end of the extracted reads with mismatches after removing the mismatched nucleotides.
#The first nucleotide downstream of the last mismatched nucleotide, which is identical with the first nucleotide of the read, is treated as the real 5' end.
dn=0
up=30

#Extract DNA sequence from the E. coli genome NC_000913.2
awk -v OFS="\t" 'BEGIN{getline;} \
{if($2=="+") {print "NC_000913.2",$3-"'$dn'"-1,$3+"'$up'",$3,1,"+"} \
else {print "NC_000913.2",$3-"'$up'"-1,$3+"'$dn'",$3,1,"-"}}' \
3.4.3_$fq_R2\_R12_uniq_1to3_mismatch_reads_RD.txt \
> 4.0_TSS+4_coor_tem.bed

bedtools getfasta -fi /Users/sz/Data_analysis/Genome/NC_000913.2.fasta \
-bed 4.0_TSS+4_coor_tem.bed \
-s -name -fo 4.0_TSS+4_coor_tem.fasta

awk "NR%2==0" 4.0_TSS+4_coor_tem.fasta \
> 4.0_TSS+4_coor_ext_tem.fasta
mv 4.0_TSS+4_coor_ext_tem.fasta 4.0_TSS+4_coor_tem.fasta

awk -v OFS="\t" 'BEGIN{print "5end_+30"} {print $0}' \
4.0_TSS+4_coor_tem.fasta > 4.1_TSS+4_coor_tem.fasta

paste -d "\t" 3.4.3_$fq_R2\_R12_uniq_1to3_mismatch_reads_RD.txt 4.1_TSS+4_coor_tem.fasta \
> 3.4.3.1_$fq_R2\_R12_uniq_1to3_mismatch_reads_RD_5end+30.txt

mv 3.4.3.1_$fq_R2\_R12_uniq_1to3_mismatch_reads_RD_5end+30.txt 3.4.3_$fq_R2\_R12_uniq_1to3_mismatch_reads_RD.txt
rm 4.0_TSS+4_coor_tem.bed 4.0_TSS+4_coor_tem.fasta 4.1_TSS+4_coor_tem.fasta

#Extract the first nucleotide of each read sequence, and the seuqences after the last mismatched nucleotide.
awk -v OFS="\t" 'NR==1 {print $1,$2,$3,$4,$5,$6,$7,"1st_string","Dn_mm_strings";next} \
{if($1 == "1M") {print $1,$2,$3,$4,$5,$6,$7,substr($5,1,1),substr($8,2)} \
if($1 == "2M") {print $1,$2,$3,$4,$5,$6,$7,substr($5,1,1),substr($8,3)} \
if($1 == "3M") {print $1,$2,$3,$4,$5,$6,$7,substr($5,1,1),substr($8,4)} \
if($1 == "1M2M") {print $1,$2,$3,$4,$5,$6,$7,substr($5,1,1),substr($8,3)} \
if($1 == "1M3M") {print $1,$2,$3,$4,$5,$6,$7,substr($5,1,1),substr($8,4)} \
if($1 == "2M3M") {print $1,$2,$3,$4,$5,$6,$7,substr($5,1,1),substr($8,4)} \
if($1 == "-1M") {print $1,$2,$3,$4,$5,$6,$7,substr($5,1,1),substr($8,2)} \
if($1 == "-2M") {print $1,$2,$3,$4,$5,$6,$7,substr($5,1,1),substr($8,3)} \
if($1 == "-3M") {print $1,$2,$3,$4,$5,$6,$7,substr($5,1,1),substr($8,4)} \
if($1 == "-1M2M") {print $1,$2,$3,$4,$5,$6,$7,substr($5,1,1),substr($8,3)} \
if($1 == "-1M3M") {print $1,$2,$3,$4,$5,$6,$7,substr($5,1,1),substr($8,4)} \
if($1 == "-2M3M") {print $1,$2,$3,$4,$5,$6,$7,substr($5,1,1),substr($8,4)} \
}' \
3.4.3_$fq_R2\_R12_uniq_1to3_mismatch_reads_RD.txt \
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
> 3.4.7_$fq_R2\_R12_uniq_1to3_mismatch_reads_RD_TSS.txt

rm 3.4.3_$fq_R2\_R12_uniq_1to3_mismatch_reads_RD.txt 3.4.4_$fq_R2\_R12_uniq_1to3_mismatch_reads_RD_substr.txt
rm 3.4.5_$fq_R2\_R12_uniq_1to3_mismatch_reads_RD_index.txt 3.4.6_$fq_R2\_R12_uniq_1to3_mismatch_reads_RD_TSS.txt

awk -v OFS="\t" 'NR==1 {print $0,"TSS_seq","TSS+1_seq","R1to2_seq","TSS1to2_seq",
"R1to3_seq", "TSS1to3_seq","R1to4_seq", "TSS1to4_seq";next} \
{print $0,substr($5,$8,1),substr($5,$8+1,1),substr($5,1,2),substr($5,$8,2), \
substr($5,1,3),substr($5,$8,3),substr($5,1,4),substr($5,$8,4)}' \
3.4.7_$fq_R2\_R12_uniq_1to3_mismatch_reads_RD_TSS.txt \
> 3.4.8_$fq_R2\_R12_uniq_1to3_mismatch_reads_RD_seq.txt

#Extract the DNA sequence from the real 5' end (0) to +4
dn=0
up=4

awk -v OFS="\t" 'BEGIN{getline;} \
{if($2=="+") {print "NC_000913.2",$7-"'$dn'"-1,$7+"'$up'",$7,1,"+"} \
else {print "NC_000913.2",-$7-"'$up'"-1,-$7+"'$dn'",-$7,1,"-"}}' \
3.4.8_$fq_R2\_R12_uniq_1to3_mismatch_reads_RD_seq.txt \
> 4.0_TSS+4_coor_tem.bed

bedtools getfasta -fi /Users/sz/Data_analysis/Genome/NC_000913.2.fasta \
-bed 4.0_TSS+4_coor_tem.bed \
-s -name -fo 4.0_TSS+4_coor_tem.fasta

awk "NR%2==0" 4.0_TSS+4_coor_tem.fasta \
> 4.0_TSS+4_coor_ext_tem.fasta
mv 4.0_TSS+4_coor_ext_tem.fasta 4.0_TSS+4_coor_tem.fasta

awk -v OFS="\t" 'BEGIN{print "TSS_+4"} {print $0}' \
4.0_TSS+4_coor_tem.fasta > 4.1_TSS+4_coor_tem.fasta

paste -d "\t" 3.4.8_$fq_R2\_R12_uniq_1to3_mismatch_reads_RD_seq.txt 4.1_TSS+4_coor_tem.fasta \
> 3.4.8.1_$fq_R2\_R12_uniq_1to3_mismatch_reads_RD_seq_TSS+4.txt

awk -v OFS="\t" 'NR==1 {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16;next} \
{print $1,$2,$3,$4,$5,$6,$7,$8,substr($17,1,1),substr($17,2,1),$11,substr($17,1,2), \
$13,substr($17,1,3),$15,substr($17,1,4)}' \
3.4.8.1_$fq_R2\_R12_uniq_1to3_mismatch_reads_RD_seq_TSS+4.txt \
> 3.4.8.2_$fq_R2\_R12_uniq_1to3_mismatch_reads_RD_seq.txt

mv 3.4.8.2_$fq_R2\_R12_uniq_1to3_mismatch_reads_RD_seq.txt 3.4.8_$fq_R2\_R12_uniq_1to3_mismatch_reads_RD_seq.txt
rm 3.4.8.1_$fq_R2\_R12_uniq_1to3_mismatch_reads_RD_seq_TSS+4.txt
rm 4.0_TSS+4_coor_tem.bed 4.0_TSS+4_coor_tem.fasta 4.1_TSS+4_coor_tem.fasta

#Assign the reads to g1, g2, g3 and g0
awk -v OFS="\t" 'BEGIN {getline; print $0,"group"} \
{if($9 == $10) {print $0,"g1"} \
else {if(($9 != $10) && ($8 >=5)  && ($15 == $16)) {print $0,"g4"} \
else {if(($9 != $10) && ($8 >= 4)  && ($13 == $14)) {print $0,"g3"} \
else {if(($9 != $10) && ($8 >= 3)  && ($11 == $12)) {print $0,"g2"} \
else {print $0,"g0"}}}}
}' \
3.4.8_$fq_R2\_R12_uniq_1to3_mismatch_reads_RD_seq.txt \
> 3.4.9_$fq_R2\_R12_uniq_1to3_mismatch_reads_RD_group.txt

rm 3.4.7_$fq_R2\_R12_uniq_1to3_mismatch_reads_RD_TSS.txt 3.4.8_$fq_R2\_R12_uniq_1to3_mismatch_reads_RD_seq.txt
mv 3.4.9_$fq_R2\_R12_uniq_1to3_mismatch_reads_RD_group.txt 3.4.4_$fq_R2\_R12_uniq_1to3_mismatch_reads_RD_adjust_5end_all.txt

awk -v OFS="\t" 'BEGIN {getline; print $0} \
{if(($8 <= 6) && ($8 >= 2)) {print $0}}' \
3.4.4_$fq_R2\_R12_uniq_1to3_mismatch_reads_RD_adjust_5end_all.txt \
> 3.4.4_$fq_R2\_R12_uniq_1to3_mismatch_reads_RD_adjust_5end.txt

#Generate the adjusted 5' end reads count for the extracted reads with mismatches
awk -v OFS="\t" 'BEGIN {getline;} {number[$7]++} END \
{for(coordinate in number) {print coordinate,number[coordinate]}}' \
3.4.4_$fq_R2\_R12_uniq_1to3_mismatch_reads_RD_adjust_5end.txt \
> 3.4.5_$fq_R2\_R12_uniq_1to3_mismatch_reads_RD_adjust_5end_count.txt

awk -v OFS="\t" '{if($1 > 0) {print $0}}' \
3.4.5_$fq_R2\_R12_uniq_1to3_mismatch_reads_RD_adjust_5end_count.txt \
> 3.4.5.1_$fq_R2\_R12_uniq_1to3_mismatch_reads_RD_adjust_5end_count_Positive.txt

awk -v OFS="\t" '{if($1 < 0) {print $0}}' \
3.4.5_$fq_R2\_R12_uniq_1to3_mismatch_reads_RD_adjust_5end_count.txt \
> 3.4.5.2_$fq_R2\_R12_uniq_1to3_mismatch_reads_RD_adjust_5end_count_Negative.txt

#Generate the 5' end reads count of each coordinate for the extracted reads with mismatches
awk -v OFS="\t" 'NR==FNR {wig[$2]=$3;next} {print $0"\t"wig[$2]}' \
$fq_R2\_R12_uniq_m3to0_Positive.wig \
../Processed_5TNET_R21/$fq_R2\_R12_uniq_Positive.wig \
> 4.1_$fq_R2\_R12_uniq_plus_m3to0_Positive.wig

awk -v OFS="\t" 'NR==FNR {wig[$2]=$3;next} {print $0"\t"wig[$2]}' \
$fq_R2\_R12_uniq_m3to0_Negative.wig \
../Processed_5TNET_R21/$fq_R2\_R12_uniq_Negative.wig \
> 4.1_$fq_R2\_R12_uniq_plus_m3to0_Negative.wig

awk -v OFS="\t" 'NR==FNR {wig[$1]=$2;cor[$1]=$1;next} \
{if($2 in cor) {print $0"\t"wig[$2]} else {print $0"\t"0}}' \
3.4.5.1_$fq_R2\_R12_uniq_1to3_mismatch_reads_RD_adjust_5end_count_Positive.txt \
4.1_$fq_R2\_R12_uniq_plus_m3to0_Positive.wig \
> 4.3.1_$fq_R2\_R12_uniq_plus_m3to0_Positive_adjust.wig

awk -v OFS="\t" 'NR==FNR {wig[-$1]=$2;cor[-$1]=$1;next} \
{if($2 in cor) {print $0"\t"wig[$2]} else {print $0"\t"0}}' \
3.4.5.2_$fq_R2\_R12_uniq_1to3_mismatch_reads_RD_adjust_5end_count_Negative.txt \
4.1_$fq_R2\_R12_uniq_plus_m3to0_Negative.wig \
> 4.3.2_$fq_R2\_R12_uniq_plus_m3to0_Negative_adjust.wig

#Generate the adjusted 5' end reads count of each coordinate for the extracted reads with mismatches
awk -v OFS="\t" '{print $1,$2,$3-$4+$5}' 4.3.1_$fq_R2\_R12_uniq_plus_m3to0_Positive_adjust.wig \
> 4.4.1_$fq_R2\_R12_uniq_adjust_Positive.wig

awk -v OFS="\t" '{print $1,$2,$3-$4+$5}' 4.3.2_$fq_R2\_R12_uniq_plus_m3to0_Negative_adjust.wig \
> 4.4.2_$fq_R2\_R12_uniq_adjust_Negative.wig

rm 3.4.5_$fq_R2\_R12_uniq_1to3_mismatch_reads_RD_adjust_5end_count.txt
rm 3.4.5.1_$fq_R2\_R12_uniq_1to3_mismatch_reads_RD_adjust_5end_count_Positive.txt 3.4.5.2_$fq_R2\_R12_uniq_1to3_mismatch_reads_RD_adjust_5end_count_Negative.txt
rm 4.1_$fq_R2\_R12_uniq_plus_m3to0_Negative.wig 4.1_$fq_R2\_R12_uniq_plus_m3to0_Positive.wig
rm 4.3.1_$fq_R2\_R12_uniq_plus_m3to0_Positive_adjust.wig 4.3.2_$fq_R2\_R12_uniq_plus_m3to0_Negative_adjust.wig
rm 0.0_$fq_R2\_R12_uniq_1to3_mismatch.sam 3.4.4_$fq_R2\_R12_uniq_1to3_mismatch_reads_RD_adjust_5end_all.txt

mv 4.4.1_$fq_R2\_R12_uniq_adjust_Positive.wig $fq_R2\_R12_uniq_adjust_Positive.wig
mv 4.4.2_$fq_R2\_R12_uniq_adjust_Negative.wig $fq_R2\_R12_uniq_adjust_Negative.wig
mv 3.4.4_$fq_R2\_R12_uniq_1to3_mismatch_reads_RD_adjust_5end.txt $fq_R2\_R12_uniq_1to3_mismatch_reads_info.txt
mv 0.0_$fq_R2\_R12_uniq_1to3_mismatch.bam $fq_R2\_R12_uniq_1to3_mismatch.bam
mv 0.0_$fq_R2\_R12_uniq_1to3_mismatch.bam.bai $fq_R2\_R12_uniq_1to3_mismatch.bam.bai
