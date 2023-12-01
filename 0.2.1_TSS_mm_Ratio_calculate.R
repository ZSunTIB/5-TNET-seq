###############################
library("data.table")
library("RcppRoll")
library("stringr")
library("dplyr")
library("MALDIquant")
library("Biostrings")
library("tidyverse")
####################################################################################################################
#This program is used to map the reads containing mismatches in the 1st, 2nd or 3rd coordinates to the corresponding TSSs
#and calculate the slippage ratio for each TSS region.
####################################################################################################################
#Import the TSSs identified by 5'TNET-seq and dRNA-seq
TSS = read.csv("./Log phase TSSs.csv", header = TRUE, sep = ",")
TSS_plus = TSS[(TSS[,2] == "+"),]
TSS_minus = TSS[(TSS[,2] == "-"),]

#Import the read counts of adjusted 5' end
wig_5end_plus = fread("./Mismatch_1to3_extraction/Raw_fastq_R2_R12_uniq_adjust_Positive.wig")
wig_5end_minus = fread("./Mismatch_1to3_extraction/Raw_fastq_R2_R12_uniq_adjust_Negative.wig")

#Import the reads containing mismatches in the 1st, 2nd or 3rd coordinates
Mismatched_reads = read.delim("./Mismatch_1to3_extraction/Raw_fastq_R2_R12_uniq_1to3_mismatch_reads_info.txt", header = TRUE, sep = "\t", row.names = NULL)
####################################################################################################################
wig_5end_plus = wig_5end_plus[, .(V2,V3)]
wig_5end_minus = wig_5end_minus[, .(V2,V3)]
####################################################################################################################
#Remove the reads containing mismatches in the 1st, 2nd or 3rd coordinates that don't have downstream TSS
Mismatched_reads = Mismatched_reads[(Mismatched_reads$TSS_index <= 6),]
Mismatched_reads = Mismatched_reads[abs(Mismatched_reads[,]$X5end_coordinate) != abs(Mismatched_reads[,]$TSS_coor),]
Mismatched_reads = Mismatched_reads[,c(1,2,7)]
Mismatched_reads[,3] = abs(Mismatched_reads[,3])

#Count the slippage reads with same 5' end
Mismatched_reads_plus = as.data.frame(table(Mismatched_reads[(Mismatched_reads[,2] == "+"),3]))
Mismatched_reads_plus[,1] = as.numeric(levels(Mismatched_reads_plus[,1]))
Mismatched_reads_plus[,3] = Mismatched_reads_plus[,1]
colnames(Mismatched_reads_plus)[3] = "Real_5end"

for(i in 1:nrow(Mismatched_reads_plus))
{
  Mismatched_reads_plus[i, 4] = min(TSS_plus[(TSS_plus[, 1] >= Mismatched_reads_plus[i, 3]), 1]);
}

Mismatched_reads_plus[, 5] = Mismatched_reads_plus[, 3] - Mismatched_reads_plus[, 4]
colnames(Mismatched_reads_plus)[4:5] = c("Near_dn_TSS","Dis_near_dn_TSS")

Mismatched_reads_minus = as.data.frame(table(Mismatched_reads[(Mismatched_reads[,2] == "-"),3]))
Mismatched_reads_minus[,1] = as.numeric(levels(Mismatched_reads_minus[,1]))
Mismatched_reads_minus[,3] = Mismatched_reads_minus[,1]
colnames(Mismatched_reads_minus)[3] = "Real_5end"

for(i in 1:nrow(Mismatched_reads_minus))
{
  Mismatched_reads_minus[i, 4] = max(TSS_minus[(TSS_minus[, 1] <= Mismatched_reads_minus[i, 3]), 1]);
}

Mismatched_reads_minus[, 5] = -(Mismatched_reads_minus[, 3] - Mismatched_reads_minus[, 4])
colnames(Mismatched_reads_minus)[4:5] = c("Near_dn_TSS","Dis_near_dn_TSS")
####################################################################################################################
#Generate count of the reads derived from each TSS
TSS_plus[,3] = wig_5end_plus[(TSS_plus[,1]),]$V3
Mismatched_reads_plus_D0 = Mismatched_reads_plus[(Mismatched_reads_plus[,5] == "0"),]
rownames(Mismatched_reads_plus_D0) = NULL
TSS_plus[,4] = 0
colnames(TSS_plus)[3:4] = c("TSS_RD","TSS_mm_RD")

for(i in 1:nrow(Mismatched_reads_plus_D0))
{
  TSS_plus[(TSS_plus[,1] == Mismatched_reads_plus_D0[i,3]),4] = Mismatched_reads_plus_D0[i,2]
}

TSS_minus[,3] = wig_5end_minus[(TSS_minus[,1]),]$V3
Mismatched_reads_minus_D0 = Mismatched_reads_minus[(Mismatched_reads_minus[,5] == "0"),]
rownames(Mismatched_reads_minus_D0) = NULL
TSS_minus[,4] = 0
colnames(TSS_minus)[3:4] = c("TSS_RD","TSS_mm_RD")

for(i in 1:nrow(Mismatched_reads_minus_D0))
{
  TSS_minus[(TSS_minus[,1] == Mismatched_reads_minus_D0[i,3]),4] = Mismatched_reads_minus_D0[i,2]
}

total_m1to3_indEL6_R1matchTSS = sum(Mismatched_reads_plus_D0[,2]) + sum(Mismatched_reads_minus_D0[,2])
##########
TSS = rbind(TSS_plus, TSS_minus)
TSS = TSS[order(TSS$TSS_coordinate),]
rownames(TSS) = NULL

#Calculate the slippage ratio
TSS[,5] = TSS[,4]/TSS[,3]*100
colnames(TSS)[5] = "TSS_mm_Ratio"
############################################################################################################################
#Remove intermediate file.
rm(Mismatched_reads_plus_D0,Mismatched_reads_minus_D0,TSS_plus,TSS_minus)
############################################################################################################################
write.table(TSS, file="./1.0.1_TSS_slippage_Ratio.csv", sep = ",", col.names = TRUE, row.names = FALSE, quote = FALSE)
