###############################
library("data.table")
library("RcppRoll")
library("stringr")
library("dplyr")
library("MALDIquant")
library("Biostrings")
library("tidyverse")
####################################################################################################################
#Import the adjusted 5'end read count and read depth
wig_plus = fread("Raw_fastq_R2_R12_uniq_adjust_Positive.wig")
wig_minus = fread("Raw_fastq_R2_R12_uniq_adjust_Negative.wig")

wig_plus_mTEX = fread("Raw_mT_fastq_R2_R12_uniq_adjust_Positive.wig")
wig_minus_mTEX = fread("Raw_mT_fastq_R2_R12_uniq_adjust_Negative.wig")


wig_plus_RD = fread("Raw_fastq_R2_R12_uniq_Positive_Allread_depth_adjust.wig")
wig_minus_RD = fread("Raw_fastq_R2_R12_uniq_Negative_Allread_depth_adjust.wig")

total_reads_pTEX = read.delim("../../1.1_6CLB_05pT_Rep2_Processed_5TNET_R21/Raw_fastq_R2_R12_after_alignment_read_length.txt", header = FALSE, sep = "\t", row.names = NULL)
total_reads_mTEX = read.delim("../../1.2_6CLB_05mT_Rep2_Processed_5TNET_R21/Raw_mT_fastq_R2_R12_after_alignment_read_length.txt", header = FALSE, sep = "\t", row.names = NULL)

total_reads_pTEX = sum(total_reads_pTEX$V2)
total_reads_mTEX = sum(total_reads_mTEX$V2)
####################################################################################################################
#Set the cutoffs to pick up TSSs
count_M = 5
score_cut = 20
ratio = 0.2
increase = 2.0
####################################################################################################################
wig_plus = wig_plus[, .(V2,V3)]
setnames(wig_plus,c("coordinate","count"))

setkey(wig_plus, coordinate)
wig_plus = wig_plus[ , Roll_median51 := roll_median(count, n = 51, weights = NULL, by = 1L, fill = NA, partial = FALSE,
                                                    align = c("center", "left", "right"), normalize = TRUE, na.rm = FALSE),]

wig_plus = wig_plus[, score := count / Roll_median51]
wig_plus = wig_plus[, score_F := score]
wig_plus = wig_plus[score %in% c("Inf", "NaN"), score_F := count]
wig_plus = wig_plus[, score_F_abs := score_F]
wig_plus = wig_plus[, Strand := "+"]
##########################
wig_plus_mTEX = wig_plus_mTEX[, .(V2,V3)]
setnames(wig_plus_mTEX,c("coordinate","count"))
wig_plus = wig_plus[, count_mTEX := wig_plus_mTEX$count]
##########################
wig_plus_RD = wig_plus_RD[, .(V2,V3)]
setnames(wig_plus_RD,c("coordinate","count"))
wig_plus_RD_fr = wig_plus_RD[,c(2)]
wig_plus_RD_fr_last = wig_plus_RD_fr[4639675,]
wig_plus_RD_fr_bind = rbind(wig_plus_RD_fr_last, wig_plus_RD_fr)
wig_plus_RD_fr_bind = wig_plus_RD_fr_bind[-4639676,]
wig_plus_RD = cbind(wig_plus_RD, wig_plus_RD_fr_bind)
setnames(wig_plus_RD,c("coordinate","count","count_pre"))
rm(wig_plus_RD_fr, wig_plus_RD_fr_last, wig_plus_RD_fr_bind)

wig_plus_RD = wig_plus_RD[, RD_increase := count/count_pre]
wig_plus_RD = wig_plus_RD[RD_increase %in% c("Inf"), RD_increase := count]
wig_plus_RD = wig_plus_RD[RD_increase %in% c("NaN"), RD_increase := 0]

wig_plus = wig_plus[, RD_increase := wig_plus_RD$RD_increase]
##########################
#Pick up the TSSs
TSS_plus = wig_plus[(count >= (total_reads_pTEX/1000000)*count_M) & (score_F >= score_cut) 
                    & ((count/total_reads_pTEX) >= ratio * (count_mTEX/total_reads_mTEX))
                    & (RD_increase >= increase)]
############################################################################################################################
wig_minus = wig_minus[, .(V2,V3)]
setnames(wig_minus,c("coordinate","count"))

setkey(wig_minus, coordinate)
wig_minus = wig_minus[ , Roll_median51 := roll_median(count, n = 51, weights = NULL, by = 1L, fill = NA, partial = FALSE,
                                                      align = c("center", "left", "right"), normalize = TRUE, na.rm = FALSE),]
wig_minus = wig_minus[, score := count / Roll_median51]
wig_minus = wig_minus[, score_F := (- score)]
wig_minus = wig_minus[score %in% c("Inf", "NaN"), score_F := (- count)]
wig_minus = wig_minus[, score_F_abs := (- score_F)]
wig_minus = wig_minus[, Strand := "-"]
##########################
wig_minus_mTEX = wig_minus_mTEX[, .(V2,V3)]
setnames(wig_minus_mTEX,c("coordinate","count"))
wig_minus = wig_minus[, count_mTEX := wig_minus_mTEX$count]
##########################
wig_minus_RD = wig_minus_RD[, .(V2,V3)]
setnames(wig_minus_RD,c("coordinate","count"))
wig_minus_RD_fr = wig_minus_RD[,c(2)]
wig_minus_RD_fr_first = wig_minus_RD_fr[1,]
wig_minus_RD_fr_bind = rbind(wig_minus_RD_fr, wig_minus_RD_fr_first)
wig_minus_RD_fr_bind = wig_minus_RD_fr_bind[-1,]
wig_minus_RD = cbind(wig_minus_RD, wig_minus_RD_fr_bind)
setnames(wig_minus_RD,c("coordinate","count","count_pre"))

rm(wig_minus_RD_fr, wig_minus_RD_fr_first, wig_minus_RD_fr_bind)

wig_minus_RD = wig_minus_RD[, RD_increase := count/count_pre]
wig_minus_RD = wig_minus_RD[RD_increase %in% c("Inf"), RD_increase := count]
wig_minus_RD = wig_minus_RD[RD_increase %in% c("NaN"), RD_increase := 0]
##########
wig_minus = wig_minus[, RD_increase := wig_minus_RD$RD_increase]
##########################
#Pick up the TSSs
TSS_minus = wig_minus[(count >= (total_reads_pTEX/1000000)*count_M) & (score_F <= -score_cut) 
                      & ((count/total_reads_pTEX) >= ratio * (count_mTEX/total_reads_mTEX))
                      & (RD_increase >= increase)]
#############################################################################################################################
TSS = rbind(TSS_plus, TSS_minus)

TSS = TSS[order(-score_F_abs)]
TSS$count = TSS$count*1000000/total_reads_pTEX
TSS$count_mTEX = TSS$count_mTEX*1000000/total_reads_mTEX
############################################################################################################################
#Remove TSSs in tRNA and rRNA regions.
trRNA = read.delim("0.0.1 trRNA_coordinate_region_NC_000913.2.txt",
                   header = FALSE, sep = "\t", fill = TRUE, row.names = NULL)

colnames(trRNA) = c("Genome","Start","End")

TSS = data.frame(TSS)
TSS[,(ncol(TSS)+1)] = 0

for (i in 1:nrow(TSS))
{
  if (any(TSS[i,1]>=trRNA[,2] & TSS[i,1]<=trRNA[,3]))
  {
    TSS[i,ncol(TSS)] = TSS[i,ncol(TSS)] + 1;
  }
  else
  {
  }
}

TSS_rm_trRNA = TSS[(TSS[,ncol(TSS)]==0),]

TSS_rm_trRNA[, 4] = TSS_rm_trRNA[, 8]
TSS_rm_trRNA = TSS_rm_trRNA[,c(1:7,9)]
colnames(TSS_rm_trRNA)[4] = "count_mTEX"

TSS[, 4] = TSS[, 8]
TSS = TSS[,c(1:7,9)]
colnames(TSS)[4] = "count_mTEX"

rownames(TSS_rm_trRNA) = NULL
write.table(TSS_rm_trRNA,  file="./Log phase TSSs_adjust.csv.csv",  sep = ",", col.names = TRUE, row.names = FALSE, quote = FALSE)



