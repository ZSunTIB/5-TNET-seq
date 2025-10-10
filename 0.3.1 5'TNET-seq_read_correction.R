###############################
library("data.table")
library("RcppRoll")
library("stringr")
library("plyr")
library("dplyr")
library("MALDIquant")
library("Biostrings")
library("tidyverse")
####################################################################################################################
#Import the Raw read depth and reads information
wig_plus_RD = fread("Raw_fastq_R2_R12_uniq_Positive_Allread_depth.wig")
wig_minus_RD = fread("Raw_fastq_R2_R12_uniq_Negative_Allread_depth.wig")

Mismatched_reads = read.delim("Raw_fastq_R2_R12_uniq_1to3_mismatch_reads_info.txt", header = TRUE, sep = "\t", row.names = NULL)
####################################################################################################################
Mismatched_reads = Mismatched_reads[(Mismatched_reads$TSS_index <= 6),]
Mismatched_reads = Mismatched_reads[abs(Mismatched_reads[,]$X5end_coordinate) != abs(Mismatched_reads[,]$TSS_coor),]
Mismatched_reads = Mismatched_reads[,c(1,2,3,7)]
mm_reads_used = nrow(Mismatched_reads)

Mismatched_reads_plus = Mismatched_reads[(Mismatched_reads$Strand == "+"),]
Mismatched_reads_minus = Mismatched_reads[(Mismatched_reads$Strand == "-"),]

#Adjust the total read depth
wig_plus_RD_adjust = wig_plus_RD
wig_minus_RD_adjust = wig_minus_RD

wig_plus_RD_adjust[, Adjust := V3]
wig_minus_RD_adjust[, Adjust := V3]

count_plus = ddply(Mismatched_reads_plus, .(Mismatched_reads_plus$X5end_coordinate,Mismatched_reads_plus$TSS_coor),nrow)
count_minus = ddply(Mismatched_reads_minus, .(Mismatched_reads_minus$X5end_coordinate,Mismatched_reads_minus$TSS_coor),nrow)

tem1 = data.frame()
tem0 = data.frame()

for(i in 1:nrow(count_plus))
{
  tem1[1:(count_plus[i,2]-count_plus[i,1]),1] = count_plus[i,1]:(count_plus[i,2]-1)
  tem1[1:(count_plus[i,2]-count_plus[i,1]),2] = count_plus[i,3]

  tem0 = rbind(tem0,tem1)
  tem1 = data.frame()
}

count_ex_plus = tem0
rm(tem0, tem1)

count_ex_plus = count_ex_plus %>% group_by(V1) %>% summarise_all(sum)

#Calculate the adjusted read depth.
wig_plus_RD_adjust[count_ex_plus[,]$V1,]$Adjust = wig_plus_RD_adjust[count_ex_plus[,]$V1,]$Adjust - count_ex_plus[,]$V2
wig_plus_RD_adjust[, V3 := NULL]

#####
tem1 = data.frame()
tem0 = data.frame()

for(i in 1:nrow(count_minus))
  {
    tem1[1:(count_minus[i,2]+count_minus[i,1]),1] = count_minus[i,1]:(-count_minus[i,2]+1)
    tem1[1:(count_minus[i,2]+count_minus[i,1]),2] = count_minus[i,3]
    
    tem0 = rbind(tem0,tem1)
    tem1 = data.frame()
  }

count_ex_minus = tem0
rm(tem0, tem1)

count_ex_minus = count_ex_minus %>% group_by(V1) %>% summarise_all(sum)

#Calculate the adjusted read depth.
wig_minus_RD_adjust[count_ex_minus[,]$V1,]$Adjust = wig_minus_RD_adjust[count_ex_minus[,]$V1,]$Adjust - count_ex_minus[,]$V2
wig_minus_RD_adjust[, V3 := NULL]

rm(wig_plus_RD, wig_minus_RD,Mismatched_reads_plus,Mismatched_reads_minus)
####################################################################################################################
write.table(wig_plus_RD_adjust, file="Raw_fastq_R2_R12_uniq_Positive_Allread_depth_adjust.wig", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(wig_minus_RD_adjust, file="Raw_fastq_R2_R12_uniq_Negative_Allread_depth_adjust.wig", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)



