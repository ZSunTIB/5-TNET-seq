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
#This program is used to adjust the slippage ratio for each TSS region after trimming the 5' nucleotides.
####################################################################################################################
#Import the slippage ratio of each TSS region.
TSS_slippage = read.csv("./1.0.1_TSS_slippage_Ratio.csv", header = TRUE, sep = ",")

#Import the reads containing mismatches in the 1st, 2nd or 3rd coordinates after 5' end trimming
Cut_reads = read.delim("./Bam_unmapped_align/Raw_fastq_R2_R12_uniq_1to3_mismatch_reads_RD_TSS_2.txt", header = TRUE, sep = "\t", row.names = NULL)

#Import the reads containing mismatches in the 1st, 2nd or 3rd coordinates
Mismatched_reads = read.delim("./Mismatch_1to3_extraction/Raw_fastq_R2_R12_uniq_1to3_mismatch_reads_info_2.txt", header = TRUE, sep = "\t", row.names = NULL)
####################################################################################################################
#Remove the reads containing mismatches in the 1st, 2nd or 3rd coordinates that don't have downstream TSS
Cut_reads = Cut_reads[!(((Cut_reads$Mismatch_type == "1M") & (Cut_reads$TSS_index == 1)) | 
                          ((Cut_reads$Mismatch_type == "2M") & (Cut_reads$TSS_index == 2)) | 
                          ((Cut_reads$Mismatch_type == "3M") & (Cut_reads$TSS_index == 3)) | 
                          ((Cut_reads$Mismatch_type == "-1M") & (Cut_reads$TSS_index == 1)) | 
                          ((Cut_reads$Mismatch_type == "-2M") & (Cut_reads$TSS_index == 2)) | 
                          ((Cut_reads$Mismatch_type == "-3M") & (Cut_reads$TSS_index == 3)) | 
                          ((Cut_reads$Mismatch_type == "1M2M") & (Cut_reads$TSS_index == 2)) | 
                          ((Cut_reads$Mismatch_type == "1M3M") & (Cut_reads$TSS_index == 3)) | 
                          ((Cut_reads$Mismatch_type == "2M3M") & (Cut_reads$TSS_index == 3)) | 
                          ((Cut_reads$Mismatch_type == "-1M2M") & (Cut_reads$TSS_index == 2)) | 
                          ((Cut_reads$Mismatch_type == "-1M3M") & (Cut_reads$TSS_index == 3)) | 
                          ((Cut_reads$Mismatch_type == "-2M3M") & (Cut_reads$TSS_index == 3))),]

NC_000913.2 = readDNAStringSet("./NC_000913.2.fasta")

#Remove the reads containing mismatches in the 1st, 2nd or 3rd coordinates that don't have downstream TSS
Mismatched_reads = Mismatched_reads[!(((Mismatched_reads$Mismatch_type == "1M") & (Mismatched_reads$TSS_index == 1)) | 
                                     ((Mismatched_reads$Mismatch_type == "2M") & (Mismatched_reads$TSS_index == 2)) | 
                                     ((Mismatched_reads$Mismatch_type == "3M") & (Mismatched_reads$TSS_index == 3)) | 
                                     ((Mismatched_reads$Mismatch_type == "-1M") & (Mismatched_reads$TSS_index == 1)) | 
                                     ((Mismatched_reads$Mismatch_type == "-2M") & (Mismatched_reads$TSS_index == 2)) | 
                                     ((Mismatched_reads$Mismatch_type == "-3M") & (Mismatched_reads$TSS_index == 3)) | 
                                     ((Mismatched_reads$Mismatch_type == "1M2M") & (Mismatched_reads$TSS_index == 2)) | 
                                     ((Mismatched_reads$Mismatch_type == "1M3M") & (Mismatched_reads$TSS_index == 3)) | 
                                     ((Mismatched_reads$Mismatch_type == "2M3M") & (Mismatched_reads$TSS_index == 3)) | 
                                     ((Mismatched_reads$Mismatch_type == "-1M2M") & (Mismatched_reads$TSS_index == 2)) | 
                                     ((Mismatched_reads$Mismatch_type == "-1M3M") & (Mismatched_reads$TSS_index == 3)) | 
                                     ((Mismatched_reads$Mismatch_type == "-2M3M") & (Mismatched_reads$TSS_index == 3))),]
####################################################################################################################
TSS_slippage_plus = TSS_slippage[(TSS_slippage$TSS_strand == "+"),]
TSS_slippage_minus = TSS_slippage[(TSS_slippage$TSS_strand == "-"),]
rownames(TSS_slippage_plus) = NULL
rownames(TSS_slippage_minus) = NULL

Cut_reads = Cut_reads[(Cut_reads$TSS_index <= 6),]
Cut_reads_plus = Cut_reads[Cut_reads$Strand == "+",]
Cut_reads_minus = Cut_reads[Cut_reads$Strand == "-",]
rownames(Cut_reads_plus) = NULL
rownames(Cut_reads_minus) = NULL

Mismatched_reads = Mismatched_reads[(Mismatched_reads$TSS_index <= 6),]
Mismatched_reads = Mismatched_reads[, 1:8]
Mismatched_reads_plus = Mismatched_reads[(Mismatched_reads$Strand == "+"),]
Mismatched_reads_minus = Mismatched_reads[(Mismatched_reads$Strand == "-"),]
##########
#Find the nearest TSS for each reads
for(i in 1:nrow(Cut_reads_plus))
{
  Cut_reads_plus[i,9] = TSS_slippage_plus[which.min(abs(TSS_slippage_plus$TSS_coordinate - Cut_reads_plus[i,7])),1]
  Cut_reads_plus[i,10] = TSS_slippage_plus[which.min(abs(TSS_slippage_plus$TSS_coordinate - Cut_reads_plus[i,7])),2]
}

colnames(Cut_reads_plus)[9:10] = c("Near_TSS","TSS_group")
cut_match_plus = Cut_reads_plus[(Cut_reads_plus$TSS_coor == Cut_reads_plus$Near_TSS),]

for(i in 1:nrow(Cut_reads_minus))
{
  Cut_reads_minus[i,9] = TSS_slippage_minus[which.min(abs(TSS_slippage_minus$TSS_coordinate + Cut_reads_minus[i,7])),1]
  Cut_reads_minus[i,10] = TSS_slippage_minus[which.min(abs(TSS_slippage_minus$TSS_coordinate + Cut_reads_minus[i,7])),2]
}

Cut_reads_minus[,9] = - Cut_reads_minus[,9]
colnames(Cut_reads_minus)[9:10] = c("Near_TSS","TSS_group")
cut_match_minus = Cut_reads_minus[(Cut_reads_minus$TSS_coor == Cut_reads_minus$Near_TSS),]
####################################################################################################################
#Adjust the slippage ratio due to the reads aligned to the genome after 5' end trimming
for(i in 1:nrow(TSS_slippage))
{
  if(TSS_slippage[i,2] == "+")
  {
    if(any(TSS_slippage[i,1] == cut_match_plus[,]$TSS_coor))
    {
      TSS_slippage[i,6] = sum((TSS_slippage[i,1] == cut_match_plus[,]$TSS_coor), na.rm = TRUE)
    }
    else
    {
      TSS_slippage[i,6] = 0
    }
  }
  else
  {
    if(any(TSS_slippage[i,1] == -cut_match_minus[,]$TSS_coor))
    {
      TSS_slippage[i,6] = sum((TSS_slippage[i,1] == -cut_match_minus[,]$TSS_coor), na.rm = TRUE)
    }
    else
    {
      TSS_slippage[i,6] = 0
    }
  }
}

colnames(TSS_slippage)[6] = "5cut_align_num"
TSS_slippage[,7] = TSS_slippage$TSS_RD + TSS_slippage$`5cut_align_num`
TSS_slippage[,8] = TSS_slippage$TSS_mm_RD + TSS_slippage$`5cut_align_num`
TSS_slippage[,9] = TSS_slippage[,8]/TSS_slippage[,7]*100
colnames(TSS_slippage)[7:9] = c("adj_TSS_RD","adj_TSS_mm_RD","adj_TSS_mm_ratio")

TSS_slippage = TSS_slippage[,c(1,2,7:9)]
####################################################################################################################
write.table(TSS_slippage,  file="./1.0.2_TSS_slippage_Ratio_adjust.csv.csv",  sep = ",", col.names = TRUE, row.names = FALSE, quote = FALSE)
