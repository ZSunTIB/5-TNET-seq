###############################
library("data.table")
library("RcppRoll")
library("stringr")
library("dplyr")
library("MALDIquant")
library("Biostrings")
library("tidyverse")
####################################################################################################################
#Read raw data of 3' end read count for each coordinate
wig_plus = fread("./Processed data/Processed data_uniq_3end_Positive.wig")
wig_minus = fread("./Processed data/Processed data_uniq_3end_Positive.wig")

total_reads = sum(wig_plus$V3 + wig_minus$V3)
####################################################################################################################
#Set the variable cutoffs to pick up Pause sites.
#Export the Pause sites satisfying the following 2 conditions:
# 1. Read count of all merged reads 3' end >= 10 reads/million;
count_M = 10
# 2. Score >= 20 (Score = 3' end count of merged reads/median 3' end count of merged reads in a surrounding 51-bp window size);
score_cut = 20
####################################################################################################################
ps_plus = wig_plus[, .(V2,V3)]
setnames(ps_plus,c("coordinate","count"))
##########################
setkey(ps_plus, coordinate)

#Calculate the rolling median in a 51-nt window size surrounding each coordinate.
ps_plus_median = ps_plus[ , Roll_median51 := roll_median(count, n = 51, weights = NULL, by = 1L, fill = NA, partial = FALSE,
                                      align = c("center", "left", "right"), normalize = TRUE, na.rm = FALSE),]

#Calculate the Score of each coordinate (Count/Roll_median51)
ps_plus_median_divide = ps_plus_median[, score := count / Roll_median51]
ps_plus_median_divide_C = ps_plus_median_divide[, score_F := score]
ps_plus_F = ps_plus_median_divide_C[score %in% c("Inf", "NaN"), score_F := count]
ps_plus_Final = ps_plus_median_divide_C[, Strand := "+"]

#Export the pause sites
pause_sites_plus = ps_plus_Final[(count >= (total_reads/1000000)*count_M) & (score_F >= score_cut)]
############################################################################################################################
ps_minus = wig_minus[, .(V2,V3)]
setnames(ps_minus,c("coordinate","count"))
##########################
setkey(ps_minus, coordinate)

#Calculate the Score of each coordinate (Count/Roll_median51)
ps_minus_median = ps_minus[ , Roll_median51 := roll_median(count, n = 51, weights = NULL, by = 1L, fill = NA, partial = FALSE,
                                        align = c("center", "left", "right"), normalize = TRUE, na.rm = FALSE),]

#Calculate the Score of each coordiante (Count/Roll_median51)
ps_minus_median_divide = ps_minus_median[, score := count / Roll_median51]
ps_minus_median_divide_C = ps_minus_median_divide[, score_F := score]
ps_minus_F = ps_minus_median_divide_C[score %in% c("Inf", "NaN"), score_F := count]
ps_minus_Final = ps_minus_median_divide_C[, Strand := "-"]

#Export the pause sites
pause_sites_minus = ps_minus_Final[(count >= (total_reads/1000000)*count_M) & (score_F >= score_cut)]
############################################################################################################################
pause_sites = rbind(pause_sites_plus, pause_sites_minus)
pause_sites = pause_sites[order(-score_F)]
############################################################################################################################
#Remove the pause sites in the tRNA and rRNA regions.
trRNA = read.delim("./trRNA_region_NC_000913.2.txt", header = FALSE, sep = "\t", fill = TRUE, row.names = NULL)
colnames(trRNA) = c("Genome","Start","End")
pause_sites = data.frame(pause_sites)
pause_sites[,(ncol(pause_sites)+1)] = 0

for (i in 1:nrow(pause_sites))
{
  if (any(pause_sites[i,1]>=trRNA[,2] & pause_sites[i,1]<=trRNA[,3]))
  {
    pause_sites[i,ncol(pause_sites)] = pause_sites[i,ncol(pause_sites)] + 1;
  }
  else
  {
  }
}
pause_sites_rm_trRNA = pause_sites[(pause_sites[,ncol(pause_sites)]==0),]


pause_sites_rm_trRNA = pause_sites_rm_trRNA[,c(1:3,5,6)]
rownames(pause_sites_rm_trRNA) = NULL
############################################################################################################################
write.table(pause_sites_rm_trRNA, file="./Processed data/Processed data_PS.csv",  sep = ",", col.names = TRUE, row.names = FALSE, quote = FALSE)
############################################################################################################################
rm(ps_minus,ps_minus_F,ps_minus_median,ps_minus_median_divide,ps_minus_median_divide_C)
rm(ps_plus,ps_plus_F,ps_plus_median,ps_plus_median_divide,ps_plus_median_divide_C)
