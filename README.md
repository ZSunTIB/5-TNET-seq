# 5-TNET-seq
Directions and custom scripts to process the 5'TNET-seq data. The raw fastq.gz files obtained on an NovaSeq 6000 platform (2 × 50 bp paired end) were processed to find the pause sites, extract the reads derived from productive reiterative initatioin and calculate the slippage ratio for each TSS region.

Following are versions of the operating system and packages used:
macOS Monterey 12.0.1, bbmap 38.22, Flash 1.2.11, bowtie 1.2.2, STAR 2.6.1, bedtools 2.27.1, R 4.4.0

Following is the procedure for processing the 5'TENT-seq data:

1. Run 0.0.1 5'TNET-seq_reads_aligner.sh to align the sequencing reads to E. coli genome;
2. Run 0.0.2 5'TNET-seq_pause site_finder.R to find the RNAP pause sites;
3. Run 0.1.1 5'TNET-seq_unmapped reads alignment.sh to extract all reads that were unable to align to the genome;
4. Run 0.1.2 5'TNET-seq_unmapped reads alignment_cycle.sh to trim the 5' end nucleotides from the unmapped reads and align the residual sequences to the genome;
5. Run 0.1.3 5'TNET-seq_1to3_mismatch_reads_extraction.sh to extract the reads that containing mismatches in the 1st, 2nd, or 3rd coordinates;
6. Run 0.1.4 5'TNET-seq_5'trimmed_1to3_mismatch_reads_extraction.sh to extract the reads that containing mismatches in the 1st, 2nd, or 3rd coordinates after 5' end trimming;
7. Run 0.2.1_TSS_mm_Ratio_calculate.R to calculate the slippage ratio for each TSS region;
8. Run 0.2.2_TSS_mm_Ratio_adjust.R to adjust the slippage ratio for each TSS region.

The expected run time depends on the sequencing depth.
