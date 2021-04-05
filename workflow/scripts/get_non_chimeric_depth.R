library(tidyverse)
library(magrittr)
library(VariantAnnotation)
library(Rsamtools)

set.seed(1)

#bam <- Sys.glob('~/amarel-scratch/TE-proj-reorg/gte21-chimeric-rnaseq/results/srt/larval_testes_cleaned_papain_01.srt.bam')

#brk_fl <- "~/amarel-scratch/TE-proj-reorg/gte21-chimeric-rnaseq/results/breakpoints/larval_testes_cleaned_papain_01.breakpoints.csv"

bam <- snakemake@input[["bam"]]

brk_fl <- snakemake@input[["csv"]]

brks <- read_csv(brk_fl)

brks2 <- brks %>% dplyr::select(chr = other_chr, start = other_brk) %>% distinct() %>% mutate(end = start)

gr <- GRanges(brks2) %>% unstrand() %>% GenomicRanges::reduce()

pup <- PileupParam(max_depth = 10e6,
                   distinguish_nucleotides = F,
                   distinguish_strands = F,min_mapq = 0, 
                   min_minor_allele_depth = -1,
                   min_nucleotide_depth = -1)

scbp <- ScanBamParam(which = gr)

pileups.df <- pileup(bam, scanBamParam = scbp, pileupParam = pup) %>%
  as_tibble(.) 
  
res <- brks %>% left_join(pileups.df, by=c(other_chr = "seqnames",other_brk="pos")) %>% 
  dplyr::rename(other_reads = count) %>% 
  dplyr::select(-which_label) %>%
  replace_na(list(other_reads=0)) %>%
  dplyr::relocate(other_reads, .after = supporting_reads)

write_csv(res, snakemake@output[["csv"]])
