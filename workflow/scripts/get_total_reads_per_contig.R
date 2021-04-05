library(Rsamtools)
library(tidyverse)

#bam <- Sys.glob('~/amarel-scratch/TE-proj-reorg/gte21-chimeric-rnaseq/results/srt/larval_testes_cleaned_papain_01.srt.bam')

bam <- snakemake@input[["bam"]]

idxstats <- Rsamtools::idxstatsBam(bam) %>% 
  as_tibble() %>% 
  dplyr::select(-unmapped, other_chr.seqlength  = seqlength, other_chr.mapped = mapped)

write_csv(idxstats, snakemake@output[["csv"]])
