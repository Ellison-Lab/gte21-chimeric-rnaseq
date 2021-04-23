library(tidyverse)
library(magrittr)

set.seed(1)

#chim_reads_fl <- Sys.glob('~/amarel-scratch/TE-proj-reorg/gte21-chimeric-rnaseq/results/star/larval_testes_cleaned_papain_01/Chimeric.out.junction')

chim_reads_fl <- snakemake@params[["junctions"]]
samp <- snakemake@wildcards[["s"]] # samp <-"larval_testes_cleaned_papain_01"

main_chroms <- c("2L", "2R","3L","3R","4","X","Y")

# import chimeric junctions file
chimeric_reads_0 <-  chim_reads_fl %>%
  read_tsv(.,comment = '#')

# filter reads that map to a single breakpoint on a main chrom and a te and have the highest aln score of all mappings
chimeric_reads_1 <- chimeric_reads_0 %>%
  filter(num_chim_aln == 1) %>%
  group_by(read_name) %>%
  filter(length(unique(chr_acceptorB)) == 1) %>%
  filter(length(unique(chr_donorA)) == 1) %>%
  filter(this_chim_aln_score == bestall_chim_aln_score) %>%
  filter(n() == 1) %>%
  ungroup() %>%
  #filter(xor(chr_acceptorB %in% te_lookup$gene_id, chr_donorA %in% te_lookup$gene_id)) %>%
  filter(xor(chr_donorA %in% main_chroms, chr_acceptorB %in% main_chroms)) %>%
  filter(chr_donorA!=chr_acceptorB)

# annotate
chimeric_reads_2 <- chimeric_reads_1 %>%
  mutate(chr = ifelse(chr_donorA %in% main_chroms,chr_donorA, chr_acceptorB)) %>%
  mutate(start = ifelse(chr_donorA %in% main_chroms,brkpt_donorA, brkpt_acceptorB)) %>%
  mutate(strand = ifelse(chr_donorA %in% main_chroms,strand_donorA, strand_acceptorB)) %>%
  mutate(other_chr = ifelse(!chr_donorA %in% main_chroms,chr_donorA, chr_acceptorB)) %>%
  mutate(other_brk = ifelse(!chr_donorA %in% main_chroms,brkpt_donorA, brkpt_acceptorB)) %>%
  mutate(other_strand = ifelse(!chr_donorA %in% main_chroms,strand_donorA, strand_acceptorB)) %>%
  arrange(chr_donorA, chr_acceptorB, brkpt_donorA, brkpt_acceptorB) %>%
  mutate(r1 = chr_donorA)

breakpt_df <- chimeric_reads_2 %>%
  dplyr::select(chr, start, strand , other_chr, other_brk, other_strand, read_name, r1) %>%
  group_by(chr, start, strand,  other_chr, other_brk, other_strand,  r1) %>%
  summarise(supporting_reads=length(unique(read_name)), .groups = "drop") %>%
  mutate(breakpoint_id = paste(chr, start, strand, other_strand, other_chr,other_brk, sep="_")) %>%
  mutate(sample = samp)

write_csv(breakpt_df, snakemake@output[["csv"]])
