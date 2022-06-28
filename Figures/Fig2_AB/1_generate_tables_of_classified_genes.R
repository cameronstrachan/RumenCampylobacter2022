library(tidyverse)

df <- read_csv("~/RumenCampylobacter2022/Processing/genomes/output/genes_against_genomes.csv")

complete_genomes <- c("JMF18_spades_pomoxis_polished_min2000", "131980_spades_pomoxis_polished_min2000")

S0 <- c("131980_spades_pomoxis_polished_min2000", "131981_min1000", "131982_min1000", "JMF_2102_8_0021", "JMF_2102_8_0023_filtered")

S1 <- c("JMF18_spades_pomoxis_polished_min2000", "JMF_2102_8_0019_filtered", "JMF_2102_8_0020_filtered", "JMF_2102_8_0024_filtered", 
        "JMF_2102_8_0013_filtered", "JMF_2102_8_0016", "JMF_2102_8_0017_filtered", "JMF_2102_8_0012_filtered")

all_genomes <- c(S0, S1)

df_hit_counts <- df %>% 
  
  filter(palign > 70) %>%
  filter(pident > 70) %>%
  
  filter(genome1 %in% complete_genomes) %>%
  filter(genome2 %in% all_genomes) %>%
  filter(genome1 != genome2) %>%
  
  group_by(genome1) %>%
  mutate(n_genes = length(unique(gene1))) %>%
  ungroup() %>%
  
  mutate(pop1 = if_else(genome1 == "131980_spades_pomoxis_polished_min2000", "F0", 
                if_else(genome1 == "JMF18_spades_pomoxis_polished_min2000", "F1", "NA"))) %>%
  
  mutate(pop2 = if_else(genome2 %in% S0, "F0", 
                if_else(genome2 %in% S1, "F1", "NA"))) %>%
  
  mutate(comparison = if_else(pop1 == pop2, "within", "between")) %>%
  
  group_by(genome1, gene1, genome2) %>%
  
  group_by(genome1, gene1) %>%
  mutate(n_genome_hits = length(unique(genome2))) %>% 
  ungroup() %>%
  
  group_by(genome1, gene1, comparison) %>%
  mutate(n_genome_hits_comp = length(unique(genome2))) %>%
  ungroup() %>%
  
  group_by(genome1, comparison) %>%
  mutate(n_genomes_comp = length(unique(genome2))) %>%
  ungroup() %>%
  
  select(genome1, gene1, comparison, n_genes, n_genome_hits, n_genome_hits_comp, n_genomes_comp) %>%
  
  distinct() %>%
  
  mutate(per_comp = (n_genome_hits_comp/n_genomes_comp)*100) 

####

df_core <- df_hit_counts %>%
  
  filter(n_genome_hits == 12) %>%
  
  group_by(genome1) %>%
  mutate(n_core_genes = length(unique(gene1))) %>%
  ungroup() %>%
  
  select(genome1, gene1, n_genes, n_core_genes) %>%
  
  distinct() %>%
  
  mutate(per_core_genes = (n_core_genes / n_genes)*100) 

write.csv(df_core, '~/RumenCampylobacter2022/Figures/Fig2_AB/intermediate_tables/all_genomes_core_nucl.csv')

df_pop_specific <- df_hit_counts %>%
  
  filter(n_genome_hits == n_genome_hits_comp) %>%
  
  filter(per_comp == 100) %>%
  
  group_by(genome1) %>%
  mutate(n_pop_specific_genes = length(unique(gene1))) %>%
  ungroup() %>%
  
  select(genome1, gene1, n_genes, n_pop_specific_genes) %>%
  
  distinct() %>%
  
  mutate(per_pop_specific_genes = (n_pop_specific_genes / n_genes)*100) 

write.csv(df_pop_specific, '~/RumenCampylobacter2022/Figures/Fig2_AB/intermediate_tables/all_genomes_pop_specific_nucl.csv')

df_flex <- df_hit_counts %>%
  
  filter(n_genome_hits < 11) %>%
  
  filter(!(gene1 %in% unique(df_pop_specific$gene1))) %>%
  
  group_by(genome1) %>%
  mutate(n_flex_genes = length(unique(gene1))) %>%
  ungroup() %>%
  
  select(genome1, gene1, n_genes, n_flex_genes) %>%
  
  distinct() %>%
  
  mutate(per_flex_genes = (n_flex_genes / n_genes)*100)

write.csv(df_flex, '~/RumenCampylobacter2022/Figures/Fig2_AB/intermediate_tables/all_genomes_flex_nucl.csv')