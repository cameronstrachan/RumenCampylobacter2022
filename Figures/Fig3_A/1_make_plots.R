library(tidyverse)

allele_count <- read_table("~/RumenCampylobacter2022/Processing/metagenomes/output/R1_clones_filtered_onlySNPs.frq.count", col_names = c("contig", "position", "n_alleles", "total_count", "allele_1", "allele_2", "allele_3", "allele_4"), skip = 1)

allele_count_gather <- allele_count %>%
  
  gather(allele, count, -contig, -position, -n_alleles, -total_count) %>%
  
  filter(!is.na(count)) %>%
  select(-allele) %>%
  separate(count, into = c("base", "count"), sep = ":") %>%
  
  filter(count != "0") %>%
  
  select(-n_alleles, -total_count, -count) 
  
  
S0_mutations <- allele_count_gather %>%
  filter(contig == "NODE_1_length_1570312_cov_317.865870") %>%
  arrange(position) %>%
  mutate(diff = position - lag(position, default = first(position)))
  
S1_mutations <- allele_count_gather %>%
  filter(contig == "NODE_1_length_1536821_cov_689.307783") %>%
  arrange(position) %>%
  mutate(diff = position - lag(position, default = first(position)))



window = 1000

start = 428381
end = 628381
len = end - start

region_grouping <- as.data.frame(cbind(seq(start, end - 1, 1), rep(1:(len/window), each=window)))
colnames(region_grouping) <- c("position", "wind")

S0_mutations_region <- S0_mutations %>%
  filter(position >= start) %>%
  filter(position <= end) %>%
  select(position) %>%
  inner_join(region_grouping) %>%
  group_by(wind) %>%
  mutate(n_mutations = n()) %>%
  ungroup() %>%
  select(-position) %>%
  distinct() %>%
  full_join(region_grouping) %>%
  mutate(n_mutations = replace_na(n_mutations, 0)) %>%
  select(-position) %>%
  mutate(position = (wind*window) + start) %>%
  distinct() %>%
  mutate(percent_snps = (n_mutations / window) *100)

library(zoo)

ggplot(S0_mutations_region, aes(x = position, y = percent_snps)) + 
  geom_point() + 
  geom_line(aes(y=rollmean(percent_snps, 200))) + 
  theme_bw() 

597184

window = 1000

start = 497184
end = 697184
len = end - start

region_grouping <- as.data.frame(cbind(seq(start, end - 1, 1), rep(1:(len/window), each=window)))
colnames(region_grouping) <- c("position", "wind")

S1_mutations_region <- S1_mutations %>%
  filter(position >= start) %>%
  filter(position <= end) %>%
  select(position) %>%
  inner_join(region_grouping) %>%
  group_by(wind) %>%
  mutate(n_mutations = n()) %>%
  ungroup() %>%
  select(-position) %>%
  distinct() %>%
  full_join(region_grouping) %>%
  mutate(n_mutations = replace_na(n_mutations, 0)) %>%
  select(-position) %>%
  mutate(position = (wind*window) + start) %>%
  distinct() %>%
  mutate(percent_snps = (n_mutations / window) *100)

pdf("~/RumenCampylobacter2022/Figures/Fig3_A/output/fig3_A.pdf", width=12, height=4.5)

ggplot(S1_mutations_region, aes(x = position, y = percent_snps)) + 
  geom_point() + 
  #geom_line(aes(y=rollmean(percent_snps, 200))) + 
  theme_bw() 

dev.off()
