library(ggplot2)
library(tidyverse)
library(ggbreak) 

core <- read.csv('~/RumenCampylobacter2022/Figures/Fig1_AB/intermediate_tables/all_genomes_core_nucl.csv') %>%
  select(gene1) %>%
  mutate(gene1 = as.character(gene1)) %>%
  separate(gene1, into=c("label", "gene_num"), sep = "_", remove = FALSE) %>%
  mutate(gene_num = as.numeric(gene_num)) %>%
  mutate(class = "core")

pop_specific <- read.csv('~/RumenCampylobacter2022/Figures/Fig1_AB/intermediate_tables/all_genomes_pop_specific_nucl.csv') %>%
  select(gene1) %>%
  separate(gene1, into=c("label", "gene_num"), sep = "_", remove = FALSE) %>%
  mutate(gene_num = as.numeric(gene_num)) %>%
  mutate(class = "pop_specific")

flex <- read.csv('~/RumenCampylobacter2022/Figures/Fig1_AB/intermediate_tables/all_genomes_flex_nucl.csv') %>%
  select(gene1) %>%
  mutate(gene1 = as.character(gene1)) %>%
  separate(gene1, into=c("label", "gene_num"), sep = "_", remove = FALSE) %>%
  mutate(gene_num = as.numeric(gene_num)) %>%
  mutate(class = "flex")

gene_categories <- bind_rows(core, pop_specific) %>%
  bind_rows(flex) %>%
  select(-gene1)

gene_categoriesX <- gene_categories %>% 
  filter(label == "KHAFFMAE") 
colnames(gene_categoriesX) <- c("label1", "gene1_num", "classX")

gene_categoriesY <- gene_categories %>% 
  filter(label == "HKMHBKFO") %>%
  rename(label2 = label) %>%
  rename(gene2_num = gene_num) %>%
  rename(classY = class) %>%
  mutate(gene2_num_adjust = if_else(gene2_num <= 1101, gene2_num+462, gene2_num-1101)) %>%
  mutate(gene2_num_adjust = abs(gene2_num_adjust - 1564)) %>%
  select(-gene2_num)

df_F0 <- read_delim("~/RumenCampylobacter2022/Processing/genomes/output/JMF18_spades_pomoxis_polished_min2000_V_131980_spades_pomoxis_polished_min2000.ffn.txt", delim = "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
df_F0$file <- "JMF18_spades_pomoxis_polished_min2000_V_131980_spades_pomoxis_polished_min2000"
df_F1 <- read_delim("~/RumenCampylobacter2022/Processing/genomes/output/131980_spades_pomoxis_polished_min2000_V_JMF18_spades_pomoxis_polished_min2000.ffn.txt", delim = "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
df_F1$file <- "131980_spades_pomoxis_polished_min2000_V_JMF18_spades_pomoxis_polished_min2000"

compiled_complete_blast_hits_nucl <- bind_rows(df_F0, df_F1)
colnames(compiled_complete_blast_hits_nucl )[1:14] <- c("gene1", "gene2", "pident", "sstart", "send", "qstart", "qend", "evalue", "bitscore", "score", "qlen", "length", "sseq", "file")

compiled_complete_blast_hits_nucl  <- compiled_complete_blast_hits_nucl  %>%
  
  select(-sseq) %>%
  separate(file, into = c("genome1", "genome2"), sep = "_V_") %>%
  mutate(palign = (length / qlen)*100) %>%
  
  group_by(gene1) %>%
  top_n(1, bitscore) %>%
  ungroup() %>%
  
  filter(palign > 70) %>%
  filter(pident > 70) 

compiled_complete_blast_hits_nucl1 <- select(compiled_complete_blast_hits_nucl, genome1, gene1, genome2, gene2, pident, palign)
colnames(compiled_complete_blast_hits_nucl1)[5] <- "pident1"
colnames(compiled_complete_blast_hits_nucl1)[6] <- "palign1"

compiled_complete_blast_hits_nucl2 <- select(compiled_complete_blast_hits_nucl, genome1, gene1, genome2, gene2, pident, palign)
colnames(compiled_complete_blast_hits_nucl2) <- c("genome2", "gene2", "genome1", "gene1", "pident2", "palign2")

reciprocal_best_hits_nucl <- inner_join(compiled_complete_blast_hits_nucl1, compiled_complete_blast_hits_nucl2) %>%
  mutate(pident = (pident1 + pident2)/2) %>%
  mutate(palign = (palign1 + palign2)/2) %>%
  select(-pident1, -pident2, -palign1, -palign2)

df_plot <- reciprocal_best_hits_nucl %>%
  
  filter(genome1 == "JMF18_spades_pomoxis_polished_min2000") %>%
  filter(genome2 == "131980_spades_pomoxis_polished_min2000") %>%
  
  select(gene1, gene2, pident) %>%
  mutate(gene1 = as.character(gene1)) %>%
  mutate(gene2 = as.character(gene2)) %>%
  
  separate(gene1, into=c("label1", "gene1_num"), sep = "_", remove = FALSE) %>%
  separate(gene2, into=c("label2", "gene2_num"), sep = "_", remove = FALSE) %>%
  
  mutate(gene1_num = as.numeric(gene1_num)) %>%
  mutate(gene2_num = as.numeric(gene2_num)) %>%
  
  mutate(gene2_num_adjust = if_else(gene2_num <= 1101, gene2_num+462, gene2_num-1101)) %>%
  mutate(gene2_num_adjust = abs(gene2_num_adjust - 1564)) %>%
  
  select(label1, gene1_num, label2, gene2_num_adjust, pident) %>%
  
  full_join(gene_categoriesX) %>%
  mutate(rm = if_else(is.na(pident) & classX == "core", "yes", "no")) %>%
  filter(rm != "yes") %>%
  
  full_join(gene_categoriesY) %>%
  mutate(rm = if_else(is.na(pident) & classY == "core", "yes", "no")) %>%
  mutate(rm = if_else(!(is.na(classX)), "no", rm)) %>%
  filter(rm != "yes")  %>%
  
  select(-rm) %>%
  
  mutate(label1 = if_else(is.na(label1), "KHAFFMAE", label1)) %>%
  mutate(gene1_num = if_else(is.na(gene1_num), 0, gene1_num)) %>%
  
  mutate(label2 = if_else(is.na(label2), "HKMHBKFO", label2)) %>%
  mutate(gene2_num_adjust = if_else(is.na(gene2_num_adjust), 0, gene2_num_adjust)) %>%
  
  mutate(pident = if_else(is.na(pident), 0, pident)) %>%
  
  mutate(classX = if_else(is.na(classX), "none", classX)) %>%
  mutate(classY = if_else(is.na(classY), "none", classY)) %>%
  
  mutate(flexX = if_else(classX == "flex", "0", "NA")) %>%
  mutate(flexY  = if_else(classY == "flex", "0", "NA")) %>%
  mutate(pop_specificX = if_else(classX == "pop_specific", "0", "NA")) %>%
  mutate(pop_specificY = if_else(classY == "pop_specific", "0", "NA"))
  
df_plot$flexX <- as.numeric(df_plot$flexX)
df_plot$flexY <- as.numeric(df_plot$flexY)
df_plot$pop_specificX <- as.numeric(df_plot$pop_specificX)
df_plot$pop_specificY <- as.numeric(df_plot$pop_specificY)

### PLOT

pdf("~/RumenCampylobacter2022/Figures/Fig1_AB/output/fig1_A.pdf", width=9, height=4.5)

ggplot(df_plot, aes(x=gene1_num, y=gene2_num_adjust)) +
  geom_point(aes(colour=pident), shape = 108, size = 2) + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) +
  scale_colour_gradient(low = "white", high = "black", na.value = NA, limits = c(75,100)) + 
  geom_point(aes(x=gene1_num, y=flexX), colour = "red", shape = 108, size = 2) +
  geom_point(aes(x=gene1_num, y=pop_specificX ), colour = "blue", shape = 108, size = 2) +
  geom_point(aes(x=flexY, y=gene2_num_adjust), colour = "red", shape = 108, size = 2) +
  geom_point(aes(x=pop_specificY, y=gene2_num_adjust), colour = "blue", shape = 108, size = 2) 
  
dev.off()

###

gene_categories_summary <- df_plot %>%
  
  select(classX, classY) %>%
  
  gather(genome, class) %>%
  
  group_by(genome, class) %>%
  mutate(total_class = n()) %>%
  ungroup() %>%
  
  distinct() %>%
  
  filter(class != "none") %>%
  
  group_by(genome) %>%
  mutate(total_genome = sum(total_class)) %>%
  ungroup() %>%
  
  mutate(per = (total_class / total_genome)*100)

gene_categories_summary$total_class <- as.numeric(gene_categories_summary$total_class)

scaleFUN <- function(x) sprintf("%.0f", x)

### PLOT

pdf("~/RumenCampylobacter2022/Figures/Fig1_AB/output/fig1_B.pdf", width=6, height=6)

ggplot(gene_categories_summary, aes(x=class, y=total_class, fill=genome)) + 
  
  geom_bar(stat="identity", position=position_dodge()) +
  
  theme_minimal() +
  
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))   

dev.off()

###

pop_specific_annotations <- read.csv("~/RumenCampylobacter2022/Processing/genomes/output/compiled_annotations.csv") %>%
  filter(locus_tag %in% pop_specific$gene1) 

write.csv(pop_specific_annotations, "~/RumenCampylobacter2022/Figures/Fig1_AB/output/pop_specific_annotations_total.csv")

pop_specific_annotations_noHyp <- pop_specific_annotations %>%
  filter(!(product == "hypothetical protein" & product2 == "hypothetical protein" & product3 == "hypothetical protein"))

write.csv(pop_specific_annotations_noHyp, "~/RumenCampylobacter2022/Figures/Fig1_AB/output/pop_specific_annotations_no_hypothetical.csv")

pop_specific_annotations$product3 <- gsub("putative protein", "hypothetical protein", pop_specific_annotations$product3)

annotation_count_summary <- pop_specific_annotations %>%
  select(file, product3) %>%
  
  group_by(file, product3) %>%
  mutate(count = n()) %>%
  ungroup() %>%
  
  group_by(file) %>%
  mutate(total = n()) %>%
  ungroup() %>%
  
  distinct() %>%
  
  mutate(per = (count / total)*100)

write.csv(annotation_count_summary, "~/RumenCampylobacter2022/Figures/Fig1_AB/output/annotation_count_summary.csv")
