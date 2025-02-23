library(tidyverse)
library(viridis)

pops <- read.csv("~/RumenCampylobacter2022/METAdata/population_classification.csv")
colnames(pops)[1] <- "genome2"

campy_genomes <- pops$genome2

# reduce the output table size of the all against all blast results so that I can upload it to github.

# reciprocal_best_hits_nucl <- read_csv("~/RumenCampylobacter2022/Processing/genomes/output/reciprocal_best_hits_nucl.csv", col_types = cols(.default = "c")) %>% 
#   filter(genome1 %in% campy_genomes) %>% 
#   filter(genome2 %in% campy_genomes)  %>%
#   filter(genome1 == "131980_spades_pomoxis_polished_min2000") %>%
#   filter(genome2 == "JMF18_spades_pomoxis_polished_min2000") %>%
#   filter(genome1 != genome2)
# 
# write.csv(reciprocal_best_hits_nucl, "~/RumenCampylobacter2022/Processing/genomes/output/reciprocal_best_hits_nucl_red.csv", row.names = FALSE)

reciprocal_best_hits_nucl <- read_csv("~/RumenCampylobacter2022/Processing/genomes/output/reciprocal_best_hits_nucl_red.csv", col_types = cols(.default = "c"))

reciprocal_best_hits_nucl$pident_nucl <- as.numeric(reciprocal_best_hits_nucl$pident_nucl)
reciprocal_best_hits_nucl$palign_nucl <- as.numeric(reciprocal_best_hits_nucl$palign_nucl)

rep_rbh_nucl <- reciprocal_best_hits_nucl %>%
  filter(palign_nucl >= 80) 



hist(rep_rbh_nucl$pident_nucl, breaks = 100)

rep_rbh_nucl <- rep_rbh_nucl %>%
  filter(pident_nucl <= 97.5)

rep_rbh_nucl$gene_num_rbh <- seq(1, nrow(rep_rbh_nucl), 1)

rep_rbh_nucl_num <- rep_rbh_nucl %>%
  select(gene1, gene2, gene_num_rbh) %>%
  gather(rm, ID, -gene_num_rbh) %>%
  select(-rm)

transcriptome_map <- read_csv("~/RumenCampylobacter2022/METAdata/sample_mapping.csv", col_types = cols(.default = "c"))

expression_counts <- read.csv("~/RumenCampylobacter2022/Processing/metatranscriptomes/output/compiled_counts_clones.csv")  %>%
  
  select(-gene, -product) %>%
  separate(reads, into = c("sample", "type"), sep = "\\.R1\\.") %>%
  
  filter(type != "metagenome") %>%
  
  inner_join(transcriptome_map) %>%
  
  inner_join(rep_rbh_nucl_num) %>%
  
  separate(ID, into = c("genome"), sep = "_") %>%
  
  unite(genome_sample, c("genome", "sample"), sep = "_") %>%
  
  select(-type, -treatment) %>%
  
  spread(genome_sample, count)

library(DESeq2)

cts <- as.matrix(expression_counts[,2:13])
rownames(cts) <- expression_counts$gene_num_rbh

coldata = matrix(c("131980_spades_pomoxis_polished_min2000", "131980_spades_pomoxis_polished_min2000", "131980_spades_pomoxis_polished_min2000", "131980_spades_pomoxis_polished_min2000", "131980_spades_pomoxis_polished_min2000", "131980_spades_pomoxis_polished_min2000", "JMF18_spades_pomoxis_polished_min2000", "JMF18_spades_pomoxis_polished_min2000", "JMF18_spades_pomoxis_polished_min2000" , "JMF18_spades_pomoxis_polished_min2000", "JMF18_spades_pomoxis_polished_min2000", "JMF18_spades_pomoxis_polished_min2000"))
colnames(coldata) <- "genome"
rownames(coldata) <- c(names(expression_counts)[2:13])

dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design= ~ genome)


dds <- DESeq(dds)
resultsNames(dds)
res <- results(dds)
deseq <- as.data.frame(res, name="genome_JMF_2102_8_0018_vs_131980_min1000")
rm(list = c("res", "dds", "cts", "coldata"))

deseq$gene_num_rbh <- rownames(deseq)

deseq <- deseq %>%
  select(gene_num_rbh, pvalue, padj, log2FoldChange) %>%
  filter(!(is.na(padj)))

rep_rbh_nucl_num$gene_num_rbh <- as.character(rep_rbh_nucl_num$gene_num_rbh)

summary <- deseq %>%
  inner_join(rep_rbh_nucl_num) %>%
  mutate(neglog10pAdj = -log10(padj)) %>%
  mutate(direction = if_else(neglog10pAdj > 5 & log2FoldChange >= 1.5, "S39", 
                             if_else(neglog10pAdj > 5 & log2FoldChange <= -1.5, "22B", "stable")))

pdf("~/RumenCampylobacter2022/Figures/Fig4_A/output/fig4_A.pdf", width=9, height=4.5)

ggplot(data = summary, 
       aes(x = log2FoldChange, 
           y = neglog10pAdj, 
           colour = direction)) +
  
  scale_color_manual(values=c("blue", "red", "grey")) +
  
  geom_point(alpha=0.4, size=3.5) +
  
  geom_vline(xintercept=c(-1.5,1.5),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = 5,lty=4,col="black",lwd=0.8) +
  
  labs(x="log2(fold change)",
       y="-log10 (adj.p-value)",
       title="Differential expression")  +
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="right", 
        legend.title = element_blank())

dev.off()

annotations <- read.csv("~/RumenCampylobacter2022/Processing/genomes/output/compiled_annotations.csv") 
colnames(annotations)[2] <- "ID"

summary_annotations <- summary %>%
  filter(direction != "stable") %>%
  inner_join(annotations) 