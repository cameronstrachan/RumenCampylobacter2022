library(tidyverse)
library(rRDP)

counts <- read.delim("output/asv-counts-merged.txt", header=FALSE)
counts <- counts[-1,]
counts[] <- lapply(counts, as.character)
names(counts) <- counts[1,]
counts <- counts[-1,]
counts[,2:ncol(counts)] <- lapply(counts[,2:ncol(counts)], as.numeric)
names(counts)[1] <- "asv"

counts_long_norm <- gather(counts, ID, count, -asv) %>%
  filter(count != 0) %>%
  group_by(ID) %>%
  mutate(reads= sum(count)) %>%
  ungroup() %>%
  mutate(counts_normalized = (count / reads)*100)



seq <- readDNAStringSet("output/asv-seqs-merged.fasta")
pred <- predict(rdp(), seq)
conf <- attr(pred, "confidence")

pred$asv <- row.names(pred)
conf <- as.data.frame(conf)
colnames(conf) <- paste(colnames(conf), "conf", sep = "_")
conf$asv <- row.names(conf)

taxa <- inner_join(pred, conf) %>%
  select(asv, domain, domain_conf, phylum, phylum_conf, class, 
         class_conf, order, order_conf, family, family_conf, genus, genus_conf)


taxa_select <- taxa %>%
  select(asv, phylum, class, order, family)

counts_long_norm_taxa_select <- inner_join(counts_long_norm, taxa_select) %>%
  
  group_by(asv) %>%
  mutate(asv_med = median(counts_normalized)) %>%
  mutate(asv_prev = length(unique(ID))) %>%
  ungroup() %>%
  
  filter((max(asv_prev)/asv_prev) > 0.9) %>%
  filter(asv_med > 1) %>%
  mutate(asv_sub = substr(asv, 1,5)) %>%
  
  unite(asv_class, c("asv_sub", "phylum", "class",  "order", "family"), sep = "_", remove=FALSE)


pdf(file = "plots/top_asv.pdf", width = 7, height = 10)

ggplot(counts_long_norm_taxa_select, aes(x=reorder(asv_class, asv_med), y=counts_normalized)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(shape=16, position=position_jitter(0.2), size = 1, alpha = 0.8) +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust=1))

dev.off()
