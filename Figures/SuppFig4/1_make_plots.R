library(tidyverse)

df_all_v_all_long <- read.delim("~/RumenCampylobacter2022/Processing/genomes/output/all_v_all_ani.txt", header=FALSE)
df_all_v_all_long$V4 <- NULL
df_all_v_all_long$V5 <- NULL
colnames(df_all_v_all_long) <- c("genome1", "genome2", "ani")

## clean up file names
df_all_v_all_long$genome1 <- gsub("/data/Unit_LMM/selberherr-group/strachan/campy/comparative_genomics/assemblies/", "", df_all_v_all_long$genome1)

df_all_v_all_long$genome2 <- gsub("/data/Unit_LMM/selberherr-group/strachan/campy/comparative_genomics/assemblies/", "", df_all_v_all_long$genome2)

df_all_v_all_long$genome1 <- gsub("/data/Unit_LMM/selberherr-group/strachan/campy/metagenomes/assemblies/selected_bins/", "", df_all_v_all_long$genome1)

df_all_v_all_long$genome2 <- gsub("/data/Unit_LMM/selberherr-group/strachan/campy/metagenomes/assemblies/selected_bins/", "", df_all_v_all_long$genome2)

df_all_v_all_long$genome1 <- gsub("\\.fasta", "", df_all_v_all_long$genome1)

df_all_v_all_long$genome2 <- gsub("\\.fasta", "", df_all_v_all_long$genome2)

df_all_v_all_wide <- df_all_v_all_long %>%
  
  spread(genome2, ani)

dist <- dist(t(as.matrix(df_all_v_all_wide[,2:ncol(df_all_v_all_wide)])))
order <- hclust(dist)$order
order_names <- names(df_all_v_all_wide[,2:ncol(df_all_v_all_wide)])[order]

df_all_v_all_long$genome1 <- factor(df_all_v_all_long$genome1, levels = order_names)
df_all_v_all_long$genome2 <- factor(df_all_v_all_long$genome2, levels = order_names)

pdf("~/RumenCampylobacter2022/Figures/SuppFig4/output/supp_fig4.pdf", width=24, height=24)

ggplot(df_all_v_all_long, aes(genome1, genome2) ) +
  geom_tile(aes(fill = ani)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 18), 
        axis.text.y = element_text(size = 20), 
        axis.title.x = element_blank(), 
        axis.title.y = element_blank()) + 
  theme(legend.key.size = unit(2, "cm"), 
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 30)) + 
  coord_fixed() +
  scale_fill_gradient2(low = "white", high = "deepskyblue2", midpoint = 90, limits = c(90,100)) + 
  geom_text(aes(label = round(ani, 1))) 

dev.off()

pdf("~/RumenCampylobacter2022/Figures/SuppFig4/output/supp_fig4_dendogram.pdf", width=6, height=6)

plot(as.dendrogram(hclust(dist)))

dev.off()
