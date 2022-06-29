library(tidyverse)

class <- read.delim("~/RumenCampylobacter2022/Processing/metagenomes/output/gtdbtk.bac120.summary.tsv")[,1:2]
colnames(class)[1] <- "bin_id"

checkm <- read_table("~/RumenCampylobacter2022/Processing/metagenomes/output/checkM.txt", col_names = FALSE, skip = 3, comment = "---")[,c(1,13:15)]
colnames(checkm) <- c("bin_id", "completeness", "contamination", "strain_heterogeneity")

map <- read.csv("~/RumenCampylobacter2022/Processing/metagenomes/output/bin_contig_map.csv", header=FALSE)
colnames(map) <- c("bin_id", "contig_id")

read_count_files <- list.files('~/RumenCampylobacter2022/Processing/metagenomes/output/mapping/', pattern = "bins_readcounts.txt")
df_list <- list()
i <- 1

for (file in read_count_files){
  file_loc <- paste("~/RumenCampylobacter2022/Processing/metagenomes/output/mapping/", file, sep = "")
  df <- read.delim(file_loc, header=FALSE)[,1:3]
  colnames(df) <- c("contig_id", "contig_len", "count")
  df$reads <- gsub("_bins_readcounts.txt", "", file)
  df_list[[i]] <- df
  i <- i + 1
}

read_count <- bind_rows(df_list)

compiled <- inner_join(class, checkm) %>%
  inner_join(map) %>%
  inner_join(read_count) 

summary <- compiled %>%

  filter(count > 3) %>%
	
  separate(reads, into=c("sample", "lane"), sep = "_") %>%
  
  group_by(bin_id, sample) %>%
  mutate(bin_len = sum(contig_len)) %>%
  mutate(sample_count = sum(count)) %>%
  ungroup() %>%
  
  group_by(bin_id) %>%
  mutate(total_count = sum(count)) %>%
  ungroup() %>%

  mutate(sample_med_count_per_kb = (sample_count/bin_len)*1000) %>%
  mutate(total_med_count_per_kb = (total_count/bin_len)*1000) %>%
  
  select(-count, -contig_id, -contig_len, -lane) %>%
  distinct() %>%
  separate(classification, into = c("domain", "phylum", "class", "order", "family", "genus"), sep = ";") %>%
  unite(id, c("bin_id", "phylum", "class", "order", "family"), sep = "_", remove = FALSE) %>%
  filter(contamination < 10) %>%
  filter(completeness > 50) 

pdf(file = "~/RumenCampylobacter2022/Figures/Fig1_AB/output/fig1_B.pdf",   # The directory you want to save the file in
    width = 7, # The width of the plot in inches
    height = 9)

ggplot(summary, aes(x=reorder(id, total_med_count_per_kb), y=sample_med_count_per_kb)) + 
  geom_boxplot(outlier.shape = NA) +  
  geom_point(shape=16, size = 1) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

dev.off()
