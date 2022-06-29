library(tidyverse)

###

#df_list <- list()
#files <- list.files("/data/Unit_LMM/selberherr-group/strachan/campy/comparative_genomics/blast_intermediate_files/genes_table_select/", pattern = ".txt")
#i <- 1
#
#for (file in files){
#  df <- read_delim(paste("/data/Unit_LMM/selberherr-group/strachan/campy/comparative_genomics/blast_intermediate_files/genes_table_select/", file, sep = ""), "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
#  file <- gsub(".txt", "", file)
#  file <- gsub(".ffn", "", file)
#  
#  df$file <- file
#  
#  df_list[[i]] <- df
#  i <- i + 1
#}
#
#compiled_blast_hits_nucl <- bind_rows(df_list)
#colnames(compiled_blast_hits_nucl)[1:14] <- c("gene1", "gene2", "pident", "sstart", "send", "qstart", "qend", "evalue", "bitscore", "score", "qlen", "length", "sseq", "file")
#
#compiled_blast_hits_nucl_trim  <- compiled_blast_hits_nucl %>%
#  
#  select(-sseq) %>%
#  separate(file, into = c("genome1", "genome2"), sep = "_V_") %>%
#  mutate(palign = (length / qlen)*100)
#
#compiled_blast_hits_nucl_trim1 <- select(compiled_blast_hits_nucl_trim, genome1, gene1, genome2, gene2, pident, palign)
#colnames(compiled_blast_hits_nucl_trim1)[5] <- "pident1"
#colnames(compiled_blast_hits_nucl_trim1)[6] <- "palign1"
#
#compiled_blast_hits_nucl_trim2 <- select(compiled_blast_hits_nucl_trim, genome1, gene1, genome2, gene2, pident, palign) 
#colnames(compiled_blast_hits_nucl_trim2) <- c("genome2", "gene2", "genome1", "gene1", "pident2", "palign2")
#
#reciprocal_best_hits_nucl <- inner_join(compiled_blast_hits_nucl_trim1, compiled_blast_hits_nucl_trim2) %>%
#  mutate(pident_nucl = (pident1 + pident2)/2) %>%
#  mutate(palign_nucl = (palign1 + palign2)/2) %>%
#  select(-pident1, -pident2, -palign1, -palign2)
#
#write.csv(reciprocal_best_hits_nucl, 'output/reciprocal_best_hits_nucl_select.csv', row.names = FALSE)

###

df_list <- list()
files <- list.files("/data/Unit_LMM/selberherr-group/strachan/campy/comparative_genomics/blast_intermediate_files/prots_table_select/", pattern = ".txt")
i <- 1

for (file in files){
  df <- read_delim(paste("/data/Unit_LMM/selberherr-group/strachan/campy/comparative_genomics/blast_intermediate_files/prots_table_select/", file, sep = ""), "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
  file <- gsub(".txt", "", file)
  file <- gsub(".faa", "", file)

  df$file <- file
  
  df_list[[i]] <- df
  i <- i + 1
}

compiled_blast_hits_prot <- bind_rows(df_list)
colnames(compiled_blast_hits_prot)[1:14] <- c("gene1", "gene2", "pident", "sstart", "send", "qstart", "qend", "evalue", "bitscore", "score", "qlen", "length", "sseq", "file")

compiled_blast_hits_prot_trim <- compiled_blast_hits_prot %>%
  
  select(-sseq) %>%
  separate(file, into = c("genome1", "genome2"), sep = "_V_") %>%
  mutate(palign = (length / qlen)*100) 

compiled_blast_hits_prot_trim1 <- select(compiled_blast_hits_prot_trim, genome1, gene1, genome2, gene2, pident, palign)
colnames(compiled_blast_hits_prot_trim1)[5] <- "pident1"
colnames(compiled_blast_hits_prot_trim1)[6] <- "palign1"

compiled_blast_hits_prot_trim2 <- select(compiled_blast_hits_prot_trim, genome1, gene1, genome2, gene2, pident, palign) 
colnames(compiled_blast_hits_prot_trim2) <- c("genome2", "gene2", "genome1", "gene1", "pident2", "palign2")

reciprocal_best_hits_prot <- inner_join(compiled_blast_hits_prot_trim1, compiled_blast_hits_prot_trim2) %>%
  mutate(pident_prot = (pident1 + pident2)/2) %>%
  mutate(palign_prot = (palign1 + palign2)/2) %>%
  select(-pident1, -pident2, -palign1, -palign2)

write.csv(reciprocal_best_hits_prot, 'output/reciprocal_best_hits_prot_select.csv', row.names = FALSE)


