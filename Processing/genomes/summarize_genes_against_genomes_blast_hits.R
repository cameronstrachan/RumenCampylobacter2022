library(tidyverse)

###

df_list <- list()
files <- list.files("/data/Unit_LMM/selberherr-group/strachan/campy/comparative_genomics/blast_intermediate_files/genes_against_genomes/", pattern = ".txt")
i <- 1

for (file in files){
  df <- read_delim(paste("/data/Unit_LMM/selberherr-group/strachan/campy/comparative_genomics/blast_intermediate_files/genes_against_genomes/", file, sep = ""), "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
  file <- gsub(".txt", "", file)
  file <- gsub(".ffn", "", file)
  
  df$file <- file
  
  df_list[[i]] <- df
  i <- i + 1
}

compiled_blast_hits_nucl <- bind_rows(df_list)
colnames(compiled_blast_hits_nucl)[1:14] <- c("gene1", "gene2", "pident", "sstart", "send", "qstart", "qend", "evalue", "bitscore", "score", "qlen", "length", "sseq", "file")

compiled_blast_hits_nucl_trim  <- compiled_blast_hits_nucl %>%
  
  select(-sseq) %>%
  separate(file, into = c("genome1", "genome2"), sep = "_V_") %>%
  mutate(palign = (length / qlen)*100)

write.csv(compiled_blast_hits_nucl_trim, 'output/genes_against_genomes.csv', row.names = FALSE)

###

