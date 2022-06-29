library(tidyverse)

count_files <- list.files('/data/Unit_LMM/selberherr-group/strachan/campy/mapping/count_tables/', pattern = ".refgenomes.sort.txt")
df_list <- list()
i <- 1

for (file in count_files){
  file_loc <- paste("/data/Unit_LMM/selberherr-group/strachan/campy/mapping/count_tables/", file, sep = "")
  df <- read.delim(file_loc, header=FALSE)
  colnames(df) <- c("ID", "gene", "product", "count")
  df$reads <- gsub(".refgenomes.sort.txt", "", file)
  df_list[[i]] <- df
  i <- i + 1
}

compiled_counts  <- bind_rows(df_list)
write.csv(compiled_counts, 'output/compiled_counts_ref_genomes.csv', row.names = FALSE)



count_files <- list.files('/data/Unit_LMM/selberherr-group/strachan/campy/mapping/count_tables/', pattern = ".clones.sort.txt")
df_list <- list()
i <- 1

for (file in count_files){
  file_loc <- paste("/data/Unit_LMM/selberherr-group/strachan/campy/mapping/count_tables/", file, sep = "")
  df <- read.delim(file_loc, header=FALSE)
  colnames(df) <- c("ID", "gene", "product", "count")
  df$reads <- gsub(".clones.sort.txt", "", file)
  df_list[[i]] <- df
  i <- i + 1
}

compiled_counts  <- bind_rows(df_list)
write.csv(compiled_counts, 'output/compiled_counts_clones.csv', row.names = FALSE)








