library(tidyverse)

df_list <- list()
files <- list.files("/data/Unit_LMM/selberherr-group/strachan/campy/comparative_genomics/annotations/tables/",
 pattern = ".tsv")
i <- 1

for (file in files){
	df_no_ref <- read.delim(paste("/data/Unit_LMM/selberherr-group/strachan/campy/comparative_genomics/annotations/tables/",
		 	file, sep = ""), header = TRUE, sep = '\t')

	df_ref1 <- read.delim(paste("/data/Unit_LMM/selberherr-group/strachan/campy/comparative_genomics/annotations/tables_ref/",
                 	file, sep = ""), header = TRUE, sep = '\t') %>%
		select(locus_tag, gene, product) %>%
		rename(gene2 = gene) %>%
		rename(product2 = product)

	df_ref2 <- read.delim(paste("/data/Unit_LMM/selberherr-group/strachan/campy/comparative_genomics/annotations/tables_ref2/",
                 	file, sep = ""), header = TRUE, sep = '\t') %>%
                select(locus_tag, gene, product) %>%
                rename(gene3 = gene) %>%
                rename(product3 = product)

	
	df_all_refs <- inner_join(df_no_ref, df_ref1) %>%
			inner_join(df_ref2)

	df_all_refs$file <- file

	df_list[[i]] <- df_all_refs
	i <- i + 1
	

}


df_compiled <- bind_rows(df_list) %>%
	select(file, locus_tag, ftype, length_bp, gene, gene2, gene3, product, product2, product3)

write.csv(df_compiled, 'output/compiled_annotations.csv', row.names = FALSE)


