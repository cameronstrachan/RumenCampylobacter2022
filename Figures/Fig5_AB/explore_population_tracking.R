library(tidyverse)
library(readxl)

qPCR_16s <- read_excel("~/Desktop/DATA_TO_ANALYZE/cleaned_population_tracking.xlsx", sheet = "16S_qPCR")
dPCR_pops <- read_excel("~/Desktop/DATA_TO_ANALYZE/cleaned_population_tracking.xlsx", sheet = "Pop_dPCR")

meta <- read.csv("~/Desktop/DATA_TO_ANALYZE/adda_meta_data_mod.csv")

meta_clean <- meta %>%
  
  gather(SCFA, Concentration, -ID, -Timepoint, -Cow, -Additive) %>%
  
  group_by(ID, SCFA) %>%
  mutate(Min = min(Concentration)) %>%
  ungroup() %>%
  
  select(-Concentration) %>%
  
  gather(Metric, Concentration, -ID, -Timepoint, -Cow, -Additive, -SCFA) %>%
  
  unite(SCFA, c("SCFA", "Metric"), sep = "_") %>%
  
  select(-Timepoint, -Cow) %>%
  
  distinct() %>%
  
  spread(SCFA, Concentration)

meta_clean$ID <- as.character(meta_clean$ID)  

qPCR_16s_clean <- qPCR_16s %>%
  distinct()

compiled <- full_join(qPCR_16s, dPCR_pops) %>%
  
  filter(!(is.na(SampleName))) %>%
  distinct() %>%
  
  mutate(`16S_total_copies` = as.numeric(gsub("No Ct", 0, `16S_total_copies`))) %>%
  separate(ID, into=c("Run", "Day", "Cow"), remove = FALSE) %>%
  select(-SampleName) 

compiled_long <- compiled %>%
  
  gather(Measurement, Copies, -ID, -Run, -Day, -Cow) %>%
  
  distinct() %>%
  
  group_by(Run, Day, Measurement) %>%
  mutate(n = n()) %>%
  mutate(Median = median(Copies)) %>%
  mutate(SD = sd(Copies)) %>%
  ungroup()

ggplot(compiled_long, aes(x=Day, y=Median)) +
  geom_bar(stat="identity", color="black", position=position_dodge()) +
  facet_grid(Measurement ~ Run, scales = "free") +
  geom_errorbar(aes(ymin=Median, ymax=Median+SD), width=.2,position=position_dodge(.9))

population_raio_variance <- compiled_long %>%
  filter(Measurement != "16S_total_copies")

ggplot(population_raio_variance, aes(x=Day, y=Copies)) +
  geom_point(aes(colour=Measurement)) +
  scale_y_continuous(trans='log2') +
  facet_grid(Cow ~ Run) + 
  theme_bw()


library(corrplot)


  
  

correlation <- inner_join(compiled, meta_clean) %>%
  select(-`16S_total_copies`, -Butyrate_uM_Min, -Caproate_uM_Min, -Valerate_uM_Min) %>%
  select(-ID, -Run, -Day, -Cow, -Additive)

colnames(correlation) <- gsub("_uM_Min", "", names(correlation))

corr_mat <- cor(correlation)
corr_mat_testRes = cor.mtest(correlation, conf.level = 0.95)

corrplot(corr_mat, p.mat=corr_mat_testRes$p, method = 'color', diag = FALSE, type = 'upper',
         sig.level = c(0.001, 0.01, 0.05), pch.cex = 0.9,
         insig = 'label_sig')


corrplot(corr_mat, type="upper", order="original", p.mat = corr_mat_testRes$p, sig.level = 0.01, insig = "blank", diag = FALSE)

corrplot(corr_mat, method="circle", col = colorRampPalette(c("red", "black"))(50),
         type="lower", order="original", 
         addCoef.col = "white", # Add coefficient of correlation
         tl.col="black", tl.srt=45, #Text label color and rotation
         # Combine with significance
         p.mat = corr_mat_testRes$p, sig.level = 0.01, insig = "pch", 
         # hide correlation coefficient on the principal diagonal
         diag=FALSE)
)

