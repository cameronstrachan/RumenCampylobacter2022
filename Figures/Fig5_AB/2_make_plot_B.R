library(tidyverse)

qPCR_16s <- read_csv("~/RumenCampylobacter2022/Raw/dPCR_feeding_trial.csv")
dPCR_pops <- read_csv("~/RumenCampylobacter2022/Raw/qPCR_feeding_trial.csv")

meta <- read.csv("~/RumenCampylobacter2022/MetaData/metadata_feeding_trial.csv")

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

library(corrplot)

correlation <- inner_join(compiled, meta_clean) %>%
  select(-`16S_total_copies`, -Butyrate_uM_Min, -Caproate_uM_Min, -Valerate_uM_Min) %>%
  select(-ID, -Run, -Day, -Cow, -Additive)

colnames(correlation) <- gsub("_uM_Min", "", names(correlation))

corr_mat <- cor(correlation)
corr_mat_testRes = cor.mtest(correlation, conf.level = 0.95)

pdf("~/RumenCampylobacter2022/Figures/Fig5_AB/output/fig5_B.pdf", width=6, height=6)

corrplot(corr_mat, method="circle", col = colorRampPalette(c("red", "black"))(50),
         type="lower", order="original", 
         addCoef.col = "white", # Add coefficient of correlation
         tl.col="black", tl.srt=45, #Text label color and rotation
         # Combine with significance
         p.mat = corr_mat_testRes$p, sig.level = 0.01, insig = "pch", 
         # hide correlation coefficient on the principal diagonal
         diag=FALSE)

dev.off()
