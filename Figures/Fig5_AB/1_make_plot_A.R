library(tidyverse)

compiled <- read.csv("~/RumenCampylobacter2022/PCRdata/dPCR_in_vitro.csv") %>%
  group_by(Treatment, Measure, Inoculation) %>%
  mutate(Avg=mean(Copies)) %>%
  mutate(Sd =sd(Copies)) %>%
  ungroup() %>%
  
  select(Treatment, Measure, Inoculation, Avg, Sd) %>%
  
  distinct()

compiled_treatment <- compiled %>%
  filter(Treatment != "Control") 

compiled_control <- compiled %>%
  filter(Treatment == "Control") %>%
  rename(Avg_control = Avg) %>%
  rename(Sd_control = Sd) %>%
  select(-Treatment)

fold_change <- inner_join(compiled_treatment, compiled_control) %>%
  mutate(fold = Avg / Avg_control) %>%
  mutate(fold_sd = Sd / Avg_control)

pdf("~/RumenCampylobacter2022/Figures/Fig5_AB/output/fig5_A.pdf", width=9, height=4.5)

ggplot(fold_change , aes(x=Treatment, y=fold, fill=Measure)) + 
  
  geom_bar(stat="identity", color="black", 
           position=position_dodge(), alpha = 0.5) +
  
  geom_errorbar(aes(ymin=fold-fold_sd, ymax=fold+fold_sd), width=.2,
                position=position_dodge(.9)) + 
  theme_bw() + 
  facet_grid(. ~ Inoculation, scales = "free_x", space="free") 

dev.off()

# changes to correspond to nature policy

compiled2 <- read.csv("~/RumenCampylobacter2022/PCRdata/dPCR_in_vitro.csv") %>%
  group_by(Treatment, Measure, Inoculation) %>%
  mutate(Avg=mean(Copies)) %>%
  mutate(Sd =sd(Copies)) %>%
  ungroup() %>%
  
  select(Treatment, Measure, Inoculation, Avg, Sd, Copies) %>%
  
  distinct()

compiled_treatment2 <- compiled2 %>%
  filter(Treatment != "Control")

compiled_control2 <- compiled2 %>%
  filter(Treatment == "Control") %>%
  rename(Avg_control = Avg) %>%
  rename(Sd_control = Sd) %>%
  select(-Treatment, -Copies) %>%
  distinct()

fold_change2 <- inner_join(compiled_treatment2, compiled_control2) %>%
  mutate(fold = Avg / Avg_control) %>%
  mutate(fold_sd = Sd / Avg_control) %>%
  mutate(fold_ind = Copies / Avg_control)

pdf("~/RumenCampylobacter2022/Figures/Fig5_AB/output/fig5_A_v2.pdf", width=9, height=4.5)

ggplot(fold_change2, aes(x=Treatment, y=fold, fill=Measure)) + 
  
  geom_bar(stat="identity", color="black", 
           position=position_dodge(), alpha = 0.5) +
  
  geom_errorbar(aes(ymin=fold-fold_sd, ymax=fold+fold_sd), width=.2,
                position=position_dodge(.9)) + 
  theme_bw() + 
  facet_grid(. ~ Inoculation, scales = "free_x", space="free") + 
  geom_jitter(aes(y=fold_ind), position = position_jitterdodge(0.5), alpha=0.4)

dev.off()
