library("tidyverse")
library("RColorBrewer")
library("dplyr")
library("ggpubr")

# Set your working directory to the "sats" directory
setwd("/Users/mangofrog7/Documents/Graduate_School/Projects/bioinformatics/TC_final/sats")

###import files 
SAT1 <- read_delim("4.5hpf_K9_AVG.BRSAT1.tab", delim = '\t', skip = 2)
control <- read_delim("4.5hpf.SAT1_control.tab", delim = '\t', skip = 2)

###sum across rows
BRSAT1 <- SAT1 %>% rowwise() %>% mutate(total = rowSums(across(where(is.numeric)))) %>%
  mutate(region = "BRSAT1")
Control_df <- control %>% rowwise() %>% mutate(total = rowSums(across(where(is.numeric)))) %>%
  mutate(region = "Control")

combined_df <- bind_rows(BRSAT1, Control_df) %>% select(total, region)

###figure 
fill_colors <- c("steelblue","gray83", "skyblue3","lightcyan3", "gray60", "lightsteelblue2", "gray40", "skyblue4", "lightskyblue2")

enrich_trimmed <- combined_df %>% 
  group_by(region) %>% 
  mutate(zPW = scale(total)) %>% 
  filter(between(zPW,-3,+3))

SAT1_enrichment <- ggplot(enrich_trimmed, aes(x = region, y = total, fill=region)) + 
  geom_violin(trim=FALSE, scale="width") + 
  labs(y = "Average H3K9me3 Signal") + 
  theme_minimal() + scale_fill_manual(values = fill_colors) +
  geom_boxplot(width=0.1, outlier.shape=NA) + theme(legend.position="none") +
  theme(axis.title.x = element_blank()) + ylim(-10000, 165000) +
  theme(text = element_text(size = 8), axis.text.y = element_text(size=5)) + scale_y_continuous(labels = function(x) format(x, scientific = TRUE)) +
  stat_compare_means(method = "t.test", label = "p.format", label.x = 1.5, size = 2)
ggsave("/Users/mangofrog7/Documents/Graduate_School/Papers/K9_Kinetics/Fig1/SAT1_enrichment.pdf", plot = SAT1_enrichment, width = 1.75, height = 2.25, units = "in")
