library("tidyverse")
library("dplyr")
library("ggpubr")
fill_colors <- c("steelblue","gray83", "skyblue3","lightcyan3", "gray60", "lightsteelblue2", "gray40", "skyblue4", "lightskyblue2", "gray83")


##read in peak height info####
peak_amplitude <- read_delim("/Users/mangofrog7/Documents/Graduate_School/Projects/bioinformatics/K9anti_CnR_5.2022/4.5h_24h_shared_peaks_matrix.tab", delim = '\t', skip = 2) %>% 
  mutate_if(is.factor, ~as.numeric(as.character(.)))

###assign unique peak numbers###
peak_amplitude$ID <- seq.int(nrow(peak_amplitude))

amp_4.5h <- select(peak_amplitude, 1:100, 202)
amp_24h <- select(peak_amplitude, 101:200, 202)

amp_4.5h_scores <- amp_4.5h %>% rowwise(ID) %>% summarise(amp = rowSums(across(where(is.numeric)))) %>% mutate(timepoint = "4.5h")
amp_24h_scores <- amp_24h %>% rowwise(ID) %>% summarise(amp = rowSums(across(where(is.numeric)))) %>% mutate(timepoint = "24h")

total_amp <- bind_rows(amp_4.5h_scores, amp_24h_scores) %>% mutate(timepoint = factor(timepoint, levels = c("4.5h", "24h"))) %>% filter(!is.na(amp))

amp_noOutliers <- total_amp %>% 
  group_by(timepoint) %>% 
  mutate(zPW = scale(amp)) %>% 
  filter(between(zPW,-3,+3)) %>%
  mutate(timepoint = factor(timepoint, levels = c("4.5h", "24h")))

K9_mat_amp_diff <- ggplot(amp_noOutliers, aes(x = timepoint, y = amp)) +
  geom_violin(trim=FALSE, scale="width", aes(fill=timepoint)) + 
  ylim(0,5000) + labs(y = "Average H3K9me3 Signal") + 
  theme_minimal() + scale_fill_manual(values = fill_colors) +
  geom_boxplot(width=0.1, outlier.shape=NA, aes(fill=timepoint)) + 
  theme(legend.position="none") +
  theme(text = element_text(size = 8)) + 
  theme(axis.title.x = element_blank()) 
ggsave("/Users/mangofrog7/Documents/Graduate_School/Papers/K9_Kinetics/FigS10/TEprops_4.5hpf_pie.pdf", plot = K9_mat_amp_diff, width = 2, height = 3, units = "in")

total_amp %>% group_by(timepoint) %>% summarise(mean_amp = mean(amp))
total_amp %>% group_by(timepoint) %>% summarise(med_amp = median(amp))

