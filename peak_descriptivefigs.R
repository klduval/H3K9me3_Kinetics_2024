library("tidyverse")
library("RColorBrewer")
library("dplyr")
library("ggpubr")

####import data####
bed_col_names <- c("Chr", "Start", "End", "name", "none", "none_1", "none_2", "none_3", "none_4")

p2hpf <- read_delim("/Users/mangofrog7/Documents/Graduate_School/Projects/bioinformatics/TC_final/peaks/2hpf_K9_final.bed", col_names = bed_col_names) %>% 
  mutate(time_point = "2hpf")
p2.5hpf <- read_delim("/Users/mangofrog7/Documents/Graduate_School/Projects/bioinformatics/TC_final/peaks/2.5hpf_K9_final.bed", col_names = bed_col_names) %>% 
  mutate(time_point = "2.5hpf")
p3hpf <- read_delim("/Users/mangofrog7/Documents/Graduate_School/Projects/bioinformatics/TC_final/peaks/3hpf_K9_final.bed", col_names = bed_col_names) %>% 
  mutate(time_point = "3hpf")
p3.5hpf <- read_delim("/Users/mangofrog7/Documents/Graduate_School/Projects/bioinformatics/TC_final/peaks/3.5hpf_K9_final.bed", col_names = bed_col_names) %>% 
  mutate(time_point = "3.5hpf")
p4hpf <- read_delim("/Users/mangofrog7/Documents/Graduate_School/Projects/bioinformatics/TC_final/peaks/4hpf_K9_final.bed", col_names = bed_col_names) %>% 
  mutate(time_point = "4hpf")
p4.5hpf <- read_delim("/Users/mangofrog7/Documents/Graduate_School/Projects/bioinformatics/TC_final/peaks/4.5hpf_K9_final.bed", col_names = bed_col_names) %>% 
  mutate(time_point = "4.5hpf")

###combining timepoints####
time_points_combined_rows <- bind_rows(p2hpf, p2.5hpf, p3hpf, p3.5hpf, p4hpf, p4.5hpf)

time_points_combined_rows <- time_points_combined_rows %>% 
  mutate(time_point = factor(time_point, levels = c("2hpf", "2.5hpf", "3hpf", "3.5hpf", "4hpf", "4.5hpf")))

time_points_combined_rows.combined <- time_points_combined_rows %>%
  mutate(peak_width = End - Start) %>% 
  group_by(time_point) %>% 
  summarise(mean_peak_width = mean(peak_width), 
            n_peaks = n())

####making descriptive figs####
npeaks_bar <- ggplot(time_points_combined_rows.combined, aes(x = time_point, y = n_peaks, fill = time_point)) + 
  geom_bar(stat = "identity", color = "grey50") + scale_fill_brewer(palette = "Blues") + theme_minimal() + 
  ylab("H3K9me3 Peaks") + theme(text = element_text(size = 8)) + theme(legend.position = "none")
ggsave("/Users/mangofrog7/Documents/Graduate_School/Papers/K9_Kinetics/Fig3/npeaks_bar.pdf", plot = npeaks_bar, width = 2, height = 2.5, units = "in")


###peak_width
pw_noOutliers <- time_points_combined_rows.combined %>% 
  group_by(time_point) %>% 
  mutate(zPW = scale(peak_width)) %>% 
  filter(between(zPW,-5,+5))

peakwidth_viol <- ggplot(pw_noOutliers, aes(x = time_point, y = peak_width, fill=time_point)) + 
  geom_violin(trim=FALSE, scale="width") + 
  labs(y = "Peak Width (bps)")+ 
  theme_minimal() + scale_fill_brewer(palette="Blues") +
  geom_boxplot(width=0.1, outlier.shape=NA) + theme(legend.position="none") +
  theme(text = element_text(size = 8)) 
ggsave("/Users/mangofrog7/Documents/Graduate_School/Papers/K9_Kinetics/Fig3/peakwidth_viol.pdf", plot = peakwidth_viol, width = 2, height = 2, units = "in")


###now going to read in the peak enrichment matrices so I can calculate peak heights at each time####
h2hpf <- read_delim("/Users/mangofrog7/Documents/Graduate_School/Projects/bioinformatics/TC_final/peaks/2hpf_peakmatrix.tab", delim = '\t', skip = 2) %>% 
    mutate_if(is.factor, ~as.numeric(as.character(.)))
h2.5hpf <- read_delim("/Users/mangofrog7/Documents/Graduate_School/Projects/bioinformatics/TC_final/peaks/2.5hpf_peakmatrix.tab", delim = '\t', skip = 2) %>%  
  mutate_if(is.factor, ~as.numeric(as.character(.)))
h3hpf <- read_delim("/Users/mangofrog7/Documents/Graduate_School/Projects/bioinformatics/TC_final/peaks/3hpf_peakmatrix.tab", delim = '\t', skip = 2) %>%  
  mutate_if(is.factor, ~as.numeric(as.character(.)))
h3.5hpf <- read_delim("/Users/mangofrog7/Documents/Graduate_School/Projects/bioinformatics/TC_final/peaks/3.5hpf_peakmatrix.tab", delim = '\t', skip = 2) %>% 
  mutate_if(is.factor, ~as.numeric(as.character(.)))
h4hpf <- read_delim("/Users/mangofrog7/Documents/Graduate_School/Projects/bioinformatics/TC_final/peaks/4hpf_peakmatrix.tab", delim = '\t', skip = 2) %>%  
  mutate_if(is.factor, ~as.numeric(as.character(.)))
h4.5hpf <- read_delim("/Users/mangofrog7/Documents/Graduate_School/Projects/bioinformatics/TC_final/peaks/4.5hpf_peakmatrix.tab", delim = '\t', skip = 2) %>% 
  mutate_if(is.factor, ~as.numeric(as.character(.)))

h2hpf_scores <- h2hpf %>% rowwise() %>% mutate(total = rowSums(across(where(is.numeric)))) %>%
  mutate(time_point = "2hpf") %>% mutate(EGA = "pre")
as_tibble(rowSums(h2hpf)) %>% mutate(time_point = "2hpf")
h2.5hpf_scores <- h2.5hpf %>% rowwise() %>% mutate(total = rowSums(across(where(is.numeric)))) %>% 
  mutate(time_point = "2.5hpf") %>% mutate(EGA = "pre")
h3hpf_scores <- h3hpf %>% rowwise() %>% mutate(total = rowSums(across(where(is.numeric)))) %>% 
  mutate(time_point = "3hpf") %>% mutate(EGA = "pre")
h3.5hpf_scores <- h3.5hpf %>% rowwise() %>% mutate(total = rowSums(across(where(is.numeric)))) %>% 
  mutate(time_point = "3.5hpf") %>% mutate(EGA = "post")
h4hpf_scores <- h4hpf %>% rowwise() %>% mutate(total = rowSums(across(where(is.numeric)))) %>% 
  mutate(time_point = "4hpf") %>% mutate(EGA = "post")
h4.5hpf_scores <- h4.5hpf %>% rowwise() %>% mutate(total = rowSums(across(where(is.numeric)))) %>% 
  mutate(time_point = "4.5hpf") %>% mutate(EGA = "post")

height_scores_combined <- bind_rows(h2hpf_scores, h2.5hpf_scores, h3hpf_scores, h3.5hpf_scores, h4hpf_scores, h4.5hpf_scores)

height_noOutliers <- height_scores_combined %>% 
  group_by(time_point) %>% 
  mutate(zPW = scale(total)) %>% 
  filter(between(zPW,-5,+5)) %>%
  mutate(time_point = factor(time_point, levels = c("2hpf", "2.5hpf", "3hpf", "3.5hpf", "4hpf", "4.5hpf")))

t_test_height <- compare_means(total ~ EGA, p.adjust.method = "bonferroni", method='t.test', data = height_noOutliers)  
t_test_height <- t_test_height %>% mutate(y.position = c(40000))
##p value is 0 but will add manually in illustrator

peak_height <- ggplot(height_noOutliers, aes(x = time_point, y = total)) + 
  geom_violin(trim=FALSE, scale="width", aes(fill=time_point)) + 
  labs(y = "Average H3K9me3 Signal") + 
  theme_minimal() + scale_fill_brewer(palette="Blues") +
  geom_boxplot(width=0.1, outlier.shape=NA, aes(fill=time_point)) + theme(legend.position="none") +
  theme(text = element_text(size = 8)) + ylim(0,45000) +
  theme(axis.title.x = element_blank()) 
  stat_pvalue_manual(t_test_height, label = "p.adj", hide.ns = TRUE, tip.length = 0.0, size = 2.5)

ggsave("/Users/mangofrog7/Documents/Graduate_School/Papers/K9_Kinetics/Fig3/peak_height.pdf", plot = peak_height, width = 2, height = 2.5, units = "in")
