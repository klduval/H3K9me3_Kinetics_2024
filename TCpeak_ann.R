library("tidyverse")
library("RColorBrewer")
library("dplyr")
library("tidyr")
library(purrr)
library(ggpubr)

####import annotation data####
col_names <- c("count", "Name")

K9_2hpf <- read_delim("/Users/mangofrog7/Documents/Graduate_School/Projects/bioinformatics/TC_final/peaks/ann2/2hpf_K9_.TEann_35_0.1filt_counts.bed", col_names = col_names) %>% 
 mutate(timepoint = "2 hpf")
K9_2.5hpf <- read_delim("/Users/mangofrog7/Documents/Graduate_School/Projects/bioinformatics/TC_final/peaks/ann2/2.5hpf_K9_.TEann_35_0.1filt_counts.bed", col_names = col_names) %>% 
  mutate(timepoint = "2.5 hpf")
K9_3hpf <- read_delim("/Users/mangofrog7/Documents/Graduate_School/Projects/bioinformatics/TC_final/peaks/ann2/3hpf_K9_.TEann_35_0.1filt_counts.bed", col_names = col_names) %>% 
  mutate(timepoint = "3 hpf")
K9_3.5hpf <- read_delim("/Users/mangofrog7/Documents/Graduate_School/Projects/bioinformatics/TC_final/peaks/ann2/3.5hpf_K9_.TEann_35_0.1filt_counts.bed", col_names = col_names) %>% 
  mutate(timepoint = "3.5 hpf")
K9_4hpf <- read_delim("/Users/mangofrog7/Documents/Graduate_School/Projects/bioinformatics/TC_final/peaks/ann2/4hpf_K9_.TEann_35_0.1filt_counts.bed", col_names = col_names) %>% 
  mutate(timepoint = "4 hpf")
K9_4.5hpf <- read_delim("/Users/mangofrog7/Documents/Graduate_School/Projects/bioinformatics/TC_final/peaks/ann2/4.5hpf_K9_.TEann_35_0.1filt_counts.bed", col_names = col_names) %>% 
  mutate(timepoint = "4.5 hpf")

combTC <- bind_rows(K9_2.5hpf, K9_2hpf, K9_3.5hpf, K9_3hpf, K9_4.5hpf, K9_4hpf)

##importing TE table####
ref_ann <- read_delim("/Users/mangofrog7/Documents/Graduate_School/Projects/bioinformatics/general/TEref_table.txt")  
TE_counts_0.1filt <- read_delim("/Users/mangofrog7/Documents/Graduate_School/Projects/bioinformatics/general/TEann_35_0.1filt_counts.txt", col_names = c("copyNumber_0.1", "Name"))  

ref_ann <- left_join(ref_ann, TE_counts_0.1filt)

##analysing annotation####
big_name_table <- left_join(ref_ann, combTC) %>% mutate(., cn_prop = count / copyNumber_0.1) %>% 
  filter(., Class %in% c("DNA", "LTR", "LINE", "SINE", "RC")) %>% filter(copyNumber_0.1 >= 10) %>% filter(isLTR == "FALSE")

K9_4.5h_only <- big_name_table %>% filter(timepoint == "4.5 hpf")
Class_4.5hpf <- K9_4.5h_only %>% group_by(Class) %>% summarise(sum = sum(count), sum_CN = sum(copyNumber_0.1))
Class_4.5hpf_1 <- Class_4.5hpf %>% mutate(prop = sum / sum_CN, prop_type = "CN")
Class_4.5hpf_2 <- Class_4.5hpf %>% mutate(prop = sum / sum(sum), prop_type = "TE")
Class_4.5hpf <- bind_rows(Class_4.5hpf_2, Class_4.5hpf_1)

####class-level figures####
##lets set a color scheme
fill_colors <- c("steelblue", "lightsteelblue2", "skyblue3", "gray70", "lightcyan3", "gray55", "skyblue4", "lightskyblue2", "gray83")

TE_Props_graph <- ggplot(filter(Class_4.5hpf, prop_type == "TE"), aes(x = Class, y = prop, fill = Class)) + 
  geom_col(position = position_dodge()) + 
  theme_minimal() +
  scale_fill_manual(values = fill_colors) +
  theme(text = element_text(size = 8)) +
  labs(y = "Prop. of H3K9me3-Enriched TEs") +
  theme(legend.position = "none") +
  theme(axis.title.x = element_blank())
ggsave("/Users/mangofrog7/Documents/Graduate_School/Papers/K9_Kinetics/Fig2/TEprops_4.5hpf.pdf", plot = TE_Props_graph, width = 2, height = 2, units = "in")

CN_Props_graph <- ggplot(filter(Class_4.5hpf, prop_type == "CN"), aes(x = Class, y = prop, fill = Class)) + 
  geom_col(position = position_dodge()) + 
  theme_minimal() +
  scale_fill_manual(values = fill_colors) +
  theme(text = element_text(size = 8)) +
  labs(y = "Prop. of Total TE Copies") +
  theme(legend.position = "none") +
  theme(axis.title.x = element_blank())
ggsave("/Users/mangofrog7/Documents/Graduate_School/Papers/K9_Kinetics/Fig2/CNprops_4.5hpf.pdf", plot = CN_Props_graph, width = 2, height = 2, units = "in")


######family level####
Family_4.5hpf <- K9_4.5h_only %>% group_by(Class, Family, Name) %>% summarise(sum = sum(count), sum_CN = sum(copyNumber_0.1))
Family_4.5hpf <- Family_4.5hpf %>% mutate(prop = sum / sum_CN)
Family_20p_cutoff <- Family_4.5hpf %>% group_by(Class) %>% summarise(yes = sum(prop >= 0.2), no = sum(prop < 0.2))

Family_20p_cutoff_grpah <- ggplot(Family_20p_cutoff, aes(x = Class, y = yes, fill = Class)) + 
  geom_col(position = position_dodge()) + 
  theme_minimal() +
  scale_fill_manual(values = fill_colors) +
  theme(text = element_text(size = 8)) +
  labs(y = "Prop. Subfamilies w/ 20% Copies Enriched for H3K9me3") +
  theme(legend.position = "none") +
  theme(axis.title.x = element_blank())
ggsave("/Users/mangofrog7/Documents/Graduate_School/Papers/K9_Kinetics/Fig2/Family_20p_cutoff.pdf", plot = Family_20p_cutoff_grpah, width = 2, height = 2, units = "in")
ggsave("/Users/mangofrog7/Documents/Graduate_School/Papers/K9_Kinetics/Fig2/Family_20p_cutoff_bigger.pdf", plot = Family_20p_cutoff_grpah, width = 3, height = 3, units = "in")

####LTR families####
LTRs_20p_Fams <- Family_4.5hpf %>% filter(Class == "LTR", Family != "Gypsy?") %>% group_by(Family) %>% summarise(yes = sum(prop >= 0.2), no = sum(prop < 0.2))
LTRs_20props <- Family_4.5hpf %>% filter(Class == "LTR", Family != "Gypsy?") %>% group_by(Family) %>% 
  summarise(count = sum(sum), CN = sum(sum_CN)) %>% mutate(prop = count / CN)

LTRs_20props_pics <- ggplot(LTRs_20props, aes(x = Family, y = prop, fill = Family)) + 
  geom_col(position = position_dodge()) + 
  theme_minimal() +
  scale_fill_manual(values = fill_colors) +
  theme(text = element_text(size = 8)) +
  labs(y = "Prop. Copies Enriched for H3K9me3") +
  theme(legend.position = "none") +
  theme(axis.title.x = element_blank())
ggsave("/Users/mangofrog7/Documents/Graduate_School/Papers/K9_Kinetics/Fig2/LTRs_20props_pics.pdf", plot = LTRs_20props_pics, width = 3.5, height = 2, units = "in")

####is proportion of copies enriched related to TE age####
cnprop_medlen_scatter <- ggplot(K9_4.5h_only, aes(x = cn_prop, y = MedLen)) + 
  geom_point(size = 0.5) + geom_smooth(method = "lm", se=F, linewidth = 0.5) +
  facet_wrap(vars(Class), nrow = 1) +  
  theme_minimal() +  theme(text = element_text(size = 8)) + stat_cor(size = 2) +
  labs(y = "Median Branch Length", x = "Proportion of Copies Enriched for H3K9me3") +
  theme(legend.position = "none") + ylim(0,0.5)
ggsave("/Users/mangofrog7/Documents/Graduate_School/Papers/K9_Kinetics/Fig2/cnprop_medlen_scatter.pdf", plot = cnprop_medlen_scatter, width = 5.5, height = 2, units = "in")
