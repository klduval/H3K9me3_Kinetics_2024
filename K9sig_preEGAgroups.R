library("tidyverse")
library("RColorBrewer")
library("dplyr")
library("ggpubr")
fill_colors <- c("steelblue","gray83", "skyblue3","lightcyan3", "gray60", "lightsteelblue2", "gray40", "skyblue4", "lightskyblue2", "gray83")

####annotations of these groups####
col_names <- c("count", "Name")

preEGA_neither <- read_delim("~/peaks/preEGA_K9peaks_noperic_nopiRNA_TEcounts.bed", col_names = col_names) %>% 
  mutate(cat = "neither")
preEGA_peric_only <- read_delim("~/peaks/preEGA_K9peaks_peric_only_TEcounts.bed", col_names = col_names) %>% 
  mutate(cat = "Peric")
preEGA_piRNA_only <- read_delim("~/peaks/preEGA_K9peaks_piRNA_only_TEcounts.bed", col_names = col_names) %>% 
  mutate(cat = "piRNA")
preEGA_both <- read_delim("~/peaks/preEGA_K9peaks_piRNA_N_peric_TEcounts.bed", col_names = col_names) %>% 
  mutate(cat = "Peric&piRNA")

total_ann <- bind_rows(preEGA_neither, preEGA_both, preEGA_peric_only, preEGA_piRNA_only)

###reftable
ref_ann <- read_delim("~/general/TEref_table.txt")  
TE_counts_0.1filt <- read_delim("~/general/TEann_35_0.1filt_counts.txt", col_names = c("copyNumber_0.1", "Name"))  

ref_ann <- left_join(ref_ann, TE_counts_0.1filt)

ref_w_preEGAcats <- left_join(ref_ann, total_ann) %>% mutate(., cn_prop = count / copyNumber_0.1) %>% 
  filter(., Class %in% c("DNA", "LTR", "LINE", "SINE")) %>% filter(copyNumber_0.1 >= 10) %>% filter(isLTR == "FALSE")

Class_4.5hpf <- ref_w_preEGAcats %>% group_by(Class) %>% 
  summarise(sum_CN = sum(copyNumber_0.1)) 

Class_4.5hpf_2 <- ref_w_preEGAcats %>% group_by(Class, cat) %>% summarise(sum = sum(count))

Class_4.5hpf <- left_join(Class_4.5hpf_2, Class_4.5hpf) %>% mutate(prop = sum / sum_CN) %>% filter(!is.na(prop))

ggplot(Class_4.5hpf, aes(x = cat, y = sum, fill = cat)) + 
  geom_col(position = position_dodge()) + facet_wrap(~Class, scales = "free", nrow = 1) +
  theme_minimal() + 
  scale_fill_manual(values = fill_colors) +
  theme(text = element_text(size = 8)) +
  labs(y = "Elements with preEGA H3K9me3") +
  theme(legend.position = "none") +
  theme(axis.title.x = element_blank())

Class_4.5hpf_totals <- Class_4.5hpf %>% group_by(cat) %>% summarise(total = sum(sum))

Class_4.5hpf_T <- left_join(Class_4.5hpf, Class_4.5hpf_totals) %>% mutate(prop = sum / total)

ggplot(Class_4.5hpf_T, aes(x = cat, y = prop, fill = Class)) + 
  geom_col() + 
  theme_minimal() + scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = fill_colors) +
  theme(text = element_text(size = 8)) +
  labs(y = "Elements with preEGA H3K9me3") +
  theme(legend.position = "right") +
  theme(axis.title.x = element_blank())

Class_4.5hpf_Tnorm <- Class_4.5hpf_T %>% mutate(prop_norm = prop / sum_CN * 10000)

ggplot(Class_4.5hpf_Tnorm, aes(x = cat, y = prop_norm, fill = Class)) + 
  geom_col(position = position_dodge()) + 
  theme_minimal() + scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = fill_colors) +
  theme(text = element_text(size = 8)) +
  labs(y = "Elements with preEGA H3K9me3") +
  theme(legend.position = "top") +
  theme(axis.title.x = element_blank())

peric_ann <- ref_w_preEGAcats %>% filter(cat %in% c("Peric", "Peric&piRNA"))

peric_tot <- sum(peric_ann$count)

peric_ann_Class <- peric_ann %>% group_by(Class) %>% summarise(sum = sum(count), prop = sum / peric_tot)
peric_Class <- left_join(peric_ann_Class, Class_4.5hpf) %>% mutate(prop2 = sum / sum_CN)

ggplot(peric_ann_Class, aes(x = Class, y = prop, fill = Class)) + 
  geom_col(position = position_dodge()) + 
  theme_minimal() + scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = fill_colors) +
  theme(text = element_text(size = 8)) +
  labs(y = "Elements with preEGA H3K9me3") +
  theme(legend.position = "none") +
  theme(axis.title.x = element_blank())

peric_ann_Family <- peric_ann %>% group_by(Family) %>% summarise(sum = sum(count), prop = sum / peric_tot)

peric_piRNA_only <- ref_w_preEGAcats %>% filter(cat == "Peric&piRNA")
peric_piRNA_Class <- peric_piRNA_only %>% group_by(Class) %>% summarise(sum = sum(count), prop = sum / peric_tot)


####peak numbers instead of TE numbers####
peak_ann <- data.frame(cat = c("neither", "neither", "neither", "neither", "neither", "neither", "peric", "peric", "peric", "peric", "peric", "peric", 
                               "both", "both", "both", "both", "both", "both", "piRNA", "piRNA", "piRNA", "piRNA", "piRNA", "piRNA"), 
                       count = c(19630, 3255, 5304, 2309, 0, 31748, 631, 110, 313, 35, 199, 1563, 854, 425, 1276, 
                                 74, 0, 1474, 3933, 1765, 4940, 543, 0, 6533), 
                       repeat_type = c("DNA", "LINE", "LTR", "SINE", "SAT1", "total", "DNA", "LINE", "LTR", "SINE", "SAT1", "total",
                                       "DNA", "LINE", "LTR", "SINE", "SAT1", "total", "DNA", "LINE", "LTR", "SINE", "SAT1", "total"))

peak_ann_tots <- peak_ann %>% filter(repeat_type == "total") %>% rename("total" = "count") %>% select(!repeat_type)

peak_ann2 <- left_join(peak_ann, peak_ann_tots)

peak_annprops <- peak_ann2 %>% mutate(prop = count / total) %>% filter(repeat_type != "total")

ggplot(peak_annprops, aes(x = cat, y = prop, fill = repeat_type)) + 
  geom_col(position = position_dodge()) +
  theme_minimal() + scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = fill_colors) +
  theme(text = element_text(size = 8)) +
  labs(y = "Elements with preEGA H3K9me3") +
  theme(legend.position = "top") +
  theme(axis.title.x = element_blank())

peric_total <- peak_ann2 %>% filter(cat %in% c("both", "peric")) %>% filter(repeat_type != "total") %>%
  filter(repeat_type != "SAT1") %>% group_by(repeat_type) %>% summarise(rep_count = sum(count), total_count = sum(total))

peric_prop <- peric_total %>% mutate(total_rep = sum(rep_count), prop = rep_count / total_rep)

peric_peak_overlap <- ggplot(peric_prop, aes(x = repeat_type, y = prop, fill = repeat_type)) + 
  geom_col(position = position_dodge()) +
  theme_minimal() + scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = fill_colors) +
  theme(text = element_text(size = 8)) + coord_flip() +
  labs(y = "pre-EGA Pericentromeric Peak Overlap") +
  theme(legend.position = "none") +
  theme(axis.title.y = element_blank())
ggsave("~/peric_peak_overlap.pdf", plot = peric_peak_overlap, width = 1.8, height = 2, units = "in")

####import K9 signal data####
setwd("~/K9_peric_enr")

# Get the list of files in the directory
files <- list.files(pattern = "\\.tab$")

# Iterate through each file
for (file in files) {
  # Extract relevant information using a regular expression
  sat <- gsub(".*K9peaks_(.*)_peakmatrix.*", "\\1", file)
  
  # Read the file into a data frame
  df <- read_delim(file, delim = '\t', skip = 2)
  
  # Name the data frame based on the specified format (excluding "AVG" prefix)
  time <- as.numeric(str_extract(file, "\\d+\\.?\\d*"))
  df_name <- paste(sat, time, "K9", sep = "_")
  assign(df_name, df)
}

periC_only_3hpf <- peric_only_3hpf_9_K9 %>% rowwise() %>% mutate(total = rowSums(across(where(is.numeric)))) %>%
  mutate(time_point = "3hpf")
periC_only_4.5hpf <- peric_only_4.5hpf_9_K9 %>% rowwise() %>% mutate(total = rowSums(across(where(is.numeric)))) %>%
  mutate(time_point = "4.5hpf")


piRNA_only_3hpf <- piRNA_only_3hpf_9_K9 %>% rowwise() %>% mutate(total = rowSums(across(where(is.numeric)))) %>%
  mutate(time_point = "3hpf")
piRNA_only_4.5hpf <- piRNA_only_4.5hpf_9_K9 %>% rowwise() %>% mutate(total = rowSums(across(where(is.numeric)))) %>%
  mutate(time_point = "4.5hpf")

piRNA_N_peric_3hpf <- piRNA_N_peric_3hpf_9_K9 %>% rowwise() %>% mutate(total = rowSums(across(where(is.numeric)))) %>%
  mutate(time_point = "3hpf")
piRNA_N_peric_4.5hpf <- piRNA_N_peric_4.5hpf_9_K9 %>% rowwise() %>% mutate(total = rowSums(across(where(is.numeric)))) %>%
  mutate(time_point = "4.5hpf")

noperic_nopiRNA_3hpf <- noperic_nopiRNA_3hpf_9_K9 %>% rowwise() %>% mutate(total = rowSums(across(where(is.numeric)))) %>%
  mutate(time_point = "3hpf")
noperic_nopiRNA_4.5hpf <- noperic_nopiRNA_4.5hpf_9_K9 %>% rowwise() %>% mutate(total = rowSums(across(where(is.numeric)))) %>%
  mutate(time_point = "4.5hpf")

postEGA_4.5hpf <- postEGApeaks_comp_4.5hpf_peakmatrix.tab_4.5_K9 %>% rowwise() %>% mutate(total = rowSums(across(where(is.numeric)))) %>%
  mutate(time_point = "4.5hpf")

###3 hpf figure####
periC_only_3hpf <- periC_only_3hpf %>% mutate(cat = "Peric.")
piRNA_only_3hpf <- piRNA_only_3hpf %>% mutate(cat = "piRNA")
piRNA_N_peric_3hpf <- piRNA_N_peric_3hpf %>% mutate(cat = "Peric. & piRNA")
noperic_nopiRNA_3hpf <- noperic_nopiRNA_3hpf %>% mutate(cat = "Neither")

total_3h <- bind_rows(periC_only_3hpf, piRNA_only_3hpf, piRNA_N_peric_3hpf, noperic_nopiRNA_3hpf) %>%
  mutate(cat = factor(cat, levels = c("Neither", "Peric.", "piRNA", "Peric. & piRNA")))

t.test_3h <- compare_means(total ~ cat, p.adjust.method = "bonferroni", data=total_3h, method = "t.test")
t.test_3h <- t.test_3h %>% mutate(y.position = c(320, 340, 370, 330, 360, 350)) %>% mutate(p.signif = ifelse(p < 0.05, "*", "ns"))

preEGA_3h_pic <- ggplot(total_3h, aes(x = fct_rev(cat), y = total)) +
  geom_boxplot(outlier.shape = NA, aes(fill = cat), notch = TRUE) +
  labs(y = "H3K9me3 Enrichment") + ylim(0,375) +
  theme_minimal() + facet_wrap(~time_point, strip.position = "right") +
  scale_fill_manual(values = fill_colors) +
  theme(legend.position = "none") +  
  theme(text = element_text(size = 8)) +
  theme(axis.title.y = element_blank()) + coord_flip() +
  stat_pvalue_manual(t.test_3h, label = "p.signif", size = 2, tip.length = 0)
ggsave("~/preEGA_3h_pic.pdf", plot = preEGA_3h_pic, width = 4.5, height = 1.25, units = "in")

###4.5h figures####
periC_only_4.5hpf <- periC_only_4.5hpf %>% mutate(cat = "Peric.")
piRNA_only_4.5hpf <- piRNA_only_4.5hpf %>% mutate(cat = "piRNA")
piRNA_N_peric_4.5hpf <- piRNA_N_peric_4.5hpf %>% mutate(cat = "Peric. & piRNA")
noperic_nopiRNA_4.5hpf <- noperic_nopiRNA_4.5hpf %>% mutate(cat = "Neither")
postEGA_4.5hpf <- postEGA_4.5hpf %>% mutate(cat = "post-EGA")

total_4.5h <- bind_rows(periC_only_4.5hpf, piRNA_only_4.5hpf, piRNA_N_peric_4.5hpf, noperic_nopiRNA_4.5hpf, postEGA_4.5hpf) %>%
mutate(cat = factor(cat, levels = c("Neither", "Peric.", "piRNA", "Peric. & piRNA", "post-EGA"))) %>% 
  mutate(time = ifelse(cat == "post-EGA", "post-EGA", "pre-EGA"))

t.test_4.5h <- compare_means(total ~ time, p.adjust.method = "bonferroni", data=total_4.5h, method = "t.test")
t.test_4.5h <- t.test_4.5h %>% mutate(y.position = 15000) %>% mutate(p.signif = ifelse(p < 0.05, "*", "ns"))

simple_4.5h <- ggplot(total_4.5h, aes(x = time, y = total)) +
  geom_boxplot(outlier.shape = NA, aes(fill = time), notch = TRUE) +
  labs(y = "H3K9me3 Enrichment") + ylim(0,34050) +
  theme_minimal() + facet_wrap(~time_point, strip.position = "right") +
  scale_fill_manual(values = fill_colors) +
  theme(legend.position = "none") + ylim(0,15500) + 
  theme(text = element_text(size = 8)) +
  theme(axis.title.y = element_blank()) + coord_flip() +
  stat_pvalue_manual(t.test_4.5h, label = "p.signif", size = 2, tip.length = 0)
ggsave("~/simple_4.5h.pdf", plot = simple_4.5h, width = 4.5, height = 1, units = "in")

preEGA_4.5h <- total_4.5h %>% filter(cat != "post-EGA")
t.test_preEGA4.5h <- compare_means(total ~ cat, p.adjust.method = "bonferroni", data=preEGA_4.5h, method = "t.test")
t.test_preEGA4.5h <- t.test_preEGA4.5h %>% mutate(y.position = c(25000, 27000, 30000, 26000, 29000, 28000)) %>% 
                                                    mutate(p.signif = ifelse(p < 0.05, "*", "ns"))

preEGA_4.5h_pic <- ggplot(preEGA_4.5h, aes(x = fct_rev(cat), y = total)) +
  geom_boxplot(outlier.shape = NA, aes(fill = cat), notch = TRUE) +
  labs(y = "H3K9me3 Enrichment") + ylim(0,30050) +
  theme_minimal() + facet_wrap(~time_point, strip.position = "right") +
  scale_fill_manual(values = fill_colors) +
  theme(legend.position = "none") +  
  theme(text = element_text(size = 8)) +
  theme(axis.title.y = element_blank()) + coord_flip() +
  stat_pvalue_manual(t.test_preEGA4.5h, label = "p.signif", size = 2, tip.length = 0)
ggsave("~/preEGA_4.5h.pdf", plot = preEGA_4.5h_pic, width = 4.5, height = 1.25, units = "in")
