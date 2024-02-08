library("tidyverse")
library("dplyr")
library("tidyr")
library(purrr)
library(ggpubr)
fill_colors <- c("steelblue","gray83", "skyblue3","lightcyan3", "gray60", "lightsteelblue2", "gray40", "skyblue4", "lightskyblue2", "gray83")

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

big_name_table <- left_join(ref_ann, combTC) %>% mutate(., cn_prop = count / copyNumber_0.1) %>% 
  filter(., Class %in% c("DNA", "LTR", "LINE", "SINE")) %>% filter(copyNumber_0.1 >= 10) %>% filter(isLTR == "FALSE")

##analysing annotation####

K9_4.5h_only <- big_name_table %>% filter(timepoint == "4.5 hpf")
Class_4.5hpf <- K9_4.5h_only %>% group_by(Class) %>% summarise(sum = sum(count), sum_CN = sum(copyNumber_0.1))
Class_4.5hpf_1 <- Class_4.5hpf %>% mutate(prop = sum / sum_CN, prop_type = "CN")
Class_4.5hpf_2 <- Class_4.5hpf %>% mutate(prop = sum / sum(sum), prop_type = "TE")
Class_4.5hpf <- bind_rows(Class_4.5hpf_2, Class_4.5hpf_1) %>% mutate(perc = prop * 100)

####class-level figures####

TE_Props_graph <- ggplot(filter(Class_4.5hpf, prop_type == "TE"), aes(x = "", y = prop, fill = Class)) + 
  geom_bar(stat = "identity", width = 1, color = "white") + coord_polar("y", start=0) +
  theme_void() + 
  scale_fill_manual(values = fill_colors) +
  theme(text = element_text(size = 8)) +
  labs(y = "Prop. of H3K9me3-Enriched TEs") +
  theme(axis.title.x = element_blank())
ggsave("/Users/mangofrog7/Documents/Graduate_School/Papers/K9_Kinetics/Fig2/TEprops_4.5hpf_pie.pdf", plot = TE_Props_graph, width = 2, height = 2, units = "in")

CN_Props_graph <- ggplot(filter(Class_4.5hpf, prop_type == "CN"), aes(x = Class, y = prop, fill = Class)) + 
  geom_col(position = position_dodge()) + 
  theme_minimal() + scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = fill_colors) +
  theme(text = element_text(size = 8)) +
  labs(y = "Elements Enriched for H3K9me3") +
  theme(legend.position = "none") +
  theme(axis.title.x = element_blank())
ggsave("/Users/mangofrog7/Documents/Graduate_School/Papers/K9_Kinetics/Fig2/CNprops_4.5hpf.pdf", plot = CN_Props_graph, width = 1.5, height = 2, units = "in")


######family level####
Family_4.5hpf <- K9_4.5h_only %>% group_by(Class, Family, Name) %>% summarise(sum = sum(count), sum_CN = sum(copyNumber_0.1))
Family_4.5hpf <- Family_4.5hpf %>% mutate(prop = sum / sum_CN)  
Family_20p_cutoff <- Family_4.5hpf %>% group_by(Class) %>% 
  summarise(yes = sum(prop >= 0.2), no = sum(prop < 0.2)) %>%
  mutate(total = yes + no, perc = yes /total)

Family_20p_cutoff_graph <- ggplot(Family_20p_cutoff, aes(x = Class, y = perc, fill = Class)) + 
  geom_col(position = position_dodge()) + 
  theme_minimal() + scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = fill_colors) +
  theme(text = element_text(size = 8)) +
  labs(y = "Families w/ 20% Elements Enriched for H3K9me3") +
  theme(legend.position = "none") +
  theme(axis.title.x = element_blank())
ggsave("/Users/mangofrog7/Documents/Graduate_School/Papers/K9_Kinetics/Fig2/Family_20p_cutoff.pdf", plot = Family_20p_cutoff_graph, width = 1.5, height = 2, units = "in")
ggsave("/Users/mangofrog7/Documents/Graduate_School/Papers/K9_Kinetics/Fig2/Family_20p_cutoff_bigger.pdf", plot = Family_20p_cutoff_graph, width = 1.5, height = 3, units = "in")

####LTR families####
LTRs_20p_Fams <- Family_4.5hpf %>% filter(Class == "LTR", Family %in% c("Copia", "DIRS", "ERV1", "Gypsy", "Ngaro", "Pao")) %>% 
  group_by(Family) %>% summarise(yes = sum(prop >= 0.2), no = sum(prop < 0.2)) %>%
  mutate(total = yes + no, perc = yes /total)
LTRs_20props <- Family_4.5hpf %>% filter(Class == "LTR", Family %in% c("Copia", "DIRS", "ERV1", "Gypsy", "Ngaro", "Pao")) %>% 
  group_by(Family) %>% summarise(count = sum(sum), CN = sum(sum_CN)) %>% mutate(prop = count / CN) %>% mutate(perc = prop * 100)

LTRs_20props_pics <- ggplot(LTRs_20props, aes(x = Family, y = prop, fill = Family)) + 
  geom_col(position = position_dodge()) + 
  theme_minimal() + scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = fill_colors) +
  theme(text = element_text(size = 8)) +
  labs(y = "Elements Enriched for H3K9me3") +
  theme(legend.position = "none") +
  theme(axis.title.x = element_blank())
ggsave("/Users/mangofrog7/Documents/Graduate_School/Papers/K9_Kinetics/Fig2/LTRs_20props_pics.pdf", plot = LTRs_20props_pics, width = 2.5, height = 2, units = "in")

LTRs_20p_cutoff <- ggplot(LTRs_20p_Fams, aes(x = Family, y = perc, fill = Family)) + 
  geom_col(position = position_dodge()) + 
  theme_minimal() + scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = fill_colors) +
  theme(text = element_text(size = 8)) +
  labs(y = "Subfamilies w/ 20% Elements Enriched for H3K9me3") +
  theme(legend.position = "none") +
  theme(axis.title.x = element_blank())
ggsave("/Users/mangofrog7/Documents/Graduate_School/Papers/K9_Kinetics/Fig2/LTRs_20p_cutoff.pdf", plot = LTRs_20p_cutoff, width = 2.5, height = 2, units = "in")
ggsave("/Users/mangofrog7/Documents/Graduate_School/Papers/K9_Kinetics/Fig2/LTRs_20p_cutoff_tall.pdf", plot = LTRs_20p_cutoff, width = 2.5, height = 3, units = "in")

####is proportion of copies enriched related to TE age####
K9_4.5h_only <- K9_4.5h_only %>% mutate(perc = cn_prop * 100)

cnprop_medlen_scatter <- ggplot(K9_4.5h_only, aes(x = perc, y = MedLen)) + 
  geom_point(size = 0.5) + geom_smooth(method = "lm", se=F, linewidth = 0.5) +
  facet_wrap(vars(Class), nrow = 1) +  
  theme_minimal() +  theme(text = element_text(size = 8)) + stat_cor(size = 2) +
  labs(y = "Median Branch Length", x = "Percent of Elements Enriched for H3K9me3") +
  theme(legend.position = "none") + ylim(0,0.5)
ggsave("/Users/mangofrog7/Documents/Graduate_School/Papers/K9_Kinetics/Fig2/cnprop_medlen_scatter.pdf", plot = cnprop_medlen_scatter, width = 5.5, height = 2, units = "in")

####now just compare 3 hpf to 4.5 hpf####
library(ggprism)

col_names <- c("count", "Name")
K9_4.5m3h <- read_delim("/Users/mangofrog7/Documents/Graduate_School/Projects/bioinformatics/TC_final/peaks/ann2/4.5m3h_TEann_counts.bed", col_names = col_names) %>%
  rename(late_hpf = count) %>% select(Name, late_hpf)
K9_3h_only <- big_name_table %>% filter(timepoint == "3 hpf") %>% rename(early_hpf = count) %>% select(Name, Family, Class, early_hpf)

K9_4.5m3h_table <- left_join(ref_ann, K9_4.5m3h) %>% filter(isLTR == "FALSE")
EGA_comp <- left_join(K9_4.5m3h_table, K9_3h_only)
EGA_comp <- EGA_comp %>% mutate(EGA_frac = early_hpf / late_hpf) %>% filter(Class %in% c("DNA", "SINE", "LINE", "LTR"))

class_values <- unique(EGA_comp$Class)

t_test_results_list <- list()

for (class_value in unique(EGA_comp$Class)) {
  subset_data <- subset(EGA_comp, Class == class_value)
  t_test_result <- t.test(subset_data$EGA_frac, mu = 1)  # Adjust mu as needed
  t_test_results_list[[as.character(class_value)]] <- t_test_result
}

annotations_df <- data.frame(
  group1 = rep(unique(EGA_comp$Class), each = 1), 
  group2 = NA,
  p = sapply(t_test_results_list, function(x) x$p.value),
  y.position = c(1, 1, 1, 1))

df_p_val <- data.frame(
  group1 = c("DNA", "LTR"),
  group2 = c(1, 3),
  x = c(1, 3),
  label = signif(annotations_df$p[annotations_df$p < 0.05], digits = 2),
  y.position = c(12, 12)
)

EGA_comp_trimmed <- EGA_comp %>% 
  group_by(Class) %>% 
  mutate(zPW = scale(EGA_frac)) %>% 
  filter(between(zPW,-5,+5))

EGA_comp_class <- ggplot(EGA_comp_trimmed, aes(x = Class, y = EGA_frac)) +
  geom_boxplot(outlier.size = 1, outlier.alpha = 0.5, aes(fill = Class)) +
  labs(y = "preEGA/postEGA H3K9me3-enriched elements") +
  theme_minimal() +
  scale_fill_manual(values = fill_colors) +
  theme(legend.position = "none") +
  theme(axis.title.x = element_blank()) +
  theme(text = element_text(size = 8)) +
  add_pvalue(df_p_val, 
             xmin = "group1", 
             x = "x", 
             label = "p = {label}",
             y.position = "y.position", 
             label.size = 2.5) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey20", linewidth = 0.5)
ggsave("/Users/mangofrog7/Documents/Graduate_School/Papers/K9_Kinetics/Fig4/EGA_comp_class.pdf", plot = EGA_comp_class, width = 2.35, height = 2.75, units = "in")

####LTR families####
EGA_comp_LTRs <- EGA_comp%>% filter(Class == "LTR") %>% filter(Family %in% c("Copia", "DIRS", "ERV1", "Gypsy", "Ngaro", "Pao"))

family_values <- unique(EGA_comp_LTRs$Family)

family_results_list <- list()

for (family_value in family_values) {
  subset_data <- subset(EGA_comp_LTRs, Family == family_value)
  t_test_result <- t.test(subset_data$EGA_frac, mu = 1)  # Adjust mu as needed
  family_results_list[[as.character(family_value)]] <- t_test_result
}

ann_df <- data.frame(
  group1 = rep(unique(EGA_comp_LTRs$Family), each = 1), 
  group2 = NA,
  p = sapply(family_results_list, function(x) x$p.value),
  y.position = c(1, 1, 1, 1, 1, 1))

df_pval <- data.frame(
  group1 = c("Gypsy", "Pao"),
  group2 = c(4, 6),
  x = c(4, 6),
  label = signif(ann_df$p[ann_df$p < 0.05], digits = 2),
  y.position = c(16, 16)
)

EGA_LTRs_trimmed <- EGA_comp_LTRs %>% 
  group_by(Family) %>% 
  mutate(zPW = scale(EGA_frac)) %>% 
  filter(between(zPW,-5,+5))

EGA_comp_LTR <- ggplot(EGA_LTRs_trimmed, aes(x = Family, y = EGA_frac)) + 
  geom_boxplot(outlier.size = 1, outlier.alpha = 0.5, aes(fill = Family)) + 
  labs(y = "preEGA/postEGA H3K9me3-enriched elements")+ 
  theme_minimal() + scale_fill_manual(values = fill_colors) +
  theme(legend.position="none") + 
  theme(axis.title.x = element_blank()) + 
  theme(text = element_text(size = 8)) + 
  add_pvalue(df_pval, 
             xmin = "group1", 
             x = "x", 
             label = "p = {label}",
             y.position = "y.position", 
             label.size = 2.5) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey20", linewidth = 0.5)
ggsave("/Users/mangofrog7/Documents/Graduate_School/Papers/K9_Kinetics/Fig4/EGA_comp_LTR.pdf", plot = EGA_comp_LTR, width = 3.1, height = 2.7, units = "in")

EGA_comp_LTR_split <- EGA_comp_LTRs %>% mutate(cat = ifelse(EGA_frac > 0.5, "pre", "post")) %>% 
  mutate(famCat = ifelse(Family %in% c("Gypsy", "Pao"), "pre", ifelse(Family %in% c("Ngaro", "DIRS"), "post", "NA")))

####import piRNA data####
piRNA_Table <- read.csv("/Users/mangofrog7/Documents/Graduate_School/Projects/bioinformatics/piRNA/piRNA_TPM.csv")
piRNAs_w3hK9 <- read_delim("/Users/mangofrog7/Documents/Graduate_School/Projects/bioinformatics/piRNA/piRNA_3hK9_TEs_counts.bed", col_names = col_names)

piRNA_big <- left_join(piRNA_Table, piRNAs_w3hK9)
LTR_piRNA <- left_join(ref_ann, piRNA_big) %>% 
  filter(Class == "LTR", isLTR == "FALSE", Family %in% c("Copia", "DIRS", "ERV1", "Gypsy", "Ngaro", "Pao"))  %>% 
  rename(copies_w3h_K9 = count)
LTR_piRNA[is.na(LTR_piRNA)] <- 0

LTR_K9_piRNA <- left_join(EGA_comp_LTRs, LTR_piRNA) %>% 
  mutate(K9_w_piRNA = copies_w3h_K9 / early_hpf, K9_w_piRNA = if_else(K9_w_piRNA > 1, 1, K9_w_piRNA)) %>% 
  filter(copyNumber_0.1 >= 10) %>% filter(early_hpf >= 10)

LTR_fam_piRNA_ttest <- compare_means(K9_w_piRNA ~ Family, data = LTR_K9_piRNA, method = "t.test", p.adjust.method = "bonferroni")
LTR_fam_piRNA_ttest<- LTR_fam_piRNA_ttest %>% mutate(y.position = c(0, 0, 0, 1.26, 0, 1.05, 1.12, 1.12, 1.19, 0))


fill_colors2 <- c("gray83", "skyblue3","lightcyan3", "gray60", "lightsteelblue2", "gray40", "skyblue4", "lightskyblue2", "gray83")

LTR_K9_piRNA_plot <- ggplot(LTR_K9_piRNA, aes(x = Family, y = K9_w_piRNA)) + 
  geom_boxplot(outlier.size = 1, outlier.alpha = 0.5, aes(fill = Family)) + 
  labs(y = "% of 3 hpf H3K9me3 Peaks targeted by piRNAs") + 
  theme_minimal() + scale_fill_manual(values = fill_colors2) +
  theme(legend.position="none") + scale_y_continuous(labels = scales::percent, breaks=c(0.2, 0.4, 0.6, 0.8, 1)) + 
  theme(text = element_text(size = 8)) +  theme(axis.title.x = element_blank()) + 
  stat_pvalue_manual(LTR_fam_piRNA_ttest, label = "p = {p.adj}", hide.ns = TRUE, size = 2.5, tip.length = 0.01)
ggsave("/Users/mangofrog7/Documents/Graduate_School/Papers/K9_Kinetics/Fig4/LTR_K9_piRNA_plot.pdf", plot = LTR_K9_piRNA_plot, width = 3.2, height = 2.7, units = "in")


####now going to make an upset plot of peaks across development####
library(ComplexUpset)

intersection_data <- c("two" = 964, "two.five" = 2049, "three" = 3580, "three.five" = 3212, "four" = 5469, "four.five" = 9078, 
                   "two&two.five" = 148, "two.five&three" = 348, "three&three.five" = 620, "three.five&four" = 2354, "four&four.five" = 6415, 
                   "two&two.five&three" = 71, "two.five&three&three.five" = 208, "three&three.five&four" = 981, "three.five&four&four.five" = 8583, 
                   "two&two.five&three&three.five" = 62, "two.five&three&three.five&four" = 402, "three&three.five&four&four.five" = 6891, 
                   "two&two.five&three&three.five&four" = 99, "two.five&three&three.five&four&four.five" = 6701, "two&two.five&three&three.five&four&four.five" = 3739)

int_data <- fromExpression(intersection_data)
colnames(int_data) <- c("2 hpf", "2.5 hpf", "3 hpf", "3.5 hpf", "4 hpf", "4.5 hpf")

time_points = colnames(int_data)

upset(int_data, time_points, name = "Timepoint")


upset <- upset(
  int_data, time_points,
  base_annotations = list(
    'Intersection size'=(
      intersection_size(counts = FALSE,
        mapping=aes(fill='bars_color'),
        text_mapping=aes(
          label=paste0(round(!!get_size_mode('exclusive_intersection')/nrow(int_data) * 100, 1), ''), 
      )
      )
      + ylab('% of Peaks')
      + scale_y_continuous(
        labels=scales::percent_format(scale=100 / nrow(int_data)),
        breaks=c(0, 5, 10, 15, 20) / 100 * nrow(int_data)
      ) + scale_fill_manual(values=c('bars_color'='steelblue'), guide='none')
      )
  ),
  set_sizes=(FALSE), 
  sort_sets=FALSE,
  name = "Timepoint", 
  matrix=(
    intersection_matrix(
      geom=geom_point(size=1.5
      ),
      outline_color=list(
        active='gray40',
        inactive='gray83'
      ), segment=geom_segment(
        color = "gray40"
      )
    )
    + scale_color_manual(
      values=c('TRUE'='gray40', 'FALSE'='gray83'),
      breaks=c('TRUE', 'FALSE'), 
      guide = "none")) + 
    theme(text = element_text(size = 8)))

ggsave("/Users/mangofrog7/Documents/Graduate_School/Papers/K9_Kinetics/Fig3/upset.pdf", plot = upset, width = 4.5, height = 3, units = "in")



####pre/postEGA peak int w peric####
##imp data##
preEGA_NULL <- read_delim("/Users/mangofrog7/Documents/Graduate_School/Projects/bioinformatics/TC_final/pericentromere/preEGA_peaks_null_ints.tsv", 
                        col_names = FALSE, delim = "/", trim_ws = TRUE) %>% distinct() %>% filter(X1 != 0) %>% select(X1, X6) %>% rename(count = X1, loc = X6)
postEGA_NULL <- read_delim("/Users/mangofrog7/Documents/Graduate_School/Projects/bioinformatics/TC_final/pericentromere/postEGA_peaks_null_ints.tsv", 
                         col_names = FALSE, delim = "/", trim_ws = TRUE) %>% distinct() %>% filter(X1 != 0) %>% select(X1, X6) %>% rename(count = X1, loc = X6)

preEGA_perc <- read_delim("/Users/mangofrog7/Documents/Graduate_School/Projects/bioinformatics/TC_final/pericentromere/preEGA_peaks.peric.bed")
postEGA_perc <- read_delim("/Users/mangofrog7/Documents/Graduate_School/Projects/bioinformatics/TC_final/pericentromere/postEGApeaks_comp.peric.bed")

##clean and combine##
preEGA_sum <- preEGA_perc %>% summarise(peric_count = n()) %>% mutate(ega = "pre")
postEGA_sum <- postEGA_perc %>% summarise(peric_count = n()) %>% mutate(ega = "post")

preEGA_NULL_neat <- data.frame(t(select(preEGA_NULL, count))) %>% mutate_all(as.numeric) %>% mutate(ega = "pre")
postEGA_NULL_neat <- data.frame(t(select(postEGA_NULL, count))) %>% mutate_all(as.numeric) %>% mutate(ega = "post")

preEGA_comb <- left_join(preEGA_sum, preEGA_NULL_neat)
postEGA_comb <- left_join(postEGA_sum, postEGA_NULL_neat)

total <- bind_rows(preEGA_comb, postEGA_comb) %>% select(-which(colSums(is.na(.)) > 0))

##permutation test##
peric_EGA_greater <- total %>%
  mutate(across(.cols = -ega, .fns = list(Greater = ~as.numeric(peric_count) > as.numeric(.)), .names = "Is_{.col}")) %>%
  select(ega, starts_with("Is_"), -Is_peric_count) %>%
  mutate_at(vars(starts_with("Is_")), ~ as.integer(.))

peric_EGA_less <- total %>%
  mutate(across(.cols = -ega, .fns = list(Less = ~as.numeric(peric_count) < as.numeric(.)), .names = "Is_{.col}")) %>%
  select(ega, starts_with("Is_"), -Is_peric_count) %>%
  mutate_at(vars(starts_with("Is_")), ~ as.integer(.))

peric_EGA_permTest <- left_join((rowwise(peric_EGA_greater) %>% 
                           mutate(gRowSum = sum(c_across(starts_with("Is_")))) %>% 
                           select( -starts_with("Is_"))), (rowwise(peric_EGA_less) %>% 
                                                             mutate(lRowSum = sum(c_across(starts_with("Is_")))) %>% 
                                                             select( -starts_with("Is_")))) %>% 
                            mutate(gPvalue = (gRowSum + 1) / 1001, lPvalue = (lRowSum +1) / 1001) %>% 
                            mutate(ega = ifelse(ega == "pre", "preEGA", ifelse(ega == "post", "postEGA")))

##calc means and sds to plot##
preEGA_NULL <- preEGA_NULL %>% mutate(count = as.numeric(count)) %>% mutate(ega = "preEGA") %>% filter(!is.na(count))
postEGA_NULL <- postEGA_NULL %>% mutate(count = as.numeric(count)) %>% mutate(ega = "postEGA") %>% filter(!is.na(count))

comb_null <- bind_rows(preEGA_NULL, postEGA_NULL)
null_summ <- comb_null %>% group_by(ega, loc) %>% summarise(sd = sd(count), count = mean(count)) %>% mutate(loc = "Control")

preEGA_peric <- preEGA_perc %>% summarise(count = n()) %>% mutate(ega = "preEGA", loc = "Pericentromere", sd = 0)
postEGA_peric <- postEGA_perc %>% summarise(count = n()) %>% mutate(ega = "postEGA", loc = "Pericentromere", sd = 0)
comb_peric <- bind_rows(preEGA_peric, postEGA_peric)

for_plot <- bind_rows(null_summ, comb_peric) %>% 
  mutate(loc = factor(loc, levels = c("Pericentromere", "Control"))) %>% mutate(ega = factor(ega, levels = c("preEGA", "postEGA")))

##format perm test for stat_pvalue function##
perm_test <- data.frame(peric_EGA_permTest %>% select(ega) %>% 
                          mutate(group2 = "NA") %>% left_join(., select(peric_EGA_permTest, ega, lPvalue)) %>% rename (group1 = ega))
perm_test <- perm_test %>% mutate(y.position = c(2500, 1600)) %>% 
  mutate(group2 = c("preEGA", "postEGA")) %>% mutate(lPvalue = signif(perm_test$lPvalue, digits = 2))

peric_permtest_plot <- ggplot(for_plot, aes(x=ega, y=count, fill=loc)) +
  geom_col(position = position_dodge()) +
  geom_errorbar(aes(ymin=count-sd, ymax=count+sd), position = position_dodge(0.9), width = 0.2) + 
  scale_fill_manual(values = fill_colors) +
  theme_minimal() + labs(y = "Number of Peaks") +
  theme(legend.position = "none") +
  theme(legend.title = element_blank()) +
  theme(axis.title.x = element_blank()) +
  theme(text = element_text(size = 8)) +
  guides(fill = guide_legend(override.aes = list(size = 0.05))) + 
  stat_pvalue_manual(perm_test, label = "p = {lPvalue}", size = 2.5, inherit.aes = FALSE)
ggsave("/Users/mangofrog7/Documents/Graduate_School/Papers/K9_Kinetics/Fig5/peric_permtest_plot.pdf", plot = peric_permtest_plot, width = 1.75, height = 2.5, units = "in")


  