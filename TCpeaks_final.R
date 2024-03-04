library("tidyverse")
library("dplyr")
library("tidyr")
library(purrr)
library(ggpubr)
fill_colors <- c("steelblue","gray83", "skyblue3","lightcyan3", "gray60", "lightsteelblue2", "gray40", "skyblue4", "lightskyblue2", "gray83")

####import annotation data####
col_names <- c("count", "Name")

K9_2hpf <- read_delim("~/peaks/ann2/2hpf_K9_.TEann_35_0.1filt_counts.bed", col_names = col_names) %>% 
 mutate(timepoint = "2 hpf")
K9_2.5hpf <- read_delim("~/peaks/ann2/2.5hpf_K9_.TEann_35_0.1filt_counts.bed", col_names = col_names) %>% 
  mutate(timepoint = "2.5 hpf")
K9_3hpf <- read_delim("~/peaks/ann2/3hpf_K9_.TEann_35_0.1filt_counts.bed", col_names = col_names) %>% 
  mutate(timepoint = "3 hpf")
K9_3.5hpf <- read_delim("~/peaks/ann2/3.5hpf_K9_.TEann_35_0.1filt_counts.bed", col_names = col_names) %>% 
  mutate(timepoint = "3.5 hpf")
K9_4hpf <- read_delim("~/peaks/ann2/4hpf_K9_.TEann_35_0.1filt_counts.bed", col_names = col_names) %>% 
  mutate(timepoint = "4 hpf")
K9_4.5hpf <- read_delim("~/peaks/ann2/4.5hpf_K9_.TEann_35_0.1filt_counts.bed", col_names = col_names) %>% 
  mutate(timepoint = "4.5 hpf")


combTC <- bind_rows(K9_2.5hpf, K9_2hpf, K9_3.5hpf, K9_3hpf, K9_4.5hpf, K9_4hpf)

##importing TE table####
ref_ann <- read_delim("~/general/TEref_table.txt")  
TE_counts_0.1filt <- read_delim("~/general/TEann_35_0.1filt_counts.txt", col_names = c("copyNumber_0.1", "Name"))  

ref_ann <- left_join(ref_ann, TE_counts_0.1filt)

big_name_table <- left_join(ref_ann, combTC) %>% mutate(., cn_prop = count / copyNumber_0.1) %>% 
  filter(., Class %in% c("DNA", "LTR", "LINE", "SINE")) %>% filter(copyNumber_0.1 >= 10) %>% filter(isLTR == "FALSE")

###making supp table 2####
int_df <- big_name_table %>% select(Name, Family, Class,copyNumber_0.1, count, timepoint)

table2 <- pivot_wider(int_df, 
            id_cols = c(Name, Family, Class,copyNumber_0.1), 
            names_from = "timepoint", 
            values_from = "count") %>% select(!"NA")
table2[is.na(table2)] <- 0

Table_S2 <- table2 %>% select(Name, Family, Class,copyNumber_0.1, "2 hpf", "2.5 hpf", "3 hpf", "3.5 hpf", "4 hpf", "4.5 hpf")

write.csv(Table_S2, "~/Table_S2.csv")


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
ggsave("~/TEprops_4.5hpf_pie.pdf", plot = TE_Props_graph, width = 2, height = 2, units = "in")

CN_Props_graph <- ggplot(filter(Class_4.5hpf, prop_type == "CN"), aes(x = Class, y = prop, fill = Class)) + 
  geom_col(position = position_dodge()) + 
  theme_minimal() + scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = fill_colors) +
  theme(text = element_text(size = 8)) +
  labs(y = "Elements Enriched for H3K9me3") +
  theme(legend.position = "none") +
  theme(axis.title.x = element_blank())
ggsave("~/CNprops_4.5hpf.pdf", plot = CN_Props_graph, width = 1.5, height = 2, units = "in")

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
ggsave("~/Family_20p_cutoff.pdf", plot = Family_20p_cutoff_graph, width = 1.5, height = 2, units = "in")
ggsave("~/Family_20p_cutoff_bigger.pdf", plot = Family_20p_cutoff_graph, width = 1.5, height = 3, units = "in")

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
ggsave("~/LTRs_20props_pics.pdf", plot = LTRs_20props_pics, width = 2.5, height = 2, units = "in")

LTRs_20p_cutoff <- ggplot(LTRs_20p_Fams, aes(x = Family, y = perc, fill = Family)) + 
  geom_col(position = position_dodge()) + 
  theme_minimal() + scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = fill_colors) +
  theme(text = element_text(size = 8)) +
  labs(y = "Subfamilies w/ 20% Elements Enriched for H3K9me3") +
  theme(legend.position = "none") +
  theme(axis.title.x = element_blank())
ggsave("~/LTRs_20p_cutoff.pdf", plot = LTRs_20p_cutoff, width = 2.5, height = 2, units = "in")
ggsave("~/LTRs_20p_cutoff_tall.pdf", plot = LTRs_20p_cutoff, width = 2.5, height = 3, units = "in")

####is proportion of copies enriched related to TE age####
K9_4.5h_only <- K9_4.5h_only %>% mutate(perc = cn_prop * 100)

cnprop_medlen_scatter <- ggplot(K9_4.5h_only, aes(x = perc, y = MedLen)) + 
  geom_point(size = 0.5) + geom_smooth(method = "lm", se=F, linewidth = 0.5) +
  facet_wrap(vars(Class), nrow = 1) +  
  theme_minimal() +  theme(text = element_text(size = 8)) + stat_cor(size = 2) +
  labs(y = "Median Branch Length", x = "Percent of Elements Enriched for H3K9me3") +
  theme(legend.position = "none") + ylim(0,0.5)
ggsave("~/cnprop_medlen_scatter.pdf", plot = cnprop_medlen_scatter, width = 5.5, height = 2, units = "in")

####now just compare 3 hpf to 4.5 hpf for pre/post EGA analysis####
col_names1 <- c("late_hpf", "Name")
K9_4.5m3h <- read_delim("~/peaks/ann2/4.5m3h_TEann_counts.bed", col_names = col_names1) 
K9_3h_only <- big_name_table %>% filter(timepoint == "3 hpf") %>% select(Name, Family, Class, count)

K9_4.5m3h_table <- left_join(ref_ann, K9_4.5m3h) %>% filter(isLTR == "FALSE")
EGA_comp <- left_join(K9_4.5m3h_table, K9_3h_only)
EGA_comp <- EGA_comp %>% mutate(EGA_frac = count / late_hpf) %>% filter(Class %in% c("DNA", "SINE", "LINE", "LTR")) %>% filter(copyNumber_0.1 > 10)

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
  y.position = c(3.8, 3.8)
)

EGA_comp_class_noOut <- ggplot(EGA_comp, aes(x = Class, y = EGA_frac)) +
  geom_boxplot(outlier.shape = NA, aes(fill = Class)) +
  labs(y = "preEGA/postEGA H3K9me3-enriched elements") +
  theme_minimal() + ylim(0,4) +
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
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey20", linewidth = 0.2)
ggsave("~/EGA_comp_class_noOut.pdf", plot = EGA_comp_class_noOut, width = 2.35, height = 2.75, units = "in")

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
  y.position = c(6.5, 6.5)
)

EGA_comp_LTR_NOout <- ggplot(EGA_comp_LTRs, aes(x = Family, y = EGA_frac)) + 
  geom_boxplot(outlier.shape = NA, aes(fill = Family)) + 
  labs(y = "preEGA/postEGA H3K9me3-enriched elements")+ 
  theme_minimal() + scale_fill_manual(values = fill_colors) +
  theme(legend.position="none") + ylim(0,7) + 
  theme(axis.title.x = element_blank()) + 
  theme(text = element_text(size = 8)) + 
  add_pvalue(df_pval, 
             xmin = "group1", 
             x = "x", 
             label = "p = {label}",
             y.position = "y.position", 
             label.size = 2.5) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey20", linewidth = 0.2)
ggsave("~/EGA_comp_LTR_NOout.pdf", plot = EGA_comp_LTR_NOout, width = 3, height = 2.7, units = "in")

EGA_comp_LTR_split <- EGA_comp_LTRs %>%
  mutate(famCat = ifelse(Family %in% c("Gypsy", "Pao"), "pre", ifelse(Family %in% c("Ngaro", "DIRS", "Copia", "ERV1"), "post", "NA")))

###TE age related to pre/post EGA K9 enrichment 
EAGfrac_MedLen_scatter <- ggplot(EGA_comp, aes(x=EGA_frac, y=MedLen)) +geom_point(size = 0.5) + 
  geom_smooth(method = "lm", se=F, linewidth = 0.5) + 
  facet_grid(~Class, scales = 'free') + theme_minimal() +
  theme(text = element_text(size = 8)) + stat_cor(size = 2) +
  labs(y = "Median Branch Length", x = "preEGA/postEGA H3K9me3-enriched elements")
ggsave("~/EAGfrac_MedLen_scatter.pdf", plot = EAGfrac_MedLen_scatter, width = 5.5, height = 2, units = "in")

####import piRNA data####
piRNA_Table <- read.csv("~/piRNA/piRNA_TPM.csv")
piRNAs_w3hK9 <- read_delim("~/piRNA/piRNA_3hK9_TEs_counts.bed", col_names = c("copies_w3h_K9", "Name"))

piRNA_big <- left_join(piRNA_Table, piRNAs_w3hK9) %>% left_join(., K9_3hpf)

K9_piRNA <- left_join(ref_ann, piRNA_big) %>% 
  filter(Class %in% c("LTR", "DNA", "SINE", "LINE"), isLTR == "FALSE", copyNumber_0.1 >= 10, count > 10)

K9_piRNA <- K9_piRNA %>% mutate(perc_K9 = copies_w3h_K9 / count) %>% select(!timepoint) %>% mutate(perc_K9 = ifelse(perc_K9 > 1, 1, perc_K9))
K9_piRNA[is.na(K9_piRNA)] <- 0

K9_piRNA_Class_t.test <- compare_means(perc_K9 ~ Class, data = K9_piRNA, method = "t.test", p.adjust.method = "bonferroni", ref.group = "LTR")
K9_piRNA_Class_t.test <- K9_piRNA_Class_t.test %>% mutate(y.position = c(1.13, 1.01, 1.07))

K9_piRNA_Class_perc_K9 <- ggplot(K9_piRNA, aes(x = Class, y = perc_K9)) + 
  geom_boxplot(outlier.shape = NA, aes(fill = Class)) + 
  labs(y = "pre-EGA H3K9me3-Marked Elements with piRNAs")+ 
  theme_minimal() + scale_fill_manual(values = fill_colors) +
  theme(legend.position="none") + scale_y_continuous(labels = scales::percent) +
  theme(axis.title.x = element_blank()) + 
  theme(text = element_text(size = 7)) +
  stat_pvalue_manual(K9_piRNA_Class_t.test, label = "p = {p.adj}", hide.ns = TRUE, size = 2.5, tip.length = 0)
ggsave("~/K9_piRNA_Class_perc_exp.pdf", plot = K9_piRNA_Class_perc_K9, width = 2.25, height = 2.15, units = "in")

K9_piRNA_summ <- K9_piRNA %>% group_by(Class) %>% summarise(total = sum(count), tot_wK9 = sum(copies_w3h_K9)) %>%
  mutate(prop = tot_wK9 / total)

K9_piRNA_col <- ggplot(K9_piRNA_summ, aes(x=Class, y=prop, fill = Class)) +
  geom_col(position = position_dodge()) + 
  labs(y = "H3K9me3-Marked Elements w piRNAs")+ 
  theme_minimal() + scale_fill_manual(values = fill_colors) +
  theme(legend.position="none") + scale_y_continuous(labels = scales::percent) +
  theme(axis.title.x = element_blank()) + 
  theme(text = element_text(size = 8))
ggsave("~/K9_piRNA_col.pdf", plot = K9_piRNA_col, width = 2.2, height = 2, units = "in")

LTR_piRNA <- K9_piRNA %>% 
  filter(Class == "LTR", isLTR == "FALSE", Family %in% c("Copia", "DIRS", "ERV1", "Gypsy", "Ngaro", "Pao")) %>%
  mutate(Family = factor(Family, levels = c("Copia", "DIRS", "ERV1", "Ngaro", "Gypsy", "Pao"))) %>%
  mutate(cat = ifelse(Family %in% c("Gypsy", "Pao"), "preskew", "postskew"))
LTR_piRNA[is.na(LTR_piRNA)] <- 0

LTR_piRNA_fam_t.test <- compare_means(perc_K9 ~ cat, data = LTR_piRNA, method = "t.test", p.adjust.method = "bonferroni")
LTR_piRNA_fam_t.test <- LTR_piRNA_fam_t.test %>% mutate(y.position = 1.1, group1 = "ERV1", group2 = "Pao")

LTR_piRNA_fam_graph <- ggplot(LTR_piRNA, aes(x = Family, y = perc_K9)) + 
  geom_boxplot(outlier.shape = NA, aes(fill = Family)) + 
  labs(y = "H3K9me3-Marked Elements w piRNAs")+ 
  theme_minimal() + scale_fill_manual(values = fill_colors) +
  theme(legend.position="none") + scale_y_continuous(labels = scales::percent) +
  theme(axis.title.x = element_blank()) + 
  theme(text = element_text(size = 7)) +
  stat_pvalue_manual(LTR_piRNA_fam_t.test, label = "p = {p.adj}", hide.ns = TRUE, size = 2.5, tip.length = 0) +
  geom_segment(aes(x = 1, y = 1.05, xend = 3, yend = 1.05)) +
  geom_segment(aes(x = 4, y = 1.05, xend = 5, yend = 1.05)) 
ggsave("~/LTR_piRNA_fam_graph.pdf", plot = LTR_piRNA_fam_graph, width = 2.25, height = 2.15, units = "in")

LTR_piRNA_summ <- LTR_piRNA %>% group_by(Family) %>% summarise(total = sum(count), tot_wK9 = sum(copies_w3h_K9)) %>%
  mutate(prop = tot_wK9 / total)

LTR_piRNA_col <- ggplot(LTR_piRNA_summ, aes(x=Family, y=prop, fill = Family)) +
  geom_col(position = position_dodge()) + 
  labs(y = "H3K9me3-Marked Elements w piRNAs")+ 
  theme_minimal() + scale_fill_manual(values = fill_colors) +
  theme(legend.position="none") + scale_y_continuous(labels = scales::percent) +
  theme(axis.title.x = element_blank()) + 
  theme(text = element_text(size = 7))
ggsave("~/LTR_piRNA_col.pdf", plot = LTR_piRNA_col, width = 2.2, height = 2, units = "in")

####now going to make an upset plot of peaks across development####
library(UpSetR)
library("ComplexUpset")

intersection_data <- c("two" = 964, "two.five" = 2049, "three" = 3580, "three.five" = 3212, "four" = 5469, "four.five" = 9078, 
                   "two&two.five" = 148, "two.five&three" = 348, "three&three.five" = 620, "three.five&four" = 2354, "four&four.five" = 6415, 
                   "two&two.five&three" = 71, "two.five&three&three.five" = 208, "three&three.five&four" = 981, "three.five&four&four.five" = 8583, 
                   "two&two.five&three&three.five" = 62, "two.five&three&three.five&four" = 402, "three&three.five&four&four.five" = 6891, 
                   "two&two.five&three&three.five&four" = 99, "two.five&three&three.five&four&four.five" = 6701, "two&two.five&three&three.five&four&four.five" = 3739)

int_data <- fromExpression(intersection_data)
colnames(int_data) <- c("2 hpf", "2.5 hpf", "3 hpf", "3.5 hpf", "4 hpf", "4.5 hpf")

time_points = colnames(int_data)

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




ggsave("~/upset.pdf", plot = upset, width = 4.5, height = 3, units = "in")


####pre/postEGA peak int w peric & piRNA####
##imp peric data##
preEGA_pericNULL <- read_delim("~/pericentromere/preEGA_peaks_null_ints.tsv", 
                        col_names = FALSE, delim = "/", trim_ws = TRUE) %>% distinct() %>% filter(X1 != 0) %>% select(X1, X6) %>% rename(count = X1, loc = X6)
postEGA_pericNULL <- read_delim("~/pericentromere/postEGA_peaks_null_ints.tsv", 
                         col_names = FALSE, delim = "/", trim_ws = TRUE) %>% distinct() %>% filter(X1 != 0) %>% select(X1, X6) %>% rename(count = X1, loc = X6)

preEGA_perc <- read_delim("~/pericentromere/preEGA_peaks.peric.bed")
postEGA_perc <- read_delim("~/pericentromere/postEGApeaks_comp.peric.bed")

##clean and combine##
preEGA_peric_sum <- preEGA_perc %>% summarise(peric_count = n(), peric_norm = peric_count / 38619) %>% mutate(ega = "pre")
postEGA_peric_sum <- postEGA_perc %>% summarise(peric_count = n(), peric_norm = peric_count / 36200) %>% mutate(ega = "post")

preEGA_pericNULL <- preEGA_pericNULL %>% mutate(count = as.numeric(count)) %>% mutate(ega = "preEGA") %>% mutate(null_norm = count / 38619) %>% filter(!is.na(count))
postEGA_pericNULL <- postEGA_pericNULL %>% mutate(count = as.numeric(count)) %>% mutate(ega = "postEGA") %>% mutate(null_norm = count / 36200) %>% filter(!is.na(count))

preEGA_periC <- preEGA_pericNULL %>% mutate(fold_norm = preEGA_peric_sum$peric_norm / null_norm)
postEGA_periC <- postEGA_pericNULL %>% mutate(fold_norm = postEGA_peric_sum$peric_norm / null_norm)

periC <- bind_rows(preEGA_periC, postEGA_periC) %>% mutate(int = "Pericentromere")

##imp piRNA data##
preEGA_piRNA_NULL <- read_delim("~/piRNA/preEGA_peaks_null_ints.tsv", 
                                col_names = FALSE, delim = "/", trim_ws = TRUE) %>% distinct() %>% filter(X1 != 0) %>% select(X1, X6) %>% rename(count = X1, loc = X6)
postEGA_piRNA_NULL <- read_delim("~/piRNA/postEGA_peaks_null_ints.tsv", 
                                 col_names = FALSE, delim = "/", trim_ws = TRUE) %>% distinct() %>% filter(X1 != 0) %>% select(X1, X6) %>% rename(count = X1, loc = X6)

preEGA_piRNA <- read_delim("~/piRNA/preEGA_peaks.piRNA.bed")
postEGA_piRNA <- read_delim("~/piRNA/postEGApeaks_comp.piRNA.bed")

##clean and combine##
preEGA_piRNA_sum <- preEGA_piRNA %>% summarise(piRNA_count = n(), piRNA_norm = piRNA_count / 38619) %>% mutate(ega = "pre")
postEGA_piRNA_sum <- postEGA_piRNA %>% summarise(piRNA_count = n(), piRNA_norm = piRNA_count / 36200) %>% mutate(ega = "post")

preEGA_piRNA_NULL <- preEGA_piRNA_NULL %>% mutate(count = as.numeric(count)) %>% mutate(ega = "preEGA") %>% mutate(null_norm = count / 38619) %>% filter(!is.na(count))
postEGA_piRNA_NULL <- postEGA_piRNA_NULL %>% mutate(count = as.numeric(count)) %>% mutate(ega = "postEGA") %>% mutate(null_norm = count / 36200) %>% filter(!is.na(count))


preEGA_piRNA <- preEGA_piRNA_NULL %>% mutate(fold_norm = preEGA_piRNA_sum$piRNA_norm / null_norm)
postEGA_piRNA <- postEGA_piRNA_NULL %>% mutate(fold_norm = postEGA_piRNA_sum$piRNA_norm / null_norm)

piRNA <- bind_rows(preEGA_piRNA, postEGA_piRNA) %>% mutate(int = "piRNA")

total <- bind_rows(periC, piRNA) %>% mutate(ega = factor(ega, levels = c("preEGA", "postEGA")))

##calc means and sds to plot##
total_stats <- total %>% group_by(ega, int) %>% summarise(sd = sd(fold_norm), fold_norm = mean(fold_norm))

int_t.test <- compare_means(fold_norm ~ ega, group.by = "int", data = total, method = "t.test", p.adjust.method = "bonferroni")
int_t.test<- int_t.test %>% mutate(y.position = c(2.3, 3.5))

pre_post_int_piRNA <- ggplot(filter(total_stats, int == "piRNA"), aes(x=ega, y=fold_norm)) +
  geom_col(position = position_dodge(), aes(fill=ega)) +
  geom_errorbar(aes(ymin=fold_norm-sd, ymax=fold_norm+sd), 
                position = position_dodge(), width = 0.05) + 
  scale_fill_manual(values = fill_colors) +
  theme_minimal() + labs(y = "piRNA Fold Enrichment over Control Regions") +
  theme(legend.position = "none") +
  theme(legend.title = element_blank()) +
  theme(axis.title.x = element_blank()) +
  theme(text = element_text(size = 8)) +
  guides(fill = guide_legend(override.aes = list(size = 0.05))) + 
  stat_pvalue_manual(filter(int_t.test, int == "piRNA"), label = "p = {p.adj}", hide.ns = TRUE, size = 2.5, tip.length = 0, 
                     xmin = "group1", xmax = "group2")
ggsave("~/piRNA_per_post.pdf", plot = pre_post_int_piRNA, width = 1.5, height = 2, units = "in")

pre_post_int_peric <- ggplot(filter(total_stats, int == "Pericentromere"), aes(x=ega, y=fold_norm)) +
  geom_col(position = position_dodge(), aes(fill=ega)) +
  geom_errorbar(aes(ymin=fold_norm-sd, ymax=fold_norm+sd), 
                position = position_dodge(), width = 0.05) + 
  scale_fill_manual(values = fill_colors) +
  theme_minimal() + labs(y = "Pericenteromeric Enrichment") +
  theme(legend.position = "none") +
  theme(legend.title = element_blank()) +
  theme(axis.title.x = element_blank()) +
  theme(text = element_text(size = 8)) +
  guides(fill = guide_legend(override.aes = list(size = 0.05))) + 
  stat_pvalue_manual(filter(int_t.test, int == "Pericentromere"), label = "p = {p.adj}", hide.ns = TRUE, size = 2.5, tip.length = 0, 
                     xmin = "group1", xmax = "group2")
ggsave("~/pre_post_int_peric.pdf", plot = pre_post_int_peric, width = 1.5, height = 2, units = "in")


