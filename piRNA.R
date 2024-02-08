library("tidyverse")
library("RColorBrewer")
library("dplyr")
library("ggpubr")
fill_colors <- c("steelblue","gray83", "skyblue3","lightcyan3", "gray60", "lightsteelblue2", "gray40", "skyblue4", "lightskyblue2", "gray83")


####import data####
##doing it w individual transcripts##
names <- c("Chr", "Start", "End", "Name", "map", "read_count")
piRNA_indTEcounts <- read_delim("piRNA_int_wTEs.bed", col_names = names)


ref_ann <- read_delim("TEref_table.txt")
TE_counts_0.1filt <- read_delim("TEann_35_0.1filt_counts.txt", col_names = c("copyNumber_0.1", "Name"))

ref_ann <- left_join(ref_ann, TE_counts_0.1filt)
###calc TPM####
piRNA_ind_counts <- piRNA_indTEcounts %>% mutate(length_kb = (End - Start) / 1000, RPK = read_count / length_kb) %>% filter(!is.na(read_count))
ind_scale <- (summarise(piRNA_ind_counts, scale = sum(RPK)) / 1000000)
piRNA_ind_TPM <- piRNA_ind_counts %>% mutate(TPM = RPK / ind_scale$scale)

new_names <- c("Chr", "Start", "End", "Name", "Repeat")
skinny_ref <- read_delim("TEann_skinny_w_coords.bed", col_names = new_names)

indTE_TPM <- left_join(skinny_ref, (select(piRNA_ind_TPM, c(Chr, Start, End, TPM))))
indTE_TPM <- indTE_TPM %>% mutate(exp = as.integer(TPM >= 0.5)) %>% filter(!is.na(exp))

expTEs <- indTE_TPM %>% filter(TPM >= 0.5)
expTEs_bed <- expTEs %>% select(Chr, Start, End)
write_delim(expTEs_bed, delim = "\t", "expTEs.bed")

name_counts <- indTE_TPM %>% group_by(Name) %>% summarise(copies = n(), exp_copies = sum(exp)) %>% mutate(perc_exp = exp_copies / copies)
mean_TPMs <- indTE_TPM %>% filter(TPM >= 0.5) %>% group_by(Name) %>% summarise(meanTPM = mean(TPM))
name_counts <- left_join(name_counts, mean_TPMs)

df <- left_join(ref_ann, name_counts) %>% filter(Class %in% c("DNA", "SINE", "LINE", "LTR")) %>% filter(isLTR == "FALSE") %>% filter(copyNumber_0.1 >= 10)

###perc_exp by class####
class_piRNA_t.test <- compare_means(perc_exp ~ Class, data = df, method = "t.test", p.adjust.method = "bonferroni", ref.group = ".all.")
class_piRNA_t.test<- class_piRNA_t.test %>%  mutate(y.position = 0.7, 0.7, 0.7, 0.7)

piRNA_byClass <- ggplot(df, aes(x = Class, y = perc_exp)) +
  geom_boxplot(outlier.size = 1, outlier.alpha = 0.5, aes(fill = Class)) +
  labs(y = "Elements with piRNAs (>0.5 TPM)")+
  theme_minimal() + scale_fill_manual(values = fill_colors) +
  theme(legend.position="none") + scale_y_continuous(labels = scales::percent) +
  theme(axis.title.x = element_blank()) +
  theme(text = element_text(size = 8)) +
  stat_pvalue_manual(class_piRNA_t.test, label = "p = {p.adj}", hide.ns = TRUE, remove.bracket = TRUE, size = 2.5)
ggsave("K9_Kinetics/Fig4/piRNA_by_class.pdf", plot = piRNA_byClass, width = 3, height = 3, units = "in")

ggplot(filter(df, Class == "LTR"), aes(x=perc_exp, y=meanTPM, color = Family)) + geom_point() + facet_wrap(~Family)

table <- df %>% select(Name, Family, Class, perc_exp, meanTPM)
write_csv(table, "piRNA_TPM.csv")

###LTRs only####
LTRs_only_table <- table %>% filter(Class == "LTR") %>% filter(Family %in% c("DIRS", "Ngaro", "ERV1", "Gypsy", "Pao", "Copia"))

LTR_fam_t.test <- compare_means(perc_exp ~ Family, data = LTRs_only_table, method = "t.test", p.adjust.method = "bonferroni", ref.group = ".all.")
LTR_fam_t.test<- LTR_fam_t.test %>%  mutate(y.position = 0.7)

ggplot(LTRs_only_table, aes(x = Family, y = perc_exp)) +
  geom_boxplot(outlier.size = 1, outlier.alpha = 0.5, aes(fill = Family)) +
  labs(y = "piRNAs mapping back to TEs (TPM)")+
  theme_minimal() + scale_fill_manual(values = fill_colors) +
  theme(legend.position="none") + scale_y_continuous(labels = scales::percent) +
  theme(axis.title.x = element_blank()) +
  theme(text = element_text(size = 8)) +
  stat_pvalue_manual(LTR_fam_t.test, label = "p = {p.adj}", hide.ns = TRUE, remove.bracket = TRUE, size = 2.5)
