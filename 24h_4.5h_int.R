library("tidyverse")
library("dplyr")
library("tidyr")
library(ggpubr)
fill_colors <- c("steelblue","gray83", "skyblue3","lightcyan3", "gray60", "lightsteelblue2", "gray40", "skyblue4", "lightskyblue2", "gray83")

####import annotation data######
col_names <- c("count", "Name")

K9_4.5hpf <- read_delim("K9abcam_4.5hpf_counts.bed", col_names = col_names) %>%
  rename(count_4.5h = count)

K9_24hpf <- read_delim("K9abcam_24hpf_counts.bed", col_names = col_names) %>%
  rename(count_24h = count)

comb_TEann <- left_join(K9_4.5hpf, K9_24hpf)

##import ref####
ref_ann <- read_delim("TEref_table.txt")
TE_counts_0.1filt <- read_delim("TEann_35_0.1filt_counts.txt", col_names = c("copyNumber_0.1", "Name"))

ref_ann <- left_join(ref_ann, TE_counts_0.1filt)

####difference in enriched TEs####
TEann_diff <- comb_TEann %>% mutate(morein24 = count_24h - count_4.5h, morein4 = count_4.5h - count_24h) %>%
  mutate(skew = ifelse(morein24 > 0, "Gained", ifelse(morein4 > 0, "Lost", "equal"))) %>% mutate(more = sqrt(morein4^2))

TEdiff_ref <- left_join(ref_ann, select(TEann_diff, Name, count_4.5h, count_24h, skew, more)) %>% filter(Class %in% c("DNA", "LTR", "SINE", "LINE")) %>%
  filter(isLTR == "FALSE") %>% filter(copyNumber_0.1 >= 10) %>% mutate(prop_skew = more / copyNumber_0.1) %>% filter(!is.na(skew))

K9_gain_loss_24h_ttest <- compare_means(prop_skew ~ skew, group.by = "Class", data = filter(TEdiff_ref, skew != "equal"), method = "t.test", p.adjust.method = "bonferroni")
K9_gain_loss_24h_ttest <- K9_gain_loss_24h_ttest %>% mutate(y.position = c(.23))


K9_gain_loss_24h <- ggplot(filter(TEdiff_ref, skew != "equal"), aes(x = Class, y = prop_skew)) +
  geom_boxplot(outlier.size = 1, outlier.alpha = 0.5, aes(fill = skew))+
  labs(y = "Differentially H3K9me3-Enriched Elements")+
  theme_minimal() + scale_fill_manual(values = fill_colors) +
  theme(legend.position="right") + theme(legend.title = element_blank()) +
  theme(axis.title.x = element_blank()) +
  theme(text = element_text(size = 8)) +
  scale_y_continuous(labels = scales::percent) +
  stat_pvalue_manual(K9_gain_loss_24h_ttest, label = "p = {p.adj}", hide.ns = TRUE, size = 2.5, tip.length = 0.01,
                     xmin = "Class", xmax = "Class")
ggsave("K9_Kinetics/Fig6/K9_gain_loss_24h.pdf", plot = K9_gain_loss_24h, width = 3, height = 2.8, units = "in")

TEdiff_ref <- TEdiff_ref %>% filter(skew != "equal")
###DNA and LTR family-level###
LTR_only <- TEdiff_ref %>% filter(Class %in% c("LTR")) %>% filter(Family %in% c("ERV1", "Gypsy", "Pao", "Ngaro")) %>% group_by(Family) %>% filter(n() > 5) %>%  ungroup()

LTR_gain_loss_24h_ttest <- compare_means(prop_skew ~ skew, group.by = "Family", data = filter(LTR_only, skew != "equal"), method = "t.test", p.adjust.method = "bonferroni")
LTR_gain_loss_24h_ttest <- LTR_gain_loss_24h_ttest %>% mutate(y.position = c(.21))

LTR_loss <- ggplot(LTR_only, aes(x = Family, y = prop_skew)) +
  geom_boxplot(outlier.size = 1, outlier.alpha = 0.5, aes(fill = skew))+
  labs(y = "Differentially H3K9me3-Enriched Elements")+
  theme_minimal() + scale_fill_manual(values = fill_colors) +
  theme(legend.position="right") + theme(legend.title = element_blank()) +
  theme(axis.title.x = element_blank()) +
  theme(text = element_text(size = 8)) +
  scale_y_continuous(labels = scales::percent) +
  stat_pvalue_manual(LTR_gain_loss_24h_ttest, label = "p = {p.adj}", hide.ns = TRUE, size = 2.5, tip.length = 0.01,
                     xmin = "Family", xmax = "Family")
ggsave("K9_Kinetics/Fig6/LTR_loss.pdf", plot = LTR_loss, width = 3, height = 2.8, units = "in")


DNA_only <- TEdiff_ref %>% filter(Class %in% c("DNA")) %>% filter(Family %in% c("TcMar-Tc1", "CMC-EnSpm", "DNA", "Kolobok-T2", "PIF-Harbinger",
                                                                                "hAT-Ac", "hAT-Tip100", "hAT-Charlie", "TcMar-ISRm11", "IS3EU")) %>%
  group_by(Family) %>% filter(n() > 5) %>% ungroup()
DNA_gain_loss_24h_ttest <- compare_means(prop_skew ~ skew, group.by = "Family", data = DNA_only, method = "t.test", p.adjust.method = "bonferroni")
DNA_gain_loss_24h_ttest <- DNA_gain_loss_24h_ttest %>% mutate(y.position = 0.12)

DNA_skinny <- DNA_only %>% group_by(Family) %>% filter(n() >50) %>% ungroup %>% filter(prop_skew < 0.12)

DNA_fam_gain <- ggplot(DNA_skinny, aes(x = Family, y = prop_skew)) +
  geom_boxplot(outlier.size = 1, outlier.alpha = 0.5, aes(fill = skew))+
  labs(y = "Differentially H3K9me3-Enriched Elements")+
  theme_minimal() + scale_fill_manual(values = fill_colors) +
  theme(legend.position="none") + theme(legend.title = element_blank()) +
  theme(axis.title.x = element_blank()) +
  theme(text = element_text(size = 8)) +
  scale_y_continuous(labels = scales::percent) +
  stat_pvalue_manual(DNA_gain_loss_24h_ttest, label = "p = {p.adj}", hide.ns = TRUE, size = 2.5, tip.length = 0.01,
                     xmin = "Family", xmax = "Family")
ggsave("K9_Kinetics/Fig6/DNA_fam_gain_noleg.pdf", plot = DNA_fam_gain, width = 3.85, height = 2.25, units = "in")

####venn diagram of peak overlap####
library(eulerr)

venn_data <-  euler(c("4.5hpf" = 33771, "24hpf" =  36575, "4.5hpf&24hpf" = 41104))

plot(venn_data,
      quantities = TRUE,
      labels = FALSE)
