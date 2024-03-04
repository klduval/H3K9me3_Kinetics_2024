library("tidyverse")
library("RColorBrewer")
library("dplyr")
library("ggpubr")
fill_colors <- c("steelblue","gray83", "skyblue3","lightcyan3", "gray60", "lightsteelblue2", "gray40", "skyblue4", "lightskyblue2", "gray83")

####import data####
##doing it w individual transcripts##
names <- c("Chr", "Start", "End", "Name", "map", "read_count")
piRNA_indTEcounts <- read_delim("~/piRNA/piRNA_int_wTEs.bed", col_names = names)


ref_ann <- read_delim("~/general/TEref_table.txt")  
TE_counts_0.1filt <- read_delim("~/general/TEann_35_0.1filt_counts.txt", col_names = c("copyNumber_0.1", "Name"))  

ref_ann <- left_join(ref_ann, TE_counts_0.1filt)
###calc TPM####
piRNA_ind_counts <- piRNA_indTEcounts %>% mutate(length_kb = (End - Start) / 1000, RPK = read_count / length_kb) %>% filter(!is.na(read_count))
ind_scale <- (summarise(piRNA_ind_counts, scale = sum(RPK)) / 1000000)
piRNA_ind_TPM <- piRNA_ind_counts %>% mutate(TPM = RPK / ind_scale$scale)

new_names <- c("Chr", "Start", "End", "Name", "Repeat")
skinny_ref <- read_delim("~/piRNA/TEann_skinny_w_coords.bed", col_names = new_names)

indTE_TPM <- left_join(skinny_ref, (select(piRNA_ind_TPM, c(Chr, Start, End, TPM))))
indTE_TPM <- indTE_TPM %>% mutate(exp = as.integer(TPM >= 0.5)) %>% filter(!is.na(exp))

expTEs <- indTE_TPM %>% filter(TPM >= 0.5)
expTEs_bed <- expTEs %>% select(Chr, Start, End)
write_delim(expTEs_bed, delim = "\t", "~/piRNA/expTEs.bed")

name_counts <- indTE_TPM %>% group_by(Name) %>% summarise(copies = n(), exp_copies = sum(exp)) %>% mutate(perc_exp = exp_copies / copies)
mean_TPMs <- indTE_TPM %>% filter(TPM >= 0.5) %>% group_by(Name) %>% summarise(meanTPM = mean(TPM))
name_counts <- left_join(name_counts, mean_TPMs)

df <- left_join(ref_ann, name_counts) %>% filter(Class %in% c("DNA", "SINE", "LINE", "LTR")) %>% filter(isLTR == "FALSE") %>% filter(copyNumber_0.1 >= 10)

###piRNAs that trigger K9 in the periC#####
piRNA_peric <- read_delim("~/piRNA/expTE.peric.bed", col_names = c("Chr", "start", "end"))
piRNA_K9 <- read_delim("~/piRNA/piRNA_3hK9_wTEs.bed", col_names = c("Chr", "start", "end", "Name", "TPM"))

##get only piRNAs that are in the peric and also overlap preEGA K9
piRNA_K9_peric <- left_join(piRNA_peric, piRNA_K9) %>% filter(!is.na(Name))

##summarise 
piRNA_K9_peric_sum <- piRNA_K9_peric %>% group_by(Name) %>% summarise(count = n())

##bring in ref
ref_ann <- read_delim("~/general/TEref_table.txt")  
TE_counts_0.1filt <- read_delim("~/general/TEann_35_0.1filt_counts.txt", col_names = c("copyNumber_0.1", "Name"))  

ref_ann <- left_join(ref_ann, TE_counts_0.1filt)

##combine
ref_w_periC_K9_piRNA <- left_join(ref_ann, piRNA_K9_peric_sum) %>% filter(!is.na(count), Class != "RC")

by_Class <- ref_w_periC_K9_piRNA %>% group_by(Class) %>% summarise(total_count = sum(count))

byClass_prop <- by_Class %>% mutate(prop = total_count / (sum(by_Class$total_count)))

piRNA_Class_pic <- ggplot(byClass_prop, aes(x = Class, y = prop, fill = Class)) + 
  geom_col(position = position_dodge()) +
  theme_minimal() + scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = fill_colors) + 
  theme(text = element_text(size = 8)) +
  labs(y = "Proportion of pre-EGA Pericentromeric piRNAs") +
  theme(legend.position = "none") +
  theme(axis.title.x = element_blank())
ggsave("~/piRNA_Class_pic.pdf", plot = piRNA_Class_pic, width = 2.15, height = 2, units = "in")

###piRNA K9 peak peric location####
piRNA_cats <- tibble(location = c("Pericentromere", "Distal"), 
                     count = c(1474, 6533)) %>% mutate(prop = count / 8007)

piRNA_cats_pie <- ggplot(piRNA_cats, aes(x = "", y = prop, fill = location)) + 
  geom_bar(stat = "identity", width = 1, color = "white") + coord_polar("y", start=0) +
  theme_void() + 
  scale_fill_manual(values = fill_colors) +
  theme(text = element_text(size = 8)) +
  theme(axis.title.x = element_blank(), legend.position = "none")
ggsave("~/piRNA_location_pie.pdf", plot = piRNA_cats_pie, width = 2, height = 2, units = "in")

