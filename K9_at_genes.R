library("tidyverse")
library("dplyr")
library("ggpubr")
fill_colors <- c("steelblue","gray83", "skyblue3","lightcyan3", "gray60", "lightsteelblue2", "gray40", "skyblue4", "lightskyblue2", "gray83")

###first venn of genic K9 peaks####
library(eulerr)

K9_at_genes <-  euler(c("TSS" = 2234, "TSS&protein_coding" =  1041, "TSS&protein_coding&notTE" = 939))

plot(K9_at_genes,
     quantities = TRUE,
     labels = FALSE)


###GO analysis fig###
GO <- read.csv("4.5hpf_K9_timecourseGOanalysis.csv")

GO <- GO %>% mutate(
  source = fct_relevel(source, "TF", "GO:BP", "REAC", "GO:CC", "HP"),
  term_name = fct_reorder(term_name, negative_log10_of_adjusted_p_value)
)

GOplot <- ggplot(GO, aes(x=term_name, y=negative_log10_of_adjusted_p_value, fill = source)) +
  geom_col(width = 0.85) +
  facet_grid(cols=vars(source), scales = "free", space = "free", switch = "y") +
  theme_minimal() + labs(y = "-log10(p.adj)") +
  theme(axis.title.x = element_blank()) +
  scale_fill_manual(values = fill_colors) +
  theme(legend.position="none") +
  theme(axis.text.x = element_text(angle = 90, hjust =1)) +
  theme(axis.text.y = element_text(angle =90)) +
  theme(text = element_text(size = 8))
ggsave("K9_Kinetics/FigS3/GOplot.pdf", plot = GOplot, width = 2.25, height = 3.75, units = "in")
