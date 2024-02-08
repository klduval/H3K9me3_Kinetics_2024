library("tidyverse")
library("RColorBrewer")
library("dplyr")


TE_35bp_scores <- read.delim("/Users/mangofrog7/Documents/Graduate_School/Projects/bioinformatics/general/TE_35bp_coords_score.bed", 
                              col.names = c("chr", "start", "end", "score"))
TE_ann <- read.delim("/Users/mangofrog7/Documents/Graduate_School/Projects/bioinformatics/general/TEann_final.bed", 
                             col.names = c("chr", "start", "end", "name"))

TEann_w35bp_mapscores <- left_join(TE_ann, TE_35bp_scores, multiple = "first")

write.table(TEann_w35bp_mapscores, "/Users/mangofrog7/Documents/Graduate_School/Projects/bioinformatics/general/TEann_w35bp_mapscores.bed", 
            row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
