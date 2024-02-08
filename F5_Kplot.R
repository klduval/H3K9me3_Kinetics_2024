###here we are going to make some karyoplots for the TE paper 
library("karyoploteR")
library("tidyverse")
library("dplyr")
library("tidyr")

genome <- read.table("/Users/mangofrog7/Documents/Graduate_School/Projects/bioinformatics/general/chrNL.txt", col.names = c("1", "2", "3"))
Z11_genome <- toGRanges(genome)
periC_coords <- read.table("/Users/mangofrog7/Documents/Graduate_School/Projects/bioinformatics/periC_ann/centromeres_sloppy2.bed", col.names = c("1", "2", "3"))
periC_coords<- periC_coords %>% mutate(gieStain = "acen")
cent_comp <- read.table("/Users/mangofrog7/Documents/Graduate_School/Projects/bioinformatics/periC_ann/cent_comp.bed", col.names = c("1", "2", "3")) %>%
  mutate(gieStain = "gvar")
periC_data <- bind_rows(periC_coords, cent_comp)
periC <- toGRanges(periC_data)


preEGA_K9 <- read.table("/Users/mangofrog7/Documents/Graduate_School/Projects/bioinformatics/TC_final/peaks/preEGA_peaks_total.bed") %>% select(., "V1", "V2", "V3")
preEGA_K9 <- toGRanges(preEGA_K9)
postEGA_K9 <- read.table("/Users/mangofrog7/Documents/Graduate_School/Projects/bioinformatics/TC_final/peaks/postEGApeaks_comp.bed") %>% select(., "V1", "V2", "V3")
postEGA_K9 <- toGRanges(postEGA_K9)

kp <- plotKaryotype(genome = Z11_genome, cytobands = periC, plot.type = 2)
kpPlotDensity(kp, data=preEGA_K9, window.size = 350000, col="steelblue")
kpPlotDensity(kp, data=postEGA_K9, window.size = 350000, col="gray83", data.panel = 2)


kp_abb <- plotKaryotype(genome = Z11_genome, cytobands = periC, plot.type = 2, chr = c("2", "8", "18", "20", "24", "25"))
kpPlotDensity(kp_abb, data=preEGA_K9, window.size = 350000, col="steelblue")
kpPlotDensity(kp_abb, data=postEGA_K9, window.size = 350000, col="gray83", data.panel = 2)
