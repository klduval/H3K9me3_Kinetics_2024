###here we are going to make some karyoplots for the paper 
library("karyoploteR")
library("tidyverse")
library("dplyr")
library("tidyr")

genome <- read.table("/Users/mangofrog7/Documents/Graduate_School/Projects/bioinformatics/general/chrNL.txt", col.names = c("1", "2", "3"))
Z11_genome <- toGRanges(genome)

pp <- getDefaultPlotParams(plot.type=2)
pp$ideogramheight <- 25

peaks_4.5hpf <- read.table("/Users/mangofrog7/Documents/Graduate_School/Projects/bioinformatics/TC_final/peaks/4.5hpf_K9_final.bed") %>% select(., "V1", "V2", "V3")
peaks_4.5hpf <- toGRanges(peaks_4.5hpf)
BRSAT1_coords <- read.table("/Users/mangofrog7/Documents/Graduate_School/Projects/bioinformatics/TC_final/pericentromere/BRSAT1_final.bed") %>% select(., "V1", "V2", "V3")
BRSAT1 <- toGRanges(BRSAT1_coords)
K9_4.5hpf_250kb_bw <- "/Users/mangofrog7/Documents/Graduate_School/Projects/bioinformatics/TC_final/bws/4.5hpf_K9_250kbBS_AVG.bw"

kp <- plotKaryotype(genome = Z11_genome, plot.type = 2, plot.params = pp, chromosomes = c("9", "12", "18", "25"))
kpAddBaseNumbers(kp)
kpPlotBigWig(kp, data=K9_4.5hpf_250kb_bw, ymax="visible.region", col="midnightblue")
kpPlotDensity(kp, data=BRSAT1, window.size = 50000, data.panel = 2, col="red4", r0 = 0, r1 = 0.75)

kp_full <- plotKaryotype(genome = Z11_genome, plot.type = 2, plot.params = pp)
kpAddBaseNumbers(kp_full)
kpPlotBigWig(kp_full, data=K9_4.5hpf_250kb_bw, ymax="visible.region", col="midnightblue")
kpPlotDensity(kp_full, data=BRSAT1, window.size = 50000, data.panel = 2, col="red4", r0 = 0, r1 = 0.75)

kp2 <- plotKaryotype(genome = Z11_genome, plot.params = pp, chromosomes = c("9", "12", "18", "25", "4"))
kpPlotDensity(kp2, data=peaks_4.5hpf, window.size = 300000, col="midnightblue")

kp2_full <- plotKaryotype(genome = Z11_genome, plot.params = pp)
kpPlotDensity(kp2_full, data=peaks_4.5hpf, window.size = 300000, col="midnightblue")
