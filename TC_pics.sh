#!/bin/bash
###lets make these bedgraphs into bigwigs for data visualization
module load ucsc/359
mkdir $BASEDIR/bws

for infile in $BASEDIR/bdgrphs/*norm.bga
do
  base=$(basename ${infile} .norm.bga)
  bedSort $infile $infile
  bedGraphToBigWig $infile $BASEDIR/genome/chrNameLength.txt $BASEDIR/bws/$base.bw
done

###lets do some broad comparisons to see what our data looks like before moving on
module load deepTools

multiBigwigSummary bins -b $BASEDIR/bws/*[1-3].bw $BASEDIR/bws/*IgG*.bw -o $BASEDIR/bwreps_summ.npz -p 24
plotCorrelation -in $BASEDIR/bwreps_summ.npz -c spearman -p heatmap -o $BASEDIR/timecourse_bwreps_summ_heatmap.pdf
plotPCA -in $BASEDIR/bwreps_summ.npz -o $BASEDIR/timecourse_bwreps_summ_PCA.pdf

###I want to make merged/average bigwig files
bigwigCompare -b1 $BASEDIR/bws/2hpf_K9_1.bw -b2 $BASEDIR/bws/2hpf_K9_2.bw --operation mean -bs 10 -p 20 -o $BASEDIR/bws/2hpf_K9_AVG.bw
bigwigCompare -b1 $BASEDIR/bws/2.5hpf_K9_1.bw -b2 $BASEDIR/bws/2.5hpf_K9_2.bw --operation add --scaleFactors 0.333:0.333 -bs 10 -p 20 -o $BASEDIR/bws/2.5hpf_K9_rep1rep2.bw
bigwigCompare -b1 $BASEDIR/bws/2.5hpf_K9_rep1rep2.bw -b2 $BASEDIR/bws/2.5hpf_K9_3.bw --operation add --scaleFactors 1:0.333 -bs 10 -p 20 -o $BASEDIR/bws/2.5hpf_K9_AVG.bw
bigwigCompare -b1 $BASEDIR/bws/3hpf_K9_1.bw -b2 $BASEDIR/bws/3hpf_K9_2.bw --operation add --scaleFactors 0.333:0.333 -bs 10 -p 20 -o $BASEDIR/bws/3hpf_K9_rep1rep2.bw
bigwigCompare -b1 $BASEDIR/bws/3hpf_K9_rep1rep2.bw -b2 $BASEDIR/bws/3hpf_K9_3.bw --operation add --scaleFactors 1:0.333 -bs 10 -p 20 -o $BASEDIR/bws/3hpf_K9_AVG.bw
bigwigCompare -b1 $BASEDIR/bws/3.5hpf_K9_1.bw -b2 $BASEDIR/bws/3.5hpf_K9_2.bw --operation add --scaleFactors 0.333:0.333 -bs 10 -p 20 -o $BASEDIR/bws/3.5hpf_K9_rep1rep2.bw
bigwigCompare -b1 $BASEDIR/bws/3.5hpf_K9_rep1rep2.bw -b2 $BASEDIR/bws/3.5hpf_K9_3.bw --operation add --scaleFactors 1:0.333 -bs 10 -p 20 -o $BASEDIR/bws/3.5hpf_K9_AVG.bw
bigwigCompare -b1 $BASEDIR/bws/4hpf_K9_1.bw -b2 $BASEDIR/bws/4hpf_K9_2.bw --operation add --scaleFactors 0.333:0.333 -bs 10 -p 20 -o $BASEDIR/bws/4hpf_K9_rep1rep2.bw
bigwigCompare -b1 $BASEDIR/bws/4hpf_K9_rep1rep2.bw -b2 $BASEDIR/bws/4hpf_K9_3.bw --operation add --scaleFactors 1:0.333 -bs 10 -p 20 -o $BASEDIR/bws/4hpf_K9_AVG.bw
bigwigCompare -b1 $BASEDIR/bws/4.5hpf_K9_1.bw -b2 $BASEDIR/bws/4.5hpf_K9_2.bw --operation add --scaleFactors 0.333:0.333 -bs 10 -p 20 -o $BASEDIR/bws/4.5hpf_K9_rep1rep2.bw
bigwigCompare -b1 $BASEDIR/bws/4.5hpf_K9_rep1rep2.bw -b2 $BASEDIR/bws/4.5hpf_K9_3.bw --operation add --scaleFactors 1:0.333 -bs 10 -p 20 -o $BASEDIR/bws/4.5hpf_K9_AVG.bw

rm $BASEDIR/bws/*rep1rep2*

####making bigwig with larger bin size for Kplots
bigwigCompare -b1 $BASEDIR/bws/4.5hpf_K9_1.bw -b2 $BASEDIR/bws/4.5hpf_K9_2.bw --operation add --scaleFactors 0.333:0.333 -bs 250000 -p 20 -o $BASEDIR/bws/4.5hpf_K9_rep1rep2.bw
bigwigCompare -b1 $BASEDIR/bws/4.5hpf_K9_rep1rep2.bw -b2 $BASEDIR/bws/4.5hpf_K9_3.bw --operation add --scaleFactors 1:0.333 -bs 250000 -p 20 -o $BASEDIR/bws/4.5hpf_K9_250kbBS_AVG.bw


###making pictures now
computeMatrix reference-point -S $BASEDIR/bws/2hpf_K9_AVG.bw $BASEDIR/bws/2.5hpf_K9_AVG.bw $BASEDIR/bws/3hpf_K9_AVG.bw $BASEDIR/bws/3.5hpf_K9_AVG.bw $BASEDIR/bws/4hpf_K9_AVG.bw $BASEDIR/bws/4.5hpf_K9_AVG.bw -R $BASEDIR/peaks/4.5hpf_K9_final.bed --referencePoint center -p max -a 10000 -b 10000 -bs 10 --missingDataAsZero -o $BASEDIR/matrices/TC_4.5peaks.gz
computeMatrix reference-point -S $BASEDIR/bws/3.5hpf_K9_AVG.bw $BASEDIR/bws/4hpf_K9_AVG.bw $BASEDIR/bws/4.5hpf_K9_AVG.bw -R $BASEDIR/peaks/4.5hpf_K9_final.bed --referencePoint center -p max -a 10000 -b 10000 -bs 10 --missingDataAsZero -o $BASEDIR/matrices/lateTC_4.5peaks.gz
computeMatrix reference-point -S $BASEDIR/bws/2hpf_K9_AVG.bw $BASEDIR/bws/2.5hpf_K9_AVG.bw $BASEDIR/bws/3hpf_K9_AVG.bw -R $BASEDIR/peaks/4.5hpf_K9_final.bed --referencePoint center -p max -a 10000 -b 10000 -bs 10 --missingDataAsZero -o $BASEDIR/matrices/earlyTC_4.5peaks.gz

for infile in $BASEDIR/matrices/*.gz
do
  base=$(basename ${infile} .gz)
  plotHeatmap -m $infile --colorMap Blues --legendLocation none --regionsLabel "Peaks" -o $BASEDIR/figs/"$base"_heatmap.pdf
done

###going to make matrices so I can compute peak "height" for quantification
computeMatrix scale-regions -S $BASEDIR/bws/2hpf_K9_AVG.bw -R $BASEDIR/peaks/2hpf_K9_final.bed --missingDataAsZero -p max --outFileNameMatrix $BASEDIR/2hpf_peakmatrix.tab -o $BASEDIR/matrices/2hpf.gz
computeMatrix scale-regions -S $BASEDIR/bws/2.5hpf_K9_AVG.bw -R $BASEDIR/peaks/2.5hpf_K9_final.bed --missingDataAsZero -p max --outFileNameMatrix $BASEDIR/2.5hpf_peakmatrix.tab -o $BASEDIR/matrices/2.5hpf.gz
computeMatrix scale-regions -S $BASEDIR/bws/3hpf_K9_AVG.bw -R $BASEDIR/peaks/3hpf_K9_final.bed --missingDataAsZero -p max --outFileNameMatrix $BASEDIR/3hpf_peakmatrix.tab -o $BASEDIR/matrices/3hpf.gz
computeMatrix scale-regions -S $BASEDIR/bws/3.5hpf_K9_AVG.bw -R $BASEDIR/peaks/3.5hpf_K9_final.bed --missingDataAsZero -p max --outFileNameMatrix $BASEDIR/3.5hpf_peakmatrix.tab -o $BASEDIR/matrices/3.5hpf.gz
computeMatrix scale-regions -S $BASEDIR/bws/4hpf_K9_AVG.bw -R $BASEDIR/peaks/4hpf_K9_final.bed --missingDataAsZero -p max --outFileNameMatrix $BASEDIR/4hpf_peakmatrix.tab -o $BASEDIR/matrices/4hpf.gz
computeMatrix scale-regions -S $BASEDIR/bws/4.5hpf_K9_AVG.bw -R $BASEDIR/peaks/4.5hpf_K9_final.bed --missingDataAsZero -p max --outFileNameMatrix $BASEDIR/4.5hpf_peakmatrix.tab  -o $BASEDIR/matrices/4.5hpf.gz

###heatmap comparing enrichment across all my 2.5hpf samples ever
#first need to pull the ranks regions from the original matrix
computeMatrix reference-point -S $BASEDIR/bws/2.5hpf_K9_AVG.bw -R $BASEDIR/peaks/2.5hpf_K9_final.bed --referencePoint center -p max -a 10000 -b 10000 -bs 10 --missingDataAsZero -o $BASEDIR/matrices/2.5hpf_peaks.gz
plotHeatmap -m $BASEDIR/matrices/2.5hpf_peaks.gz --colorMap Blues --legendLocation none --regionsLabel "2.5 hpf Peaks" --outFileSortedRegions $BASEDIR/2.5hpf_heatmap_regions.bed -o $BASEDIR/figs/2.5hpf_heatmap.pdf

antiDIR="/scratch/kld57880/antiK9_5.2022"
computeMatrix reference-point -S $antiDIR/bws/K9abcam_2.5hpf_AVG.bw -R $BASEDIR/2.5hpf_heatmap_regions.bed --referencePoint center -p max -a 10000 -b 10000 -bs 10 --sortRegions keep --missingDataAsZero -o $BASEDIR/matrices/2.5hpf_abcam2.gz
computeMatrix reference-point -S $antiDIR/bws/K9active_2.5hpf_AVG.bw -R $BASEDIR/2.5hpf_heatmap_regions.bed --referencePoint center -p max -a 10000 -b 10000 -bs 10 --sortRegions keep --missingDataAsZero -o $BASEDIR/matrices/2.5hpf_active.gz
computeMatrix reference-point -S $antiDIR/bws/K9diag_2.5hpf_AVG.bw -R $BASEDIR/2.5hpf_heatmap_regions.bed --referencePoint center -p max -a 10000 -b 10000 -bs 10 --sortRegions keep --missingDataAsZero -o $BASEDIR/matrices/2.5hpf_diag.gz

plotHeatmap -m $BASEDIR/matrices/2.5hpf_abcam2.gz --colorMap Blues --legendLocation none --regionsLabel "2.5 hpf Peaks" --sortRegions keep -o $BASEDIR/figs/2.5hpf_abcam2.pdf
plotHeatmap -m $BASEDIR/matrices/2.5hpf_active.gz --colorMap Blues --legendLocation none --regionsLabel "2.5 hpf Peaks" --sortRegions keep -o $BASEDIR/figs/2.5hpf_active.pdf
plotHeatmap -m $BASEDIR/matrices/2.5hpf_diag.gz --colorMap Blues --legendLocation none --regionsLabel "2.5 hpf Peaks" --sortRegions keep -o $BASEDIR/figs/2.5hpf_diag.pdf

###make a heatmap comparing seeded and non-seeded peaks
computeMatrix reference-point -S $BASEDIR/bws/3hpf_K9_AVG.bw -R $BASEDIR/peaks/preEGA_peaks_total.bed $BASEDIR/peaks/postEGApeaks_comp.bed --referencePoint center -p max -a 10000 -b 10000 --missingDataAsZero -bs 10 -o $BASEDIR/matrices/3hpf_preVpostEGA_peaks.gz
plotHeatmap -m $BASEDIR/matrices/3hpf_preVpostEGA_peaks.gz --colorMap Blues --yMax 5 --legendLocation none --whatToShow "plot and heatmap" --regionsLabel "preEGA Peaks" "postEGA Peaks" -o $BASEDIR/figs/3hpf_preVpostEGA_peaks.pdf
