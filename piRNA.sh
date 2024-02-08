#!/bin/bash
###downloaded trimmed fastq files from Kaaij et al 2011
ml STAR
mkdir $OUTDIR/bams

for infile in $BASEDIR/trimmed2/*.gz
do
  base=$(basename ${infile} .fastq.gz)
  STAR --runThreadN 24 --genomeDir $BASEDIR/genome --outFileNamePrefix $BASEDIR/bams/"$base" \
  --readFilesCommand zcat --readFilesIn $infile --outSAMtype BAM SortedByCoordinate \
  --outSAMmultNmax 1 --alignEndsType EndToEnd --alignIntronMax 1
done

##make bed files just for quick visualization
ml BEDTools

for infile in $BASEDIR/bams/*Aligned.sortedByCoord.out.bam
do
  base=$(basename ${infile} Aligned.sortedByCoord.out.bam)
  bedtools bamtobed -i $infile > $BASEDIR/$base.bed
done

###now I want to do a Q1 filter and also make these into bigwigs I think
ml SAMtools

for file in $BASEDIR/bams/*sortedByCoord.out.bam
do
  base=$(basename ${file} Aligned.sortedByCoord.out.bam)
  samtools view -bq1 $file | samtools sort - > $BASEDIR/bams/"$base"_q1.bam
done

samtools merge $BASEDIR/bams/merged_piRNA.bam $BASEDIR/bams/*q1.bam
samtools index $BASEDIR/bams/merged_piRNA.bam

###going to intersect the sRNA at TEs
ml BEDTools
bedtools intersect -a /work/mglab/kld/TEanns/TEann_35_0.1filt.bed -b $BASEDIR/bams/merged_piRNA.bam -c -F 0.3 > $BASEDIR/piRNA_int_wTEs.bed

###in R identified the TEs that had a TPM of over 0.5 and called those "expressed" piRNAs###
###now going to take the coordinates of those TEs and intersect with K9 peaks###
bedtools intersect -a $TCDIR/peaks/4.5hpf_K9_final.bed -b $BASEDIR/expTEs.bed -wa > $BASEDIR/4.5hpfK9peaks_w_piRNA.bed
bedtools intersect -a $TCDIR/peaks/4.5hpf_K9_final.bed -b $BASEDIR/expTEs.bed -v > $BASEDIR/4.5hpfK9peaks_NO_piRNA.bed

bedtools intersect -a $BASEDIR/expTEs.bed -b $TCDIR/peaks/4.5hpf_K9_final.bed -v > $BASEDIR/piRNA_noK9.bed
bedtools intersect -a $BASEDIR/expTEs.bed -b $TCDIR/peaks/4.5hpf_K9_final.bed -wa > $BASEDIR/piRNA_yesK9.bed

bedtools intersect -a $BASEDIR/expTEs_2.bed -b $TCDIR/peaks/3hpf_K9_final.bed -u > $BASEDIR/piRNA_3hK9.bed
bedtools intersect -a /work/mglab/kld/TEanns/TEann_35_0.1filt.bed -b $BASEDIR/piRNA_3hK9.bed -u > $BASEDIR/piRNA_3hK9_wTEs.bed
awk '{print $4}' $BASEDIR/piRNA_3hK9_wTEs.bed | sort - | uniq -c | awk '{print $1 "\t" $2}'> $BASEDIR/piRNA_3hK9_TEs_counts.bed

###pictures at K9 peaks
mkdir $BASEDIR/matrices

ml deepTools
bamCoverage -b $BASEDIR/bams/merged_piRNA.bam --filterRNAstrand forward --normalizeUsing BPM -o $BASEDIR/piRNA_fwd.bw
bamCoverage -b $BASEDIR/bams/merged_piRNA.bam --scaleFactor "-1" --filterRNAstrand reverse --normalizeUsing BPM -o $BASEDIR/piRNA_rev2.bw

computeMatrix scale-regions -S $BASEDIR/piRNA_rev_1.bw $TCDIR/bws/4.5hpf_K9_AVG.bw -R $BASEDIR/piRNA_yesK9.bed $BASEDIR/piRNA_noK9.bed --missingDataAsZero -p max -a 5000 -b 5000 -bs 1 -o $BASEDIR/matrices/piRNAcov_AND_K9cov2.gz
plotHeatmap -m $BASEDIR/matrices/piRNAcov_AND_K9cov2.gz --zMax 3 150 --colorMap Reds Blues --whatToShow "heatmap and colorbar" --heatmapHeight 20 --heatmapWidth 3 --regionsLabel "piRNA targets with H3K9me3" "piRNA targets without H3K9me3" --legendLocation none --startLabel " " --endLabel " " -o $BASEDIR/piRNAcov_AND_K9cov2.pdf

computeMatrix scale-regions -S $BASEDIR/piRNA_fwd_1.bw $BASEDIR/piRNA_rev_1.bw -R $TCDIR/peaks/3hpf_K9_final.bed $TCDIR/peaks/postEGA_peaks.bed --missingDataAsZero -p max -p max -a 500 -b 500 -bs 1 -o $BASEDIR/matrices/K9_3h_4.5h.gz
plotProfile -m $BASEDIR/matrices/K9_3h_4.5h.gz --colors blue red --yMax 1.5 --samplesLabel "fwd strand piRNA" "rev strand piRNA" --regionsLabel "3hpf H3K9me3 peaks" "4.5hpf H3K9me3 peaks" --numPlotsPerRow 1 --startLabel " " --endLabel " " --perGroup -o $BASEDIR/K9_3h_4.5h_profile.pdf

computeMatrix scale-regions -S $BASEDIR/piRNA_rev_1.bw $TCDIR/bws/4.5hpf_K9_AVG.bw -R $BASEDIR/piRNA_yesK9.bed $BASEDIR/piRNA_noK9.bed $BASEDIR/4.5hpfK9peaks_NO_piRNA.bed --missingDataAsZero -p max -a 5000 -b 5000 -bs 1 -o $BASEDIR/matrices/piRNA_K9_triple.gz
plotHeatmap -m $BASEDIR/matrices/piRNA_K9_triple.gz --zMax 3 150 --colorMap Reds Blues --whatToShow "heatmap and colorbar" --heatmapHeight 20 --heatmapWidth 3 --legendLocation none --startLabel " " --endLabel " " -o $BASEDIR/piRNA_K9_triple.pdf
