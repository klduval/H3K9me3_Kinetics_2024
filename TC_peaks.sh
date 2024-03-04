#!/bin/bash
#SBATCH --job-name=K9_timecourse
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=40gb
#SBATCH --time=12:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kld57880@uga.edu

BASEDIR="/scratch/kld57880/TC_final"
TOOLDIR='/home/kld57880/Git2/toolbox'

###peak calling
module load Homer
mkdir $BASEDIR/peaks

for infile in $BASEDIR/bdgrphs/*.norm.bga
  do base=$(basename ${infile} .norm.bga)
  cat $infile | awk '{print $1 "\t" $2 "\t" $3 "\t" "+" "\t" "+" "\t" "+"}' > $BASEDIR/peaks/$base.bgato.bed
done

for infile in $BASEDIR/peaks/*bgato.bed
  do base=$(basename ${infile} .bgato.bed)
  makeTagDirectory $BASEDIR/peaks/$base.BtB.tagdir $infile -format bed
done

for infile in $BASEDIR/peaks/*K9*.tagdir
  do base=$(basename ${infile} .BtB.tagdir)
  findPeaks $infile -style histone -minDist 1000 -gsize 1.5e9 -F 4 -i $BASEDIR/peaks/*mIgG*.tagdir -o $BASEDIR/peaks/$base.txt
done

for infile in $BASEDIR/peaks/*.txt
do
  base=$(basename ${infile} .txt)
  sed '/^#/d' $infile | awk '{print $2 "\t" $3 "\t" $4 "\t" $1 "\t" $8 "\t" $5 "\t" $6 "\t" $12 "\t" "-1"}' | sed 's/\.000000//g' > $BASEDIR/peaks/$base.peaks.bed
done

module load ChIP-R/1.1.0-foss-2019b-Python-3.7.4

chipr -i $BASEDIR/peaks/2hpf_K9_1.peaks.bed $BASEDIR/peaks/2hpf_K9_2.peaks.bed -m 2 -o $BASEDIR/peaks/2hpf_K9_repPeaks
chipr -i $BASEDIR/peaks/2.5hpf_K9_1.peaks.bed $BASEDIR/peaks/2.5hpf_K9_2.peaks.bed $BASEDIR/peaks/2.5hpf_K9_3.peaks.bed -m 2 -o $BASEDIR/peaks/2.5hpf_K9_repPeaks
chipr -i $BASEDIR/peaks/3hpf_K9_1.peaks.bed $BASEDIR/peaks/3hpf_K9_2.peaks.bed $BASEDIR/peaks/3hpf_K9_3.peaks.bed -m 2 -o $BASEDIR/peaks/3hpf_K9_repPeaks
chipr -i $BASEDIR/peaks/3.5hpf_K9_1.peaks.bed $BASEDIR/peaks/3.5hpf_K9_2.peaks.bed $BASEDIR/peaks/3.5hpf_K9_3.peaks.bed -m 2 -o $BASEDIR/peaks/3.5hpf_K9_repPeaks
chipr -i $BASEDIR/peaks/4hpf_K9_1.peaks.bed $BASEDIR/peaks/4hpf_K9_2.peaks.bed $BASEDIR/peaks/4hpf_K9_3.peaks.bed -m 2 -o $BASEDIR/peaks/4hpf_K9_repPeaks
chipr -i $BASEDIR/peaks/4.5hpf_K9_1.peaks.bed $BASEDIR/peaks/4.5hpf_K9_2.peaks.bed $BASEDIR/peaks/4.5hpf_K9_3.peaks.bed -m 2 -o $BASEDIR/peaks/4.5hpf_K9_repPeaks

###make a blacklist file
findPeaks $BASEDIR/peaks/tagdirs/mIgG.tagdir -style factor -o $BASEDIR/peaks/IgG.txt
sed '/^#/d' $BASEDIR/peaks/IgG.txt | awk '{print $2 "\t" $3 "\t" $4 "\t" $1 "\t" "1" "\t" $5 "\t" $6 "\t" $12 "\t" "-1"}' > $BASEDIR/peaks/blacklist.bed

ml BEDTools
###intersect the peaks with the blacklist file to make sure we aren't looking at sticky regions before this step
for infile in $BASEDIR/peaks/*all.bed
do
  base=$( basename ${infile} _repPeaks_all.bed)
  bedtools intersect -a $infile -b $BASEDIR/peaks/blacklist.bed -v > $BASEDIR/peaks/"$base"_final.bed
done

###peak annotation####
curl -s ftp://ftp.ensembl.org/pub/release-98/gtf/danio_rerio/Danio_rerio.GRCz11.98.gtf.gz | gunzip -c > $BASEDIR/refann.gtf
mkdir $BASEDIR/peaks/ann

for infile in $BASEDIR/peaks/*final.bed
do
  base=$( basename ${infile} final.bed)
  annotatePeaks.pl $infile danRer11 -gtf $BASEDIR/refann.gtf > $BASEDIR/peaks/ann/$base.maskann.txt
done

for infile in $BASEDIR/peaks/ann/*maskann.txt
do
  base=$(basename ${infile} .maskann.txt)
  awk -F'\t' 'sqrt($10*$10) <=1000' $infile > $BASEDIR/peaks/ann/$base.1000bp_ann.txt
done

for infile in $BASEDIR/peaks/ann/*maskann.txt
do
  base=$(basename ${infile} .maskann.txt)
  awk -F'\t' 'sqrt($10*$10) >=1000' $infile | awk '{print $2 "\t" $3 "\t" $4 }' > $BASEDIR/peaks/ann/$base.MOREthan1000bp.bed
done

for infile in $BASEDIR/peaks/*final.bed
do
  base=$( basename ${infile} final.bed)
    bedtools intersect -a /work/mglab/kld/TEanns/*0.1*.bed -b $infile -f 0.50 -u > $BASEDIR/peaks/ann4/$base.TEann.txt
done

for infile in $BASEDIR/peaks/ann4/*_.TEann.txt
do
  base=$(basename ${infile} _.TEann.txt)
  awk '{print $4}' $infile | sort - | uniq -c | awk '{print $1 "\t" $2}'> $BASEDIR/"$base"_TEcounts2.bed
done

####timcecourse intersection####
sh $TOOLDIR/peaks_over_time.sh -o $BASEDIR/peaks -m merged -t 2hpf -t 2.5hpf -t 3hpf -t 3.5hpf -t 4hpf -t 4.5hpf $BASEDIR/peaks/2hpf_K9_final.bed $BASEDIR/peaks/2.5hpf_K9_final.bed $BASEDIR/peaks/3hpf_K9_final.bed $BASEDIR/peaks/3.5hpf_K9_final.bed $BASEDIR/peaks/4hpf_K9_final.bed $BASEDIR/peaks/4.5hpf_K9_final.bed

cat $BASEDIR/peaks/2hpf_K9_final.bed $BASEDIR/peaks/2.5hpf_K9_final.bed $BASEDIR/peaks/3hpf_K9_final.bed | bedtools sort | bedtools merge -i - > $BASEDIR/peaks/preEGA_peaks_total.bed
cat $BASEDIR/peaks/3.5hpf_K9_final.bed $BASEDIR/peaks/4hpf_K9_final.bed $BASEDIR/peaks/4.5hpf_K9_final.bed | bedtools sort | bedtools merge -i - > $BASEDIR/peaks/postEGA_peaks_total.bed
bedtools intersect -a $BASEDIR/peaks/postEGA_peaks_total.bed -b $BASEDIR/peaks/preEGA_peaks_total.bed -v > $BASEDIR/peaks/postEGApeaks_comp.bed

####EGA annotation
bedtools intersect -a /work/mglab/kld/TEanns/*0.1*.bed -b $BASEDIR/peaks/postEGApeaks_comp.bed -f 0.50 -u > $BASEDIR/peaks/postEGA_comp.TEann.txt
awk '{print $4}' $BASEDIR/peaks/postEGA_comp.TEann.txt | sort - | uniq -c | awk '{print $1 "\t" $2}'> $BASEDIR/postEGA_comp_TEcounts.bed

bedtools intersect -a /work/mglab/kld/TEanns/*0.1*.bed -b $BASEDIR/peaks/preEGA_peaks_total.bed -f 0.50 -u > $BASEDIR/peaks/preEGA_peaks_total.TEann.txt
awk '{print $4}' $BASEDIR/peaks/preEGA_peaks_total.TEann.txt | sort - | uniq -c | awk '{print $1 "\t" $2}'> $BASEDIR/preEGA_peaks_total_TEcounts.bed

###peric int
mkdir $BASEDIR/peric/peak_int2
bedtools intersect -a $BASEDIR/peaks/preEGA_peaks_total.bed -b $BASEDIR/peric/centromeres_sloppy2.bed -wa > $BASEDIR/peric/peak_int2/preEGA_peaks.peric.bed
bedtools intersect -a $BASEDIR/peaks/postEGApeaks_comp.bed -b $BASEDIR/peric/centromeres_sloppy2.bed -wa > $BASEDIR/peric/peak_int2/postEGApeaks_comp.peric.bed

for infile in $BASEDIR/peric/null/NULLcoords2*.bed
do
 base=$(basename ${infile} .bed)
  bedtools intersect -a $BASEDIR/peaks/preEGA_peaks_total.bed -b $infile -wa > $BASEDIR/peric/null/$base.preEGA.bed
  bedtools intersect -a $BASEDIR/peaks/postEGApeaks_comp.bed -b $infile -wa > $BASEDIR/peric/null/$base.postEGA.bed
done

for infile in $BASEDIR/peric/null/*.preEGA.bed
do
  wc -l $infile >> $BASEDIR/peric/preEGA_peaks_null_ints.tsv
done

for infile in $BASEDIR/peric/null/*.postEGA.bed
do
  wc -l $infile >> $BASEDIR/peric/postEGA_peaks_null_ints.tsv
done

#####overlap with piRNA targets
mkdir $BASEDIR/piRNA_int
bedtools intersect -a $BASEDIR/peaks/preEGA_peaks_total.bed -b /scratch/kld57880/piRNA_Kaaij_etal/expTEs_2.bed -wa > $BASEDIR/piRNA_int/preEGA_peaks.piRNA.bed
bedtools intersect -a $BASEDIR/peaks/postEGApeaks_comp.bed -b /scratch/kld57880/piRNA_Kaaij_etal/expTEs_2.bed -wa > $BASEDIR/piRNA_int/postEGApeaks_comp.piRNA.bed

##make null set of piRNA coords
mkdir $BASEDIR/piRNA_int/null
for ((i=1;i<=1000;i++));
do
  bedtools shuffle -chrom -noOverlapping -i /scratch/kld57880/piRNA_Kaaij_etal/expTEs_2.bed -excl /scratch/kld57880/piRNA_Kaaij_etal/expTEs_2.bed -g $BASEDIR/genome/chrNameLength.txt > $BASEDIR/piRNA_int/null/piRNAnull_"$i".bed
done

for infile in $BASEDIR/piRNA_int/null/piRNAnull_*.bed
do
 base=$(basename ${infile} .bed)
  bedtools intersect -a $BASEDIR/peaks/preEGA_peaks_total.bed -b $infile -wa > $BASEDIR/piRNA_int/null/$base.preEGA.bed
  bedtools intersect -a $BASEDIR/peaks/postEGApeaks_comp.bed -b $infile -wa > $BASEDIR/piRNA_int/null/$base.postEGA.bed
done

for infile in $BASEDIR/piRNA_int/null/*.preEGA.bed
do
  wc -l $infile >> $BASEDIR/piRNA_int/preEGA_peaks_null_ints.tsv
done

for infile in $BASEDIR/piRNA_int/null/*.postEGA.bed
do
  wc -l $infile >> $BASEDIR/piRNA_int/postEGA_peaks_null_ints.tsv
done

#####
bedtools intersect -a $BASEDIR/peaks/preEGA_peaks_total.bed -b $BASEDIR/peric/centromeres_sloppy2.bed /scratch/kld57880/piRNA_Kaaij_etal/expTEs_2.bed -v > $BASEDIR/peaks/preEGA_K9peaks_noperic_nopiRNA.bed
bedtools intersect -a $BASEDIR/piRNA_int/preEGA_peaks.piRNA.bed -b $BASEDIR/peric/centromeres_sloppy2.bed -v > $BASEDIR/peaks/preEGA_K9peaks_piRNA_only.bed
bedtools intersect -a $BASEDIR/piRNA_int/preEGA_peaks.piRNA.bed -b $BASEDIR/peric/centromeres_sloppy2.bed -wa > $BASEDIR/peaks/preEGA_K9peaks_piRNA_N_peric.bed
bedtools intersect -a $BASEDIR/peric/peak_int2/preEGA_peaks.peric.bed -b $BASEDIR/piRNA_int/preEGA_peaks.piRNA.bed -v > $BASEDIR/peaks/preEGA_K9peaks_peric_only.bed

bedtools intersect -a $BASEDIR/peaks/postEGApeaks_comp.bed -b $BASEDIR/peric/centromeres_sloppy2.bed /scratch/kld57880/piRNA_Kaaij_etal/expTEs_2.bed -v > $BASEDIR/peaks/postEGA_K9peaks_noperic_nopiRNA.bed
bedtools intersect -a $BASEDIR/piRNA_int/postEGApeaks_comp.piRNA.bed -b $BASEDIR/peric/centromeres_sloppy2.bed -v > $BASEDIR/peaks/postEGA_K9peaks_piRNA_only.bed
bedtools intersect -a $BASEDIR/piRNA_int/postEGApeaks_comp.piRNA.bed -b $BASEDIR/peric/centromeres_sloppy2.bed -wa > $BASEDIR/peaks/postEGA_K9peaks_piRNA_N_peric.bed
bedtools intersect -a $BASEDIR/peric/peak_int2/postEGApeaks_comp.peric.bed -b $BASEDIR/piRNA_int/postEGApeaks_comp.piRNA.bed -v > $BASEDIR/peaks/postEGA_K9peaks_peric_only.bed

for infile in $BASEDIR/peaks/preEGA_K9peaks_*.bed
do
  base=$(basename ${infile} .bed)
  bedtools intersect -a /work/mglab/kld/TEanns/*0.1*.bed -b $infile -f 0.50 -u > $BASEDIR/peaks/$base.TEann.txt
done

for infile in $BASEDIR/peaks/preEGA_K9peaks_*.TEann.txt
do
  base=$(basename ${infile} .TEann.txt)
  awk '{print $4}' $infile | sort - | uniq -c | awk '{print $1 "\t" $2}'> $BASEDIR/"$base"_TEcounts.bed
done

for infile in $BASEDIR/peaks/preEGA_K9peaks_*.bed
do
  base=$(basename ${infile} .bed)
  bedtools intersect -a $infile -b /work/mglab/kld/TEanns/*0.1*.bed -F 0.50 -u -wa -wb > $BASEDIR/peaks/$base.peakTEann.txt
done

for infile in $BASEDIR/peaks/preEGA_K9peaks_*.bed
do
  base=$(basename ${infile} .bed)
  bedtools intersect -a $infile -b $BASEDIR/LTRsonly.bed -u > $BASEDIR/peaks/$base.LTRpeaks.bed
  bedtools intersect -a $infile -b $BASEDIR/LINEsonly.bed -u > $BASEDIR/peaks/$base.LINEpeaks.bed
  bedtools intersect -a $infile -b $BASEDIR/DNAsonly.bed -u > $BASEDIR/peaks/$base.DNApeaks.bed
  bedtools intersect -a $infile -b $BASEDIR/SINEsonly.bed -u > $BASEDIR/peaks/$base.SINEpeaks.bed
done
