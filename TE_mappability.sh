#!/bin/bash
module load Anaconda3/5.0.1
source activate genmap_kd

mkdir $BASEDIR/map_35
genmap map -K 35 -E 3 -T 24 -I $BASEDIR/genmap_index -O $BASEDIR/map_35 -t -w -bg

mkdir $BASEDIR/map_75
genmap map -K 75 -E 3 -T 24 -I $BASEDIR/genmap_index -O $BASEDIR/map_75 -t -w -bg

ml ucsc

wigToBigWig $BASEDIR/map_35/danio_refseq.genmap.wig /scratch/kld57880/MO_CnR/genome/chrNameLength.txt $BASEDIR/map_35bp.bw
wigToBigWig $BASEDIR/map_75/danio_refseq.genmap.wig /scratch/kld57880/MO_CnR/genome/chrNameLength.txt $BASEDIR/map_75bp.bw

###going to now compute the mappability of all TEs in the genome
module load deepTools/3.3.1-intel-2019b-Python-3.7.4

computeMatrix scale-regions -S $BASEDIR/map_35bp.bw -R $BASEDIR/TEfiles_Chang_etal/TEann_total.bed -p 24 -o $BASEDIR/TE_35bp_mappability.gz --outFileNameMatrix $BASEDIR/TE_35bp_mappability.tab --outFileSortedRegions $BASEDIR/TE_35bp_map_regions.bed
computeMatrix scale-regions -S $BASEDIR/map_75bp.bw -R $BASEDIR/TEfiles_Chang_etal/TEann_total.bed -p 24 -o $BASEDIR/TE_75bp_mappability.gz --outFileNameMatrix $BASEDIR/TE_75bp_mappability.tab --outFileSortedRegions $BASEDIR/TE_75bp_map_regions.bed

###in R, averaged each row across columns to give average mappability score to each TE, length normalized
###then combined back with the coordinates to make bed file with chr start end map score
###now need to go back and re-add the annotation to that
module load BEDTools/2.29.2-GCC-8.3.0

bedtools intersect -a $BASEDIR/TEfiles_Chang_etal/TEann.gtf -b $BASEDIR/TE_75bp_coords_score.bed -f 0.95 -wa -wb | sort | uniq > $BASEDIR/TEann_75mappabilityScore.bed

sed 's:"::g' $BASEDIR/TEann_75mappabilityScore.bed | sed 's:gene_id::g' | sed 's:transcript_id::g' | sed 's:family_id::g' | sed 's:class_id::g' | sed 's:Cluster::g' | sed 's:TE_classification::g' | sed 's:;::g' > $BASEDIR/TEann_75mappabilityScore1.bed
awk '{print $1 "\t" $4 "\t" $5 "\t" $9 "\t" $13 "\t" $14 "\t" $19 "\t" $20 "\t" $21 "\t" $22}' $BASEDIR/TEann_75mappabilityScore1.bed > $BASEDIR/TEann_75mapScore_neat.bed
