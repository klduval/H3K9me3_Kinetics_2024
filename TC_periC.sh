#!/bin/bash

module load BLAST+/2.12.0-gompi-2020b
makeblastdb -in $BASEDIR/danio_refseq.fa -dbtype nucl

blastn -query /work/mglab/kld/main_sat1_sequence_1992.fasta -task blastn -db $BASEDIR/danio_refseq.fa -outfmt 6

#columns: Qaccession Saccession %IdenticalMatches Length #Mismatches #Gaps Qstart Qend Sstart Send eValue Bitscore

cat sat1_blast.txt | awk '{print $2 "\t" $9 "\t" $10}' > sat1.bed
cat sat1_blast_trimmed.txt | awk '{print $2 "\t" $9 "\t" $10}' > sat1_trimmed.bed

####manually created bed file with rough calls of centromeres based on where the highest density of SAT1 repeats were called
###this failed for chrs 1 and 14, also could be more convincing for chrs 6, 15, and 21
###to broaded these regions I'm going to use bedtools slop

ml BEDTools

bedtools slop -b 1.5 -pct -i $BASEDIR/centromere_coords_raw.bed -g $BASEDIR/genome/chrNameLength.txt > $BASEDIR/centromeres_sloppy.bed

#and then to make the control regions
bedtools makewindows -b $BASEDIR/centromeres_sloppy.bed -2 > $BASEDIR/centromeres_split.bed
bedtools shuffle -chrom -noOverlapping -i $BASEDIR/centromeres_sloppy.bed -excl $BASEDIR/centromeres_sloppy.bed -g $BASEDIR/genome/chrNameLength.txt > $BASEDIR/centromereNULLcoords1.bed
bedtools shuffle -chrom -noOverlapping -i $BASEDIR/centromeres_sloppy.bed -excl $BASEDIR/centromeres_sloppy.bed -g $BASEDIR/genome/chrNameLength.txt > $BASEDIR/centromereNULLcoords2.bed
bedtools shuffle -chrom -noOverlapping -i $BASEDIR/centromeres_sloppy.bed -excl $BASEDIR/centromeres_sloppy.bed -g $BASEDIR/genome/chrNameLength.txt > $BASEDIR/centromereNULLcoords3.bed
