#!/bin/bash
module load BLAST+/2.12.0-gompi-2020b
makeblastdb -in $BASEDIR/danio_refseq.fa -dbtype nucl

blastn -query /work/mglab/kld/main_sat1_sequence_1992.fasta -task blastn -db $BASEDIR/danio_refseq.fa -outfmt 6

columns: Qaccession Saccession %IdenticalMatches Length #Mismatches #Gaps Qstart Qend Sstart Send eValue Bitscore

cat sat1_blast.txt | awk '{print $2 "\t" $9 "\t" $10}' > sat1.bed
cat sat1_blast_trimmed.txt | awk '{print $2 "\t" $9 "\t" $10}' > sat1_trimmed.bed

####manually created bed file with rough calls of centromeres based on where the highest density of SAT1 repeats were called
###to broaded these regions I'm going to use bedtools slop
ml BEDTools

for ((i=1;i<=1000;i++));
do
  bedtools shuffle -chrom -noOverlapping -i $BASEDIR/peric/centromeres_sloppy2.bed -excl $BASEDIR/peric/centromeres_sloppy2.bed -g $BASEDIR/genome/chrNameLength.txt > $BASEDIR/peric/null/NULLcoords2_"$i".bed
done
