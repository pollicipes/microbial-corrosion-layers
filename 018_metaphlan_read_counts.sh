#!/bin/bash
# ------------------------------------------------------------------
# [Author] Title
#          Description
# ------------------------------------------------------------------

WD="/home/ngp704/data/iron_reducing_bacteria/metaphlan/out/"

# ${WD}/${sample}.reads_profile.txt;
cd $WD;


for i in `ls *.profile.txt`; do smp=${i%.txt}; cx=$(grep 'reads' ${i} | cut -f1 -d' ' | cut -f2 -d'#'); echo -e ${smp}"\t"${cx}; done > reads_sequenced.txt

# Get the read counts 
for i in `ls *.reads_profile.txt`; do
    smp=${i%.reads_profile.txt};
    cx=$(grep 'estimated_reads_mapped_to_known_clades' ${i} | cut -f2 -d':');
    echo -e ${smp}"\t"${cx};
done > read_counts.txt # reads_sequenced.txt # THE OUTPUT TO SAVE THE TOTAL READS PER LIBRARY
