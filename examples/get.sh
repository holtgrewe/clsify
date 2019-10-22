#!/bin/bash

SAMPLES="
MF041971-CA_SanJoaquinValley_Tm4_16S_LsoF
MF041968-CA_SanJoaquinValley_Tm3_16S_LsoF
MF041970-CA_SanJoaquinValley_Tm4_50S_CL514F
MF041967-CA_SanJoaquinValley_Tm3_50S_CL514F
"

for sample in $SAMPLES; do
    acc=$(echo $sample | cut -d - -f 1)
    file=$(echo $sample | cut -d - -f 2-)
    curl https://www.ncbi.nlm.nih.gov/search/api/sequence/$acc/?report=fasta \
    > $file.fasta
done
