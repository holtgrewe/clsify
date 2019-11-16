#!/bin/bash

REGIONS="16S 50S"
TSV=Grimm2018

for region in $REGIONS; do
    accs=$(awk -v region=$region -F $'\t' \
             'BEGIN { OFS = FS } (NR > 1 && $3 == region) { print $4 }' \
             $TSV.tsv)
    for acc in $accs; do
        region=$(grep -w $acc $TSV.tsv | cut -f 3)
        mkdir -p $TSV/$region

        curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=${acc}&rettype=fasta&retmode=txt" \
        | grep -v '^$' \
        > $TSV/$region/$acc.fasta
    done
done

for x in $TSV/*; do
    cat $x/*.fasta \
    > ${x/\//_}.fasta
done
