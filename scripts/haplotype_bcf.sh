#!/bin/bash

# $1: gff file containing the locations of the cds of the gene of interest
# $2: BCF file with phased variants
# $3: output directory for temp files
# $3: ID of gene of interest

echo $1
echo $2
echo $3
echo $4

grep $4 $1 | awk '{print $1 ":" $4 "-" $5}' | xargs -i -n 1000 echo | sed 's/ /,/g' | xargs -i -n 1 bcftools view -Oz -r {} $2 > $3/$4.vcf.gz
