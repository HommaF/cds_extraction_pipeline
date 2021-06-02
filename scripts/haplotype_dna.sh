#!/bin/bash

# $1: ID of gene of interest
# $2: gff file containing the locations of the cds of the gene of interest
# $3: BCF file with phased variants
# $4: haplotype 1 or 2
# $5: state of phasing of gene of interest
# $6: reference genome

grep $1 ref_genome/ITAG4.0_cds_models.gff | awk '{print $1 ":" $4 "-" $5}' | xargs -i -n 1 bcftools view -Ov -r {} $3 > $1_$4$5_tmp.vcf
bcftools consensus -H $5 -f $6 $1_$4$5_tmp.vcf > $7/$4_$1_H$5.fasta


