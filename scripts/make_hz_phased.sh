#!/bin/bash


bcftools view -Ov $1 | grep $'^SL' | cut -f1,2,3,4,5,6,7,8,9 > $1_first_nine.vcf
bcftools view -Ov $1 | grep $'^SL' | cut -f10 | sed 's/^1\/1/1\|1/g' > $1_tenth.vcf
bcftools view -Ov $1 | head -100 | grep $'#' > $1_headers.vcf

paste $1_first_nine.vcf $1_tenth.vcf >> $1_headers.vcf
bcftools view -Ob $1_headers.vcf > $2

bcftools index $2

rm $1_headers.vcf
rm $1_first_nine.vcf
rm $1_tenth.vcf

