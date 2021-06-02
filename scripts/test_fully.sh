#!/bin/bash

view=$(grep $1 ref_genome/ITAG4.0_cds_models.gff | awk '{print $1 ":" $4 "-" $5}' | xargs -i -n 100 echo | sed 's/ /,/g')
echo $view
bcftools view -Ov -r $view ref_genome/ITAG_cds_models.gff | grep $'^SL' > $1_tmp

phased=$(bcftools view -p $1_tmp | grep $'^SL' | wc -l)
tot=$(bcftools view $1_tmp | grep $'^SL' | wc -l)
state=$(bcftools view -g het $1_tmp | grep $'^SL' | cut -f10 | cut -d ":" -f3 | sort | uniq | wc -l)

echo 'phased: '$phased', tot: '$tot', state: '$state

#awk -v phased="$phased" -v tot="$tot" -v prot="$1" 'BEGIN {if (tot != 0 && phased/tot == 1) print prot}' >> $2
#awk -v phased="$phased" -v tot="$tot" -v state="$state" -v prot="$1" 'BEGIN {if (tot != 0 && phased/tot == 1 && state == 1) print prot}' >> $3
#awk -v phased="$phased" -v tot="$tot" -v state="$state" -v prot="$1" 'BEGIN {if (tot - phased == 1) print prot}' >> $4

#rm $1_tmp

