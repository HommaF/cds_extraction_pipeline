#!/bin/bash

view=$(grep $1 $2 | awk '{print $1 ":" $4 "-" $5}' | xargs -i -n 100 echo | sed 's/ /,/g')

bcftools view -Ov -r $view $3 > $7_$1_tmp.bcf

phased=$(bcftools view -p $7_$1_tmp.bcf | grep $'^SL' | wc -l)
tot=$(bcftools view $7_$1_tmp.bcf | grep $'^SL' | wc -l)
state=$(bcftools view -p -g het $7_$1_tmp.bcf | grep $'^SL' | cut -f10 | cut -d ":" -f3 | sort | uniq | wc -l)


awk -v phased="$phased" -v tot="$tot" -v prot="$1" 'BEGIN {if (tot != 0 && phased/tot == 1) print prot}' >> $4
awk -v phased="$phased" -v tot="$tot" -v state="$state" -v prot="$1" 'BEGIN {if (tot != 0 && phased/tot == 1 && state == 1) print prot}' >> $5
awk -v phased="$phased" -v tot="$tot" -v state="$state" -v prot="$1" 'BEGIN {if (tot - phased == 1) print prot}' >> $6
awk -v phased="$phased" -v tot="$tot" -v prot="$1" 'BEGIN {if (tot == 0) print prot}' >> $8
