#!/bin/bash

exonerate --model protein2genome --bestn 1 --refine full -t $2/$1_h$7.fasta -q $3/ITAG4.0_proteins.id_$1_nostop.fasta --showalignment False --showvulgar False --showtargetgff True | grep $'\tcds\t' > $4/$1_h$7.gff

mysize=$(find "$4/$1_h$7.gff" -printf "%s")

if [[ $mysize > 0 ]]; then

	bedtools getfasta -fi $2/$1_h$7.fasta -bed $4/$1_h$7.gff -s > $4/$1_h$7_ms.fasta
	union -filter -sequence $4/$1_h$7_ms.fasta > $5/tmp_$6_$1_h$7.fasta
	seqtk rename $5/tmp_$6_$1_h$7.fasta $6_$1_h$7 > $5/$6_$1_h$7.fasta
	rm $5/tmp_$6_$1_h$7.fasta

else
	echo $1_h$7 >> $5/not_found.txt
fi

