#!/bin/bash

bedtools getfasta -fi $2 -bed $3/$1.bed > $4/$1_h$5.fasta
samtools faidx $4/$1_h$5.fasta
