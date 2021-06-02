#!/bin/bash

if [[ $1 != Solyc10g005125.1.1 ]]; then
	grep $1 $2 | awk '{print $1 "\t" $4-100 "\t" $5+100}' > $3/$1.bed
else
	grep $1 $2 | awk '{print $1 "\t" $4-100 "\t" $5-10}' > $3/$1.bed
fi
