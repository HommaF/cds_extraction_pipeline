#!/bin/bash


cat $1 | xargs -i -n 1 -P 4 scripts/phased_genes.sh {} $2 $3 $4 $5 $6
