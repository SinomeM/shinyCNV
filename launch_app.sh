#!/bin/bash

#wkdir=$1
#cnvs=$2
#samples=$3
#snps=$4

wkdir=/users/sm/Documents/GitHub/shinyCNV
cnvs=/users/sm/Documents/GitHub/shinyCNV/data/cnvs.txt
samples=/users/sm/Documents/GitHub/shinyCNV/data/samples.txt
snps=/users/sm/Documents/GitHub/shinyCNV/data/hd_1kG_hg19.snppos.filtered.gz

cd $wkdir

Rscript app.R $wkdir $cnvs $samples $snps