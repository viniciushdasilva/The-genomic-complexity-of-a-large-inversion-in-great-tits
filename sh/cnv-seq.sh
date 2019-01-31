#!/bin/bash

FILES=/lustre/nobackup/WUR/ABGC/shared/great_tit/BAM_30birds/*.bam
for f in $FILES;
do

test=$f
test_prefix=`echo $test | cut -d'.' -f2 `
ref_prefix=`echo $ref | rev | cut -d'/' -f1 | rev`

echo $test_prefix

## test
samtools view -F 4 $f |perl -lane 'print "$F[2]\t$F[3]"' >$test_prefix.hits
## control
samtools view -F 4 $2 7 |perl -lane 'print "$F[2]\t$F[3]"' >$ref_prefix.hits

## Run cnv-seq
cnv-seq.pl --test $test_prefix.hits --ref control.hits --genome-size 71365269

## Plot CNV
Rscript plotCNV.Rscript $test_prefix.hits*$ref_prefix.hits*.cnv $title

done
