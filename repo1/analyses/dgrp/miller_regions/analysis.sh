#!/bin/bash

echo "Balancer SNPs to table..."
less ../../../data/variants/SNVs2/wgs.freebayes-k.filter.norm.decomposed.onlySNPs.vcf.gz \
    | awk '
        BEGIN{OFS="\t";
            print "chrom","pos","ref","alt","gt"
        } 
        !/^#/ && $1~/chr[23X][LR]?/ \
        && $11~/0\/0|1\/1/ \
        && $10~/0\/1|1\/1/ \
        && !($10~/1\/1/ && $11~/0\/0/) {
            print $1,$2,$4,$5,substr($10,0,3) "_" substr($11,0,3)
        }' \
    > SNPs_wgs.txt
echo "Quick number check:"
awk 'NR>1 {print substr($1,0,4) "\t" $5}' SNPs_wgs.txt \
    | sort \
    | uniq -c


echo ""
echo "DGRP SNPs to table"
bcftools view -H ../vcf/dgrp2.dm6.bcf \
    | awk 'BEGIN  {
            OFS="\t";
            print "chrom","pos","ref","alt", "af"
        } 
        length($4)==1 && length($5)==1 {
            print $1, $2, $4, $5, gensub(/AF=([0-9\.]+);.*/, "\\1", 1, $8);
        }' \
    > SNPs_dgrp.txt

Rscript analysis.R

