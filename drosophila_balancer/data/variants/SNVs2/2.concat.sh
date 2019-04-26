#!/bin/bash
bcftools-1.2 concat wgs.freebayes-k.filter.norm.chr*.vcf.gz \
	| bgzip-1.2.1 > wgs.freebayes-k.filter.norm.vcf.gz
zcat wgs.freebayes-k.filter.norm.vcf.gz | vcfbiallelic | vcfsnps | bgzip-1.2.1 > wgs.freebayes-k.filter.norm.SNVonly.vcf.gz

tabix-1.2.1 wgs.freebayes-k.filter.norm.vcf.gz &
tabix-1.2.1 wgs.freebayes-k.filter.norm.SNVonly.vcf.gz &

for x in wgs.freebayes-k.filter.norm.chr*.vcf.gz*;
	do
	mv $x _delete_$x
	done

