#!/bin/bash

# SAMPLE[1] is cross
# SAMPLE[2] is vrg

zcat ../../../data/variants/SNVs2/wgs.freebayes-k.filter.norm.SNVonly.vcf.gz \
	| awk -f ../../../data/variants/vcfinfo.awk \
		-W source '$0~/^#/ {print} \
				   $0!~/^#/ && SAMPLE[2,"GT"]~/0\/0/ && SAMPLE[1,"GT"]~/0\/1/ \
				   		{print} \
				   $0!~/^#/ && SAMPLE[2,"GT"]~/1\/1/ && SAMPLE[1,"GT"]~/0\/1/ \
						{OFS="\t"; print $1,$2,$3,$5,$4,$6,$7,$8}' \
	| cut -f 1-8 \
	| bgzip-1.2.1 \
	> variants.vcf.gz
tabix-1.2.1 variants.vcf.gz


zcat ../../../data/variants/SNVs2/wgs.freebayes-k.filter.norm.decomposed.onlySNPs.vcf.gz \
	| awk -f ../../../data/variants/vcfinfo.awk \
		-W source '$0~/^#/ {print} \
				   $0!~/^#/ && SAMPLE[2,"GT"]~/0\/0/ && SAMPLE[1,"GT"]~/0\/1/ \
				   		{print} \
				   $0!~/^#/ && SAMPLE[2,"GT"]~/1\/1/ && SAMPLE[1,"GT"]~/0\/1/ \
						{OFS="\t"; print $1,$2,$3,$5,$4,$6,$7,$8}' \
	| cut -f 1-8 \
	| bgzip-1.2.1 \
	> variants.new.vcf.gz
tabix-1.2.1 variants.new.vcf.gz
