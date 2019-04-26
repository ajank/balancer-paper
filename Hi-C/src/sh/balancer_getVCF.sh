#!/bin/bash

# SAMPLE[1,"GT"] is cross
# SAMPLE[2,"GT"] is vrg

#
#  /g/korbel/shared/projects/drosophila_balancer/analyses/readSeparation/getVCF/variants.new.vcf.gz was generated as follows:
#

# zcat /g/korbel/shared/projects/drosophila_balancer/data/variants/SNVs2/wgs.freebayes-k.filter.norm.decomposed.onlySNPs.vcf.gz \
# 	| awk -f /g/korbel/shared/projects/drosophila_balancer/data/variants/vcfinfo.awk \
# 		-W source '$0~/^#/ {print} \
# 				   $0!~/^#/ && SAMPLE[2,"GT"]~/0\/0/ && SAMPLE[1,"GT"]~/0\/1/ \
# 				   		{print} \
# 				   $0!~/^#/ && SAMPLE[2,"GT"]~/1\/1/ && SAMPLE[1,"GT"]~/0\/1/ \
# 						{OFS="\t"; print $1,$2,$3,$5,$4,$6,$7,$8}' \
# 	| cut -f 1-8 \
# 	| bgzip \
# 	> variants.new.vcf.gz
# tabix variants.new.vcf.gz

#
#  extract common SNVs
#

zcat /g/korbel/shared/projects/drosophila_balancer/data/variants/SNVs2/wgs.freebayes-k.filter.norm.decomposed.onlySNPs.vcf.gz \
	| awk -f /g/korbel/shared/projects/drosophila_balancer/data/variants/vcfinfo.awk \
		-W source '$0~/^#/ {print} \
				   $0!~/^#/ && SAMPLE[2,"GT"]~/1\/1/ && SAMPLE[1,"GT"]~/1\/1/ \
            {OFS="\t"; print $1,$2,$3,$5,$5,$6,$7,$8}' \
	| cut -f 1-8 \
	| bgzip \
	> analysis/balancer/variants-common.new.vcf.gz
tabix analysis/balancer/variants-common.new.vcf.gz

#
#  convert to .tab files for reading in R
#

infasta="/g/furlong/genome/D.melanogaster/Dm6/fasta/dm6.UCSC.noMask.fa"

src/sh/balancer_vcf_to_tab.sh /g/korbel/shared/projects/drosophila_balancer/analyses/readSeparation/getVCF/variants.new.vcf.gz \
  $infasta analysis/balancer/SNVs.tab

src/sh/balancer_vcf_to_tab.sh analysis/balancer/variants-common.new.vcf.gz \
  $infasta analysis/balancer/SNVs_common.tab
