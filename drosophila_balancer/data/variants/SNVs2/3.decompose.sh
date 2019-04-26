#/bin/bash
vt decompose_blocksub wgs.freebayes-k.filter.norm.vcf.gz \
    | vt uniq - \
	| vcfsnps \
	| vcfbiallelic \
	| bgzip > wgs.freebayes-k.filter.norm.decomposed.onlySNPs.vcf.gz
tabix wgs.freebayes-k.filter.norm.decomposed.onlySNPs.vcf.gz
