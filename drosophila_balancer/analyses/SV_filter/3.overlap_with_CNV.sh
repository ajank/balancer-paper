#!/bin/bash

CNV=../../data/variants/CNVs/

function delly2bed {
	 vcf-sort | grep -v '^#' \
			  | cut -f 1,2,3,8 \
			  | sed 's/\t[^\t]*;END=/\t/' \
			  | sed 's/;.*$//' \
			  | awk '{print $1"\t"$2"\t"$4"\t"$3;}'
}

function cnv2bed {
	awk '{OFS="\t"; split($2,c,":"); split(c[2],p,"-"); print c[1], p[1], p[2], ($1~/deletion/?"<DEL>":"<DUP>"), $4, $5, $6, $7, $8, $9, $10}' 
}

echo "Deletions: "
bedtools intersect -r -f 0.5 -wb \
	-a <(cat DEL/??.cnv.vcf | delly2bed) \
	-b <(cat $CNV/*.calls | grep deletion | cnv2bed) | cut -f-8
echo "Duplications: "
bedtools intersect -r -f 0.5 -wb \
	-a <(cat DUP/??.dup.vcf | delly2bed) \
	-b <(cat $CNV/*.calls | grep duplication | cnv2bed) | cut -f-8

