#!/bin/bash

echo "bcftools concat somatic variants..."
echo "cross" > header
echo "vrg"  >> header

#bcftools-1.2 concat -a MP.delly071.???.vcf.gz | bgzip-1.2.1 > MP.delly071.X.vcf.gz

for x in ??.delly07?.???.vcf.gz;
do
	echo "Renaming $x"
	bcftools-1.2 reheader -s header $x > tmp.vcf.gz
	if [ "$(zcat $x | wc -l)" -eq "$(zcat tmp.vcf.gz | wc -l)" ]; then
		echo "$(zcat tmp.vcf.gz | wc -l) lines"
		echo "$(zcat $x | wc -l) lines"
		mv tmp.vcf.gz $x
		tabix-1.2.1 -p vcf $x
	else
		echo "Problem with reheader"
		exit 1
	fi
done

rm header
