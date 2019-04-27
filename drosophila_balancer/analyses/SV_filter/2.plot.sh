# Plot the filtered DELs and DUPs
# (inversion filtering has not yet yielded acceptable results)

cat DEL/MP.cnv.vcf DEL/PE.cnv.vcf DUP/MP.dup.vcf DUP/PE.dup.vcf \
	| grep -v '^#' \
	| sort -k1,1 -k2,2n \
	> tmp.vcf 

dm6vis tmp.vcf dup_del
rm tmp.vcf

