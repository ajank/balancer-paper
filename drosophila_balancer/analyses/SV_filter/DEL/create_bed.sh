function toBed {
    awk -f ../vcfinfo.awk -W source '!/^#/ {OFS="\t"; print $1, $2-1, INFO["END"], SAMPLE[1,"GT"] "_" SAMPLE[2,"GT"] "_" $3}'
}

for x in DEL.50_160.highConf.vcf DEL.above160.highConf.vcf DEL.above160.lowConf.vcf DEL.below50.highConf.vcf
do
	echo $x
	cat $x | toBed > ${x%.vcf}.bed
done

cat DEL.50_160.highConf.bed DEL.above160.highConf.bed DEL.above160.lowConf.bed DEL.below50.highConf.bed | sort -k1,1 -k2,2n > DEL.all.bed
