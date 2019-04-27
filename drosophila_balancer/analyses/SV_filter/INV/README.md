# INV filtering

Actually I cannot do much to filter them properly.

```
vcfinfo=/g/korbel/shared/projects/drosophila_balancer/data/variants/vcfinfo.awk
bcftools view /g/korbel/shared/projects/drosophila_balancer/data/variants/SVs3/delly0.7.5.INV.bcf \
    | awk -f $vcfinfo -W source '/^#/ || \
		$7~/PASS/ && $1~/^chr[23][LR]$/ && SAMPLE[2,"GT"] != "./." && SAMPLE[1,"GT"] == "0/1" && INFO["MAPQ"]>=40' \
	| bcftools view -Ob - \
	> delly0.7.5.INV.filter1.bcf
```

plot

```
bcftools view -H delly0.7.5.INV.filter1.bcf \
	| awk -f $vcfinfo -W source '{OFS="\t"; split($8,precise,";"); \
		print $1, $2, INFO["END"], $3, precise[1], SAMPLE[1,"GT"], SAMPLE[2,"GT"]}' \
	| Rscript 2.plot.R
```

