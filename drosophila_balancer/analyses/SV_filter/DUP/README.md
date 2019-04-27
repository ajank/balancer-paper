# New approach to call tandem DUPs

Based on Delly 0.7.5 calls

## Pre-filtering:

```
vcfinfo=/g/korbel/shared/projects/drosophila_balancer/data/variants/vcfinfo.awk
bcftools view /g/korbel/shared/projects/drosophila_balancer/data/variants/SVs3/delly0.7.5.DUP.bcf \
	| awk -f $vcfinfo -W source '/^#/ || $7~/PASS/ && $1~/^chr[23][LR]$/ && \
		SAMPLE[1,"GT"] == "0/1" && SAMPLE[2,"GT"]!="./."' \
	| bcftools view -Ob - > delly0.7.5.DUP.filter1.bcf
```

## rd-ratio

The following lines plot size vs. read-depth ratio, stratified by genotypes:

```
bcftools view -H delly0.7.5.DUP.filter1.bcf \
	| awk -f $vcfinfo -W source '!/^#/ {OFS="\t"; print $3, \
		SAMPLE[1,"GT"] "_" SAMPLE[2,"GT"], \
		INFO["END"] - $2, \
		SAMPLE[1,"RC"]/(SAMPLE[1,"RCR"]+SAMPLE[1,"RCL"]), \
		SAMPLE[2,"RC"]/(SAMPLE[2,"RCR"]+SAMPLE[2,"RCL"]) }' \
	| Rscript 2.plot.R delly0.7.5.DUP.filter1.pdf
```

From this plot we estimate that all calls > 100kb must be wrong.

## Manual annotation based on bi-allelic frequency 

Plots (for all <100kb) are generated in `analyses/SV_validation/DUP/bi_allelic_freq/new_plots`.

Here are the results from the manual classification

ID       | classification
-------- | ---
???      | false

## **Final** tandem DUP

Subset bcf based on this list of IDs

```
(bcftools view -h delly0.7.5.DUP.filter1.bcf; bcftools view delly0.7.5.DUP.filter1.bcf | grep -f <(awk '$2~/TRUE|likely/' 3.manual_classification.txt | cut -f1)) | bcftools view -Ob - > delly0.7.5.DUP.final.bcf
```

Then I can plot the rd-ratio plot again:

```
bcftools view -H delly0.7.5.DUP.final.bcf \
    | awk -f $vcfinfo -W source '!/^#/ {OFS="\t"; print $3, \
        SAMPLE[1,"GT"] "_" SAMPLE[2,"GT"], \
        INFO["END"] - $2, \
        SAMPLE[1,"RC"]/(SAMPLE[1,"RCR"]+SAMPLE[1,"RCL"]), \
        SAMPLE[2,"RC"]/(SAMPLE[2,"RCR"]+SAMPLE[2,"RCL"]) }' \
    | Rscript 2.plot.R delly0.7.5.DUP.final.pdf
```

