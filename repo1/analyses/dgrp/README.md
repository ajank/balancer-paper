# DGRP vcf

Downloaded a DGRP version for `dm6` from 

```
https://zenodo.org/record/155396#.WCWWL9w8SuA
```

Then I annotate AF like this (*very slow*):

```
less vcf/dgrp2_dm6_dbSNP.vcf.gz \
| vcf-annotate --fill-AC-AN \
| awk -f ../common/vcf/vcfinfo.awk -W source '/^#/ {print}  \
         !/^#/ {gsub("ALTCOUNT=", "AF=" INFO["AC"]/INFO["AN"] ";ALTCOUNT="); print}' \
> vcf/dgrp2_dm6_dbSNP.annotateAF.vcf
```

Note that I had to add header lines manually afterwards:

```
##INFO=<ID=AF,Number=1,Type=Float,Description="Allele frequency as AC/AN">
##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes">
##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
```

# Balancer vcf

Run `1.get_vcf` to collect bi-allelic SNPs and indels.


# `bcftools stats` to get the numbers and overlaps

```
mkdir plots 2> /dev/null
for type in bal-spec vrg-spec common; do
    bcftools stats vcf/dgrp2.dm6.vcf.gz vcf/${type}.vcf.gz > plots/${type}.chk
    plot-vcfstats -p plots/${type} plots/${type}.chk
done
```

# DEL_comparison

Deletions were unfortunately not included in the lifted-over file `vcf/dgrp2.dm6.bcf`.
So I had to start off the dm3 version (`../../../data/external/DGRP_freeze2/dgrp2.bcf`),
select the *287,376* simple deletions via `awk`

```
bcftools view -H ../../../data/external/DGRP_freeze2/dgrp2.bcf \
    | awk 'length($4)>1 && length($5)==1 {OFS="\t"; print "chr" $1, $2-1, $2-1+length($4)}' \
    > dels.txt
```

and then lift over this BED file using UCSC genome browser. The conversion 
failed only for 3 entries ==> `hglft_genome_5eb3_ffb3b0.bed`.

After that I ran `plots.R` to compare these Deletions to our balancer deletions.

### Same for DUPs

Duplications were not specifically marked within the VCF. So I took the category "INS"
and subsetted for INS > 200. Assuming these would be duplications (I only checked the first one
by Blatting it) I created a BED file like this:

```
bcftools view -H ../../../data/external/Dfreeze2/dgrp2.bcf | awk '$3~/INS/ && length($4)==1 && length($5)>100 {print "chr" $1 "\t" $2-length($5) "\t" $2}' > dups.txt
```

and then lifted those over to dm6 --> `hglft_genome_52e6_1d0620.bed`.

Then I ran `plots_dup.R`.

# Miller *et al* regions

See `miller_regions/analysis.sh` and `miller_regions/analysis.R` for documentation

