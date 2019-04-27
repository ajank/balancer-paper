## Exons

```
GFF=../../../data/dm6_annotation/dm6-all-filtered-r6.05.UCSC_names.genes.gff3
awk '$3~/exon/ {OFS="\t"; print $1, $4, $5, substr($9,9,11)}' $GFF \
    | sort -k1,1 -k2,2n -k3,3 \
    | uniq \
    > exons.bed
```



### Function to merge (flatten) overlapping exons

I might use the flattened exon for a quick ASE annotation as shown in the next section

```
function flatten_by_gene 
{
	dir=$(mktemp -d)
	# split into many files
	sort -k4,4 $1 | awk '{print $0 > "'$dir'/"$4}'
	# merge each file into non-overlapping intervals
	for x in $dir/*; do 
		y=$(basename $x); 
		sort -k1,1 -k2,2n $x | bedtools merge | awk '{print $0 "\t'$y'"}'; 
	done
	rm -r $dir
}
```



### Overlap with deletions

First, flatten exons (see above).

Then, the total exonic size per gene 

```
awk '{print $3-$2+1 "\t" $4}' exons_flattened.bed | sort -k2,2 | q -t "select c2,sum(c1) from - group by c2" > exon_sizes.txt
```

Next, calculate overlap of exons with DEL, annotate gene length and get relative deletion length:

```
bedtools intersect -a exons_flattened.bed -b ../deletions/min20.merged.bal-spec.bed -wo \
	| cut -f4,8 \
	| sort -k1,1 \
	| join - exon_sizes.txt \
	| awk '{print $1 "\t" $2*100/$3}' \
	> gene_del_ovlp.bal-spec.txt
```

Repeat the same for vrg-spec DELs.



### Get ASE exons to be loaded as a track in IGB

```
awk 'NR>1 && $7<0.05 {print $1}' ../../ase/deseq/DESeq.N1_heads.standardFormat.txt \
	| sort | grep -f - exons_flattened.bed \
	> exons_flattened.ASE.heads.bed
awk 'NR>1 && $7<0.05 {print $1}' ../../ase/deseq/DESeq.N1_6-8h.standardFormat.txt \
	| sort | grep -f - exons_flattened.bed \
	> exons_flattened.ASE.6-8h.bed
```

