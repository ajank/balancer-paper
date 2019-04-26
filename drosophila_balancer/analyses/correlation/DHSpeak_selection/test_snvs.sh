SNP_VCF=../../../data/variants/SNVs2/wgs.freebayes-k.filter.norm.SNVonly.vcf.gz
SNPS_BAL=$(mktemp --suffix ".bed")
SNPS_VRG=$(mktemp --suffix ".bed")
NUM_SHUFFLES=100

echo "Obtaining SNPs"
zcat $SNP_VCF | awk -f ../../common/vcf/vcfinfo.awk -W source '/^chr[23][LR]\t/ && SAMPLE[1,"GT"]=="0/1" && SAMPLE[2,"GT"]=="0/0" {print $1 "\t" $2-1 "\t" $2}' > $SNPS_BAL
zcat $SNP_VCF | awk -f ../../common/vcf/vcfinfo.awk -W source '/^chr[23][LR]\t/ && SAMPLE[1,"GT"]=="0/1" && SAMPLE[2,"GT"]=="1/1" {print $1 "\t" $2-1 "\t" $2}' > $SNPS_VRG

get_genome_cov() {
    bedtools genomecov -g genome.chr23.fai -i $1 -max 1 | awk '$1=="genome" && $2==0 {print $5}'
}

snps_vs_background() {
    # $1 = SNP bed
    # $2 = DHS bed
    echo -e "real\t$(bedtools intersect -a ${1} -b ${2} | wc -l)\t$(get_genome_cov ${2})"
    TMP=$(mktemp --suffix ".bed")
    for i in $(seq 1 $NUM_SHUFFLES)
    do
        bedtools shuffle -g genome.chr23.fai -i ${2} | sort -k1,1 -k2,2n > $TMP
        echo -e "rand\t$(bedtools intersect -a ${1} -b $TMP | wc -l)\t$(get_genome_cov $TMP)"
    done
}

ggplot_string='d[, V4 := 1-V4] # fraction of the genome covered by DHS regions
    norm_factors = d[V2=="real",.(V1,norm = V4)]
    d = merge(d, norm_factors, by=c("V1"))
    d[, V3 := V3/V4*norm] # normalize by how many bases are covered by the shuffled intervals
    d[, c("V4","norm"):= .(NULL,NULL)]
    d[, rank := frank(V3), by=V1]
    pval = d[, .(pval = (rank[V2=="real"]+0)/(length(rank[V2=="rand"])+1)), by=V1]
    ggplot(d[V2=="rand"]) + aes(V3) + 
        geom_histogram(binwidth=500) +
        geom_vline(data=d[V2=="real"], aes(xintercept=V3), col="red") + 
        facet_grid(V1~.) + 
        xlab("Total number of SNVs in intervals") +
        geom_label(x=Inf, y=Inf, aes(label=round(pval,3)), data=pval, hjust=1, vjust=1)'

echo "Counting SNPs in DHS peaks ($NUM_SHUFFLES times)"
DHS=$(mktemp --suffix ".bed")
grep -P '^chr[23][LR]' ../../tracks/enhancer/macs2_idr0.1_merge.dm6.bed | sort -k1,1 -k2,2n > $DHS
(snps_vs_background $SNPS_BAL $DHS | awk '{print "bal\t" $0}';
 snps_vs_background $SNPS_VRG $DHS | awk '{print "vrg\t" $0}') \
    | ggplot -o SNPs_DHS.pdf -W 6 -H 4 \
    "$ggplot_string + ggtitle(\"SNVs in DHS peaks and randomly matched intervals\")" \
    stdin
rm $DHS


echo "Counting SNPs in DHS summits ($NUM_SHUFFLES times)"
SUMMITS=$(mktemp --suffix ".bed")
grep -P '^chr[23][LR]' ../../tracks/enhancer/macs_summits_p50_slop50_merge25_clean.htseq.dm6.bed | sort -k1,1 -k2,2n > $SUMMITS
(snps_vs_background $SNPS_BAL $SUMMITS | awk '{print "bal\t" $0}';
 snps_vs_background $SNPS_VRG $SUMMITS | awk '{print "vrg\t" $0}') \
    | ggplot -o SNPs_summits.pdf -W 6 -H 4 \
    "$ggplot_string + ggtitle(\"SNVs in DHS 100bp summits and randomly matched intervals\")" \
    stdin
rm $SUMMITS



echo "Counting SNPs in CAD regions ($NUM_SHUFFLES times)"
CAD=../../tracks/enhancer/CAD4_plus_vienna_dm6.core80percent.onlyChr2and3.bed
(snps_vs_background $SNPS_BAL $CAD | awk '{print "bal\t" $0}';
 snps_vs_background $SNPS_VRG $CAD | awk '{print "vrg\t" $0}') \
    | ggplot -o SNPs_CAD.pdf -W 6 -H 4 \
    "$ggplot_string + ggtitle(\"SNVs in CAD regions and randomly matched intervals\")" \
    stdin



rm $SNPS_BAL
rm $SNPS_VRG