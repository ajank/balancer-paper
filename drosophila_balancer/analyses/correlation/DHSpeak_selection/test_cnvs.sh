del=../../SV_filter/FINAL/DEL.final.bed
delvrg=DEL.vrg-spec.bed
delbal=DEL.bal-spec.bed
awk 'substr($4,1,7)=="0/1_0/0"' $del > $delbal
awk 'substr($4,1,7)=="0/1_1/1"' $del > $delvrg

# Prepare DHS files
DHS=$(mktemp --suffix ".bed")
grep -P '^chr[23][LR]' ../../tracks/enhancer/DNase_HS_sites_stages9-11_HotSpot_peaks_FDR-1pc_liftedToDm6.bed | sort -k1,1 -k2,2n > $DHS

SUMMITS=$(mktemp --suffix ".bed")
grep -P '^chr[23][LR]' ../../tracks/enhancer/macs_summits_p50_slop50_merge25_clean.htseq.dm6.bed | sort -k1,1 -k2,2n > $SUMMITS

# Prepare exons
EXONS=$(mktemp --suffix ".bed")
grep -P '^chr[23][LR]' ../../tracks/exons.bed | sort -k1,1 -k2,2n > $EXONS


# Prepare plotting function
ggplot_string='
    pval = (frank(d, V2)[d$V1=="real"] + 0) / (nrow(d[V1=="rand",])+1);
    bw = ifelse(max(d$V2)-min(d$V2) > 100, 10, 5)
    ggplot(d[d$V1=="rand",]) + 
        aes(V2) + 
        geom_histogram(binwidth = bw) + 
        geom_vline(data = d[d$V1=="real",], aes(xintercept = V2), color="red") +
        geom_label(x=Inf, y=Inf, hjust=1, vjust=1, label=paste0("hat(p)==",round(pval,3)), inherit.aes=F, parse=T)'

overlap_vs_background() {
    # $1 = bed file a
    # $2 = bed file b
    # $3 = -f parameter
    flags=""
    if [[ $# -eq 3 ]]; then 
        flags="-f $3"
    fi
    echo -en "real\t"
    bedtools intersect -a ${1} -b ${2} ${flags} | wc -l
    for i in $(seq 1 500)
    do
        echo -en "rand\t"
        cat ${1} | bedtools shuffle -g genome.chr23.fai -i - | sort -k1,1 -k2,2n | bedtools intersect -a - -b ${2} ${flags} | wc -l
    done
}

echo "Exons vs. DEL, vrg"
overlap_vs_background $EXONS $delvrg 0.5 | ggplot -o exons_delvrg_min50pc.pdf -W 6 -H 4 \
    "$ggplot_string + ggtitle(\"Exons hit by deletions (50%) in Virginizer\")" \
    stdin
echo "Exons vs. DEL, bal"
overlap_vs_background $EXONS $delbal 0.5 | ggplot -o exons_delbal_min50pc.pdf -W 6 -H 4 \
    "$ggplot_string + ggtitle(\"Exons hit by deletions (50%) in Balancer\")" \
    stdin

echo "DHS peaks vs. DEL, vrg"
overlap_vs_background $DHS $delvrg 0.5 | ggplot -o DHS_delvrg_min50pc.pdf -W 6 -H 4 \
    "$ggplot_string + ggtitle(\"DHS peaks hit by deletions (50%) in Virginizer\")" \
    stdin
overlap_vs_background $DHS $delvrg 0.05 | ggplot -o DHS_delvrg_min5pc.pdf -W 6 -H 4 \
    "$ggplot_string + ggtitle(\"DHS peaks hit by deletions (5%) in Virginizer\")" \
    stdin
overlap_vs_background $DHS $delvrg | ggplot -o DHS_delvrg.pdf -W 6 -H 4 \
    "$ggplot_string + ggtitle(\"DHS peaks hit by deletions (1bp) in Virginizer\")" \
    stdin

echo "DHS peaks vs. DEL, bal"
overlap_vs_background $DHS $delbal 0.5 | ggplot -o DHS_delbal_min50pc.pdf -W 6 -H 4 \
    "$ggplot_string + ggtitle(\"DHS peaks hit by deletions (50%) in Balancer\")" \
    stdin
overlap_vs_background $DHS $delbal 0.05 | ggplot -o DHS_delbal_min5pc.pdf -W 6 -H 4 \
    "$ggplot_string + ggtitle(\"DHS peaks hit by deletions (5%) in Balancer\")" \
    stdin
overlap_vs_background $DHS $delbal | ggplot -o DHS_delbal.pdf -W 6 -H 4 \
    "$ggplot_string + ggtitle(\"DHS peaks hit by deletions (1bp) in Balancer\")" \
    stdin


echo "DHS summits vs. DEL"
overlap_vs_background $SUMMITS $delvrg 0.05 | ggplot -o summits_delvrg_min5pc.pdf -W 6 -H 4 \
    "$ggplot_string + ggtitle(\"DHS summits hit by deletions (5%) in Virginizer\")" \
    stdin
overlap_vs_background $SUMMITS $delbal 0.05 | ggplot -o summits_delbal_min5pc.pdf -W 6 -H 4 \
    "$ggplot_string + ggtitle(\"DHS summits hit by deletions (5%) in Balancer\")" \
    stdin



CAD=../../tracks/enhancer/CAD4_plus_vienna_dm6.core80percent.onlyChr2and3.bed
echo "CAD4 intervals vs. DEL, vrg"
overlap_vs_background $CAD $delvrg 0.5 | ggplot -o CAD_delvrg_min50pc.pdf -W 6 -H 4 \
    "$ggplot_string + ggtitle(\"CAD peaks hit by deletions (50%) in Virginizer\")" \
    stdin
overlap_vs_background $CAD $delvrg 0.05 | ggplot -o CAD_delvrg_min5pc.pdf -W 6 -H 4 \
    "$ggplot_string + ggtitle(\"CAD peaks hit by deletions (5%) in Virginizer\")" \
    stdin
overlap_vs_background $CAD $delvrg | ggplot -o CAD_delvrg.pdf -W 6 -H 4 \
    "$ggplot_string + ggtitle(\"CAD peaks hit by deletions (1bp) in Virginizer\")" \
    stdin

echo "CAD peaks vs. DEL, bal"
overlap_vs_background $CAD $delbal 0.5 | ggplot -o CAD_delbal_min50pc.pdf -W 6 -H 4 \
    "$ggplot_string + ggtitle(\"CAD regions hit by deletions (50%) in Balancer\")" \
    stdin
overlap_vs_background $CAD $delbal 0.05 | ggplot -o CAD_delbal_min5pc.pdf -W 6 -H 4 \
    "$ggplot_string + ggtitle(\"CAD regions hit by deletions (5%) in Balancer\")" \
    stdin
overlap_vs_background $CAD $delbal | ggplot -o CAD_delbal.pdf -W 6 -H 4 \
    "$ggplot_string + ggtitle(\"CAD regions hit by deletions (1bp) in Balancer\")" \
    stdin



DUP=$(mktemp --suffix ".bed")

echo "DHS vs. DUPs, vrg"
awk '$4=="vrg-spec"' ../../tracks/SVs/DUP.bed > $DUP
overlap_vs_background $DHS $DUP 0.05 | ggplot -o DHS_dupvrg_min5pc.pdf -W 6 -H 4 \
    "$ggplot_string + ggtitle(\"DHS peaks hit by duplication (5%) in Virginizer\")" \
    stdin
overlap_vs_background $DHS $DUP 0.5 | ggplot -o DHS_dupvrg_min50pc.pdf -W 6 -H 4 \
    "$ggplot_string + ggtitle(\"DHS peaks hit by duplication (50%) in Virginizer\")" \
    stdin

echo "DHS vs. DUPs, bal"
awk '$4=="bal-spec"' ../../tracks/SVs/DUP.bed > $DUP
overlap_vs_background $DHS $DUP 0.05 | ggplot -o DHS_dupbal_min5pc.pdf -W 6 -H 4 \
    "$ggplot_string + ggtitle(\"DHS peaks hit by duplication (5%) in Balancer\")" \
    stdin
overlap_vs_background $DHS $DUP 0.5 | ggplot -o DHS_dupbal_min50pc.pdf -W 6 -H 4 \
    "$ggplot_string + ggtitle(\"DHS peaks hit by duplication (50%) in Balancer\")" \
    stdin



echo "CAD vs. DUPs, vrg"
awk '$4=="vrg-spec"' ../../tracks/SVs/DUP.bed > $DUP
overlap_vs_background $CAD $DUP 0.05 | ggplot -o CAD_dupvrg_min5pc.pdf -W 6 -H 4 \
    "$ggplot_string + ggtitle(\"CAD regions hit by duplication (5%) in Virginizer\")" \
    stdin
overlap_vs_background $CAD $DUP 0.5 | ggplot -o CAD_dupvrg_min50pc.pdf -W 6 -H 4 \
    "$ggplot_string + ggtitle(\"CAD regions hit by duplication (50%) in Virginizer\")" \
    stdin

echo "CAD vs. DUPs, bal"
awk '$4=="bal-spec"' ../../tracks/SVs/DUP.bed > $DUP
overlap_vs_background $CAD $DUP 0.05 | ggplot -o CAD_dupbal_min5pc.pdf -W 6 -H 4 \
    "$ggplot_string + ggtitle(\"CAD regions hit by duplication (5%) in Balancer\")" \
    stdin
overlap_vs_background $CAD $DUP 0.5 | ggplot -o CAD_dupbal_min50pc.pdf -W 6 -H 4 \
    "$ggplot_string + ggtitle(\"CAD regions hit by duplication (50%) in Balancer\")" \
    stdin

rm $DUP

# clean tmp files
rm $DHS
rm $SUMMITS

