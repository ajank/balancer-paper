REF=/g/korbel/shared/datasets/refgenomes/fly/dm6.fa.fai
our_tads=/g/korbel/shared/projects/drosophila_balancer/analyses/tracks/TADs/TAD_calls/HiC_DB_6-8h_All_dm6_5000_IS100000.bed

get_flanks()
{
    # $1 = bed file
    # $2 = flank radius 
    if [ ! -f genome.fai ]; then
        cut -f-2 $REF | grep -P '^chr[23X][LR]?\t' > genome.fai
    fi
    bedtools flank -i $1 -g genome.fai -b 1 \
        | bedtools slop -i stdin -b $2 -g genome.fai \
        | bedtools sort -i stdin \
        | uniq
}


get_flanks ../../tracks/TADs/external/eagon_kornberg.bed 5000 > boundaries/eagon_kornberg.10kb.bed
get_flanks ../../tracks/TADs/external/sexton_cavalli.bed 5000 > boundaries/sexton_cavalli.10kb.bed
get_flanks ../../tracks/TADs/external/hou_corces.bed     5000 > boundaries/hou_corces.10kb.bed

get_flanks ../../tracks/TADs/external/eagon_kornberg.bed 25000 > boundaries/eagon_kornberg.50kb.bed
get_flanks ../../tracks/TADs/external/sexton_cavalli.bed 25000 > boundaries/sexton_cavalli.50kb.bed
get_flanks ../../tracks/TADs/external/hou_corces.bed     25000 > boundaries/hou_corces.50kb.bed

# our own TAD calls
awk '{OFS="\t"; print $1, ($2+$3)/2, ($2+$3)/2}' $our_tads | bedtools slop -i stdin -g genome.fai -b 5000  > boundaries/ours.10kb.bed
awk '{OFS="\t"; print $1, ($2+$3)/2, ($2+$3)/2}' $our_tads | bedtools slop -i stdin -g genome.fai -b 25000 > boundaries/ours.50kb.bed




test_intersection()
{
    # $1 = bed file of TAD boundaries
    # $2 = bed file of Chip peaks for arc. proteins
    # $3 = number of shuffles
    if [ ! -f genome.fai ]; then
        cut -f-2 $REF | grep -P '^chr[23X][LR]?\t' > genome.fai
    fi
    echo -en "actual\t"; bedtools intersect -a $2 -b $1 -wa | wc -l;

    for i in $(seq 1 $3)
    do
        echo -en "shuffle\t"; 
        bedtools shuffle -i $1 -g genome.fai \
            | bedtools sort -i stdin \
            | bedtools intersect -a $2 -b stdin -wa \
            | wc -l
    done
}

NUM_SIM=500
arcP_file_=( /g/korbel/shared/projects/drosophila_balancer/analyses/tracks/TADs/chipData/ChipSeq.CTCF.embryo_stage2.Boyle_2014.bed
            /g/korbel/shared/projects/drosophila_balancer/analyses/tracks/TADs/chipData/ChipChip.BEAF-32.embryo_0-12h.White_K_21.bed
            /g/korbel/shared/projects/drosophila_balancer/analyses/tracks/TADs/chipData/ChipChip.CP190.embryo_0-12h.White_K_22.bed
            /g/korbel/shared/projects/drosophila_balancer/analyses/tracks/TADs/chipData/ChipChip.SuHw.embryo_0-12h.White_K_27.bed
            /g/korbel/shared/projects/drosophila_balancer/analyses/tracks/TADs/chipData/ChipChip.ZW5.embryo_2-4h.Karpen_G_5265.bed )
arcP_name_=(CTCF Beaf-32 CP190 SuHw ZW5)



# Test for 2 different boundary sizes
for bnd_size in 10 50
do
    resultFile="results.${bnd_size}kb.txt"
    echo -e "TADcalls\tArcProt\ttest\toverlaps" > $resultFile

    # 4 different TAD boundary calls
    for tadcalls in hou_corces sexton_cavalli eagon_kornberg ours
    do

        # 5 different architectural proteins
        for i in $(seq 0 4)
        do
            arcP_file=${arcP_file_[$i]}
            arcP_name=${arcP_name_[$i]}

            echo "$bnd_size kb:   $tadcalls vs. $arcP_name"
            test_intersection boundaries/${tadcalls}.${bnd_size}kb.bed $arcP_file $NUM_SIM \
                | awk '{print "'${tadcalls}'\t'${arcP_name}'\t" $0}' \
                >> $resultFile
        done
    done
    Rscript plot_enrichment.R $resultFile boundary_enrichment.${bnd_size}kb.pdf
done


