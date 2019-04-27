#!/bin/bash
set -e
DELLY_all=../DEL/DEL.all.bed
FB_VCF=../../../data/variants/SNVs2/wgs.freebayes-k.filter.norm.vcf.gz
FB_del=small_dels.bed
FINAL=DEL.final.bed
tmp=delly.rest.bed

DUP_bcf=../DUP/delly0.7.5.DUP.final.bcf
DUP_manual=../MANUAL/manual_calls.txt
DUP_FINAL=DUP.final.bed


function smallDels2bed {
	awk -f ../vcfinfo.awk -W source '!/^#/ && $1~/^chr[23X][LR]?$/ \
                                             && length($4) > 1 && length($5)==1 \
	                                         && (SAMPLE[1,"GT"] == "0/1" && SAMPLE[2,"GT"] == "1/1" \
										      || SAMPLE[1,"GT"] == "0/1" && SAMPLE[2,"GT"] == "0/0" \
										      || SAMPLE[1,"GT"] == "1/1" && SAMPLE[2,"GT"] == "1/1") \
										   {OFS="\t"; print $1,$2-1,$2+length($4)-1,substr($10,1,3) "_" substr($11,1,3)}' \
								     | sort -k1,1 -k2,2n
}

# Only keep Delly calls which are not in the FB VCF
echo "Get small deletion (FreeBayes)"
less $FB_VCF | smallDels2bed > $FB_del
echo "Overlap with Delly DEL calls"
bedtools intersect -a $DELLY_all -b $FB_del -loj -f 0.5 -r \
	| awk '$5=="." && $6=="-1"' \
	| cut -f -4 > $tmp
echo "Create DEL.final.bed"
cat $FB_del $tmp | sort -k1,1 -k2,2n | awk '$3-$2>=15' > $FINAL
rm $tmp


echo "DUPLICATIONS"
cat <(bcftools view $DUP_bcf \
         | awk -f ../vcfinfo.awk -W source \
           '!/^#/ {OFS="\t"; print $1,$2-1,INFO["END"], substr($10,1,3) "_" substr($11,1,3) "_" $3}') \
     <(cut -f-4 $DUP_manual) \
     | sort -k1,1 -k2,2n \
     > $DUP_FINAL


## Part 2: statistics

echo "Plot Small DEL size distribution"
awk '{print $3-$2}' $FB_del \
    | ggplot -o small_dels.pdf -W 4 -H 3 \
      "ggplot(d) + aes(V1) + geom_histogram(binwidth=1) + theme_minimal() + xlab('size')"\
      stdin
awk '$3-$2>14 {print $3-$2}' $FB_del \
    | ggplot -o small_dels.above14.pdf -W 4 -H 3 \
      "ggplot(d) + aes(V1) + geom_histogram(binwidth=1) + theme_minimal() + xlab('size')"\
      stdin

echo "Plot Large DEL size distribution"
awk '{print $3-$2}' $DELLY_all \
    | ggplot -o DEL.all.pdf -W 4 -H 3 \
      "ggplot(d) + aes(V1) + geom_histogram(binwidth=10) + theme_minimal() + xlab('size')"\
      stdin
awk '$3-$2 < 100 {print $3-$2}' $DELLY_all \
    | ggplot -o DEL.all.below100.pdf -W 4 -H 3 \
      "ggplot(d) + aes(V1) + geom_histogram(binwidth=1) + theme_minimal() + xlab('size')"\
      stdin

echo "GTs FreeBayes (19-49bp)"
awk '$3-$2 >18 && $3-$2<50 {print substr($4,1,7)}' $FB_del | sort | uniq -c 
echo "GTs Delly (19-49bp)"
awk '$3-$2 >18 && $3-$2<50 {print substr($4,1,7)}' $DELLY_all | sort | uniq -c 

echo "GTs in final call set (<50bp):"
awk '$3-$2<50 {print substr($4,1,7)}' $FINAL | sort | uniq -c
echo "GTs in final call set (50-159bp):"
awk '$3-$2>=50 && $3-$2<160 {print substr($4,1,7)}' $FINAL | sort | uniq -c
echo "GTs in final call set (>=160bp):"
awk '$3-$2>=160 {print substr($4,1,7)}' $FINAL | sort | uniq -c




