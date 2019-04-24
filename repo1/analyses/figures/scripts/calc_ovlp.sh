# $1 = bed file a
# $2 = bed file b
# $3 = -f parameter
# $4 = genome.chr23.fai
flags="-f $3"

echo -en "real\t"
bedtools intersect -a ${1} -b ${2} ${flags} | wc -l
for i in $(seq 1 500)
do
    echo -en "rand\t"
    cat ${1} | bedtools shuffle -g $4 -i - | sort -k1,1 -k2,2n | bedtools intersect -a - -b ${2} ${flags} | wc -l
done
