rscript=/g/korbel/shared/projects/drosophila_balancer/analyses/common/R/overlap.R


if [ "$#" -ne "3" ]
then
    echo "Usage:   overlap.sh query.bed subject.bed outfile.pdf"
    exit 1
fi


# 1 = enhancer bed
# 2 = deletion bed
# 3 = out.pdf
tmp=$(mktemp -d)
for x in $(seq 0.05 0.05 1)
	do echo -en "$x\t"; 
	bedtools intersect -a $1 -b $2 -wb -f $x | wc -l
done | Rscript $rscript $3
rm -r $tmp

