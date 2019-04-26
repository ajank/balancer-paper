#!/bin/bash

if [ $# -lt 3 ]; then
  echo 1>&2 "Usage:  $0 infile.vcf.gz infile.fa outfile.tab"
  exit 1
fi

invcf="$1"
infasta="$2"
outfile="$3"

tmpdir=`mktemp -d`
mkfifo $tmpdir/i1 $tmpdir/i2 $tmpdir/o1 $tmpdir/o2

zcat $invcf | grep -v "^#" | tee $tmpdir/i1 > $tmpdir/i2 &

awk ' BEGIN { FS="\t"; OFS="\t" } { if (length($4) != length($5)) { print "sequence length mismatch in line " NR; exit 1 } print $1, $2-1, $2-1+length($4) }' $tmpdir/i1 |
  bedtools getfasta -fi $infasta -bed - -bedOut > $tmpdir/o1 &

cut -f 4,5 $tmpdir/i2 > $tmpdir/o2 &

(echo -e "chrom\tstart\tend\tREF\tVRG\tBAL"; paste $tmpdir/o1 $tmpdir/o2) > $outfile

wait
rm $tmpdir/i2 $tmpdir/o1 $tmpdir/o2
rm -rf $tmpdir
