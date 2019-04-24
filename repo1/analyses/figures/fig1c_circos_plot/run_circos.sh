#!/bin/bash

if [ "$#" -ne 5 ];
        then
        echo ""
        echo "Usage: $0 DEL DUP TRA Karyotype OUT"
        echo ""
        echo "All arguments point to files"
        echo ""
        exit 1
fi

DEL=$1
DUP=$2
TRA=$3
KAR=$4
OUT=$5

P=$(dirname $OUT)
#out=$(mktemp --suffix=.svg)

# Prepare config files
sed -e 's|DELETIONINPUTFILE|'$DEL'|' \
    -e 's|DUPLICATIONINPUTFILE|'$DUP'|' \
    -e 's|TRANSLOCATIONINPUTFILE|'$TRA'|' \
    -e 's|KARYOTYPEFILE|'$KAR'|' \
    $P/circos.conf \
    > $P/circos.conf.tmp

# run circos
/g/korbel/meiers/tools/circos-0.67-7/bin/circos \
    -noparanoid \
    -conf $P/circos.conf.tmp \
    -outputdir $(dirname $OUT) \
    -outputfile $(basename $OUT) \
    > /dev/null

# Add a label
# Add custom information:
#head -n 3 $out > $OUT
#echo "<text x=\"10\" y=\"30\" font-size=\"26.0px\" font-family=\"CMUBright-Roman\" style=\"text-anchor:start;fill:rgb(115,115,115);\">$1</text>" >> $OUT
#echo "<text x=\"10\" y=\"60\" font-size=\"26.0px\" font-family=\"CMUBright-Roman\" style=\"text-anchor:start;fill:rgb(115,115,115);\">min minDV: $DV</text>" >> $outfile.svg
#echo "<text x=\"10\" y=\"90\" font-size=\"26.0px\" font-family=\"CMUBright-Roman\" style=\"text-anchor:start;fill:rgb(115,115,115);\">min mapQ: $MAPQ</text>" >> $outfile.svg
#tail -n+4 $out >> $OUT

# convert to png
#rsvg-convert -h 1000 $out > $OUT




