#!/bin/bash
mkdir mapped_to_dm6 2> /dev/null
for file in wgs/*.bam
    do
    out=${file#wgs/}
    out=${out%.bam}
    echo "$file --> $out"
    samtools view $file | cut -f-4 | sort --parallel=2 -k1,1 | gzip > mapped_to_dm6/${out}.txt.gz
done
