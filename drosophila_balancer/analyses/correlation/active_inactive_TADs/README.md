# `joined_TADs.R`

Add explanation later...


# Do TADs of different class join more often than expected?

First we need a background model of where random breakpoints would fall. 
This depends on how much of the genome is flagged as active/inactive/NUll.

Unclassified regions (centromers/telomers):
```
grep -P '^chr[23][LR]' ../../TADs/data/TADs_Sexton_2012_embryo_16-18h_dm6.bed | bedtools complement -g ../../correlation/tad_boundaries/genome.fai -i - | grep -P '^chr[23][LR]' | awk '{sum+=$3-$2} END {print sum}'
```

Sum of active / Inactive regions
```
awk '$1~/^chr[23][LR]/ {sum[$4]+= $3-$2} END {for (x in sum) {print x "\t" sum[x]}}' ../../TADs/data/TADs_Sexton_2012_embryo_16-18h_dm6.b
ed
```

### Expectation:

The code above yields these numbers

type            | size            | frequency
---             | ---:            | ---:
Unclassified    | 15,046,056      | 13.80%
Active          | 19,572,600      | 17.96%
HP1_centromeric | 4,903,768       |  4.50%
PcG             | 12,931,282      | 11.86%
Null            | 56,536,500      | 51.87%
**total**       | **108,990,206** |100.00%

This leads to respective frequencies for a random junction to connect two regions of the same class of `17.96*17.96 + 4.50*4.50 + 11.86*11.86 + 51.87*51.87 = 31.74%`.

The frequency that one or two ends are not classified is `13.80*1 + (100-13.80)*13.80 = 25.7%`.

Consequently, the frequency of ending up in different TADs is `42.65%`.

### A "p-value" for the number of TADs joined 

Note that `n=15` due to 3 unclassified merges.

```
> 1 - pbinom(8, 15, 0.4265)                                                                                                                                            
[1] 0.1365004
```