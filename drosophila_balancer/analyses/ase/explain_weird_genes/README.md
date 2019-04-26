# Strong ASE signal probably due to TE insertions

Note: Gene `FBgn0031414` (eys) on the plus strand and gene `FBgn0031414` (CG9967)
on the minus strand are affected by the same TE insertion around 
chr2L:2,355,358-2,355,391.

## Results

These are the results of manually going through the list of high-ASE genes and
checking for MEIs:

| Gene | +/- | ASE lfc (bal/vrg) | Explanation | TE family | MEI pos | Comment |
| ---- | --- | ----------------- | ----------- | --------- | ------- | ------- |
| FBgn0265959 | + | 5.42  | MEI     | gb/AY180917/roo | downstream, + | |
| FBgn0034085 | - | 5.10  | MEI     | gb/AY180917/roo | downstream, - | |
| FBgn0051116 | - | 4.91  | MEI     | gb/AY180917/roo | downstream, - | |
| FBgn0015038 | - | 4.59  | MEI     | gb/X02599/copia | downstream, - | |
| FBgn0031461 | + | 4.36  | UNCLEAR | gb/AJ278684/Rt1a| very far upstream, - | MEI within exon of Syt1(-); poly-A; also low mapp. region close? |
| FBgn0036403 | + | 4.31  | UNCLEAR | gb/AY180917/roo | upstream, +, within 3'UTR of CG6650 (+) | Fusion or MEI?? Only in N1, but not in F1?? | 
| FBgn0038783 | - | 4.16  | UNCLEAR | gb/nnnnnnnn/invader4 or gb/AC005453/1360 | | DEL from 20308159 to ~20326000 might create fusion gene |
| FBgn0047351 | - | 3.76  | MEI     | gb/AY180917/roo | downstream, - | in 5'UTR; likely also expresses the subsequent gene (FBgn0043005) |
| FBgn0031414 | + | 3.74  | MEI     | gb/AY180917/roo | downstream, + | Hits intron of eys (FBgn0031414,+) and intron of CG9967 (FBgn0031413,-) |
| FBgn0052985 | + | 3.63  | MEI     | gb/V00246/FB    | downstream, ? | also strong drop in coverage as if there was an SV |
| FBgn0039380 | + | 3.65  | UNCLEAR | gb/AY180917/roo x2 | Deletion between two inserted MEs in the 3'UTR of the gene ||	
| FBgn0265912 | + | 3.65  | UNCLEAR |                 | | insertion far upstream, didn't match any TE class |
| FBgn0036224 | + | 3.20  | MEI     | gb/AY180917/roo | upstream, + | |
| FBgn0051164 | + | 3.09  | MEI     | gb/AY180917/roo | downstream, + | |
| FBgn0050044 | + | 2.80  | MEI     | gb/AY180917/roo | upstream, + | Additionally hits intron/5'UTR of fdl gene |
| FBgn0033469 | + | 2.66  | MEI     | gb/AY180917/roo | upstream, + | Hits intron of Pal1 gene and causes ASE in last exon of Pal1 |
| FBgn0024150 | - | 2.41  | MEI     | gb/AY180917/roo | downstream, - | Hits an exon |
| FBgn0027552 | - | 2.39  | MEI     | gb/X01472/17.6  | upstream, - | |
| FBgn0043005 | - | 2.34  | MEI     | *gb/AY180917/roo* | upstream, - | same as FBgn0047351 |
| FBgn0038784 | - | 2.19  | MEI     | gb/AC005453/1360 or gb/nnnnnnnn/invader4 | upstream, - | ME might have caused a DEL, need to check |
| FBgn0031219 | - | 2.11  | MEI     | gb/AY180917/roo | upstream, -  | upstream even of the previous gene (FBgn031220) |
| FBgn0083141 | + | 2.01  | MEI     | gb/AY180917/roo | upstream, + | hits 5'UTR of FBgn0264895 (-) |
| FBgn0039104 | - | 2.01  | MEI     | *gb/AY180917/roo* | upstream, - | same as in FBgn0043005 and FBgn0047351 |
| FBgn0031220 | - | 2.00  | MEI     | *gb/AY180917/roo* | upstream, - | same as in FBgn0031219 |
| FBgn0000527 | - | 2.00  | MEI     | gb/nnnnnnnn/412 | downstream, - | |
| FBgn0034766 | + | 1.26  | UNCLEAR | gb/AY180917/roo (gb/nnnnnnnn/412) | 20kb upstream, - (downstream,+) | both cannot explain expression |
| FBgn0031631 | - | 0.95  | MEI     | gb/AY180917/roo | upstream, - | within intron of FBgn0031632(+) |
| FBgn0039613 | + | -1.25 | MEI     | gb/AY180917/roo | upstream, + | likely also relevant for FBgn0027578(+) |
| FBgn0036454 | + | -1.32 | MEI     | gb/X03431/297   | downstream, + | |
| FBgn0033792 | + | -1.39 | MEI     | gb/AY180917/roo | upstream, + | |
| FBgn0265851 | + | -1.44 | MEI     | gb/AY180917/roo | upstream, + | |
| FBgn0267936 | + | -1.82 | MEI     | gb/X02599/copia | upstream, + | |
| FBgn0035187 | - | -1.88 | MEI     | gb/AY180917/roo | upstream, - | |
| FBgn0259219 | + | -1.98 | MEI     | gb/X02599/copia | upstream, + | |
| FBgn0250910 | - | -2.10 | MEI     | gb/V00246/FB    | downstream, ? | |
| FBgn0026592 | + | -2.44 | MEI     | gb/AY180917/roo | upstream, + | |
| FBgn0266782 | + | -2.60 | MEI     | gb/AY180917/roo | upstream, + | | 
| FBgn0011230 | - | -2.78 | MEI     | gb/AY180917/roo | upstream, + | |
| FBgn0035610 | - | -2.84 | MEI     | gb/AY180917/roo | downstream, - | |
| FBgn0267130 | + | -2.81 | MEI     | gb/AY180917/roo | upstream, ? | 8kb upstream, slighltly expressed |
|      *same* |   |       | MEI     | gb/AY180917/roo | downstream, + | within exon |
| FBgn0034467 | - | -3.78 | MEI     | gb/AY180917/roo | upstream, - | |
| FBgn0262020 | + | -3.88 | MEI     | *gb/AY180917/roo* | upstream, + | same as FBgn0266782 |
| FBgn0035638 | - | -6.14 | MEI     | gb/X02599/copia | upstream, - | |

*mentioned before* = Do not count


## Known problems

Since the original BAM file went through rmdup, there are some unpaired reads in them.
Those will be noted by `hts_fetch`, but not removed, thus they have to be removed manually
or the `bam2fastx` step will fail.

## TE source

http://www.fruitfly.org/p_disrupt/TE.html

