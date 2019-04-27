# Read separator

These scripts go through a **name-sorted** bam file mapped to the ref genome
and report whether it overlaps a het SNP. Then the read is classified by the
set of SNPs to either ref, alt, errorneous or ambigous.

## Snakemake

The `snakemake` pipeline takes all BAM files from

```
/g/korbel/shared/projects/drosophila_balancer/data/mapping/rna/bam
```

and separates them based on the SNPs in the `getVCF` directory.
Output files are written into the folder `rna`.

Command: `snakemake --cluster "bsub -M2000" -j 18`

#### Snakefile.ATAC

There is a 2nd pipeline that separates ATAC-seq data into the folder `atac`.

### Folder `counts`

Here I simply count ref/alt RNA reads using `htseq-count` to generate plots of allelic imbalance 
on the gene level.

### Read counts from May 2016

sample                  | allele | count
---                     | ---    | ---
F0_VRG_4-8h_rep1        | alt    | 42042
F0_VRG_4-8h_rep1        | amb    | 45924712
F0_VRG_4-8h_rep1        | err    | 7037
F0_VRG_4-8h_rep1        | ref    | 8454971
F0_VRG_4-8h_rep2        | alt    | 54947
F0_VRG_4-8h_rep2        | amb    | 47935195
F0_VRG_4-8h_rep2        | err    | 9096
F0_VRG_4-8h_rep2        | ref    | 8765086
F0_VRG_heads_rep1       | alt    | 43114
F0_VRG_heads_rep1       | amb    | 48697239
F0_VRG_heads_rep1       | err    | 10305
F0_VRG_heads_rep1       | ref    | 9130952
F0_VRG_heads_rep2       | alt    | 92030
F0_VRG_heads_rep2       | amb    | 39760469
F0_VRG_heads_rep2       | err    | 10386
F0_VRG_heads_rep2       | ref    | 7308051
F1_CyOTM3_heads_rep1    | alt    | 5469325
F1_CyOTM3_heads_rep1    | amb    | 60187768
F1_CyOTM3_heads_rep1    | err    | 20049
F1_CyOTM3_heads_rep1    | ref    | 6103804
F1_CyOTM3_heads_rep2    | alt    | 4314262
F1_CyOTM3_heads_rep2    | amb    | 49530475
F1_CyOTM3_heads_rep2    | err    | 12153
F1_CyOTM3_heads_rep2    | ref    | 4726294
N1_CyO_heads_rep1       | alt    | 2170933
N1_CyO_heads_rep1       | amb    | 44380315
N1_CyO_heads_rep1       | err    | 10104
N1_CyO_heads_rep1       | ref    | 7155282
N1_CyO_heads_rep2       | alt    | 2191268
N1_CyO_heads_rep2       | amb    | 45213380
N1_CyO_heads_rep2       | err    | 12139
N1_CyO_heads_rep2       | ref    | 7218683
N1_TM3_heads_rep1       | alt    | 2381502
N1_TM3_heads_rep1       | amb    | 49741243
N1_TM3_heads_rep1       | err    | 17068
N1_TM3_heads_rep1       | ref    | 8060927
N1_TM3_heads_rep2       | alt    | 2036091
N1_TM3_heads_rep2       | amb    | 42861982
N1_TM3_heads_rep2       | err    | 13021
N1_TM3_heads_rep2       | ref    | 6858172
N1_pool_4-8h_rep1       | alt    | 2206633
N1_pool_4-8h_rep1       | amb    | 84380836
N1_pool_4-8h_rep1       | err    | 18583
N1_pool_4-8h_rep1       | ref    | 11392876
N1_pool_4-8h_rep2       | alt    | 2229679
N1_pool_4-8h_rep2       | amb    | 86213770
N1_pool_4-8h_rep2       | err    | 15469
N1_pool_4-8h_rep2       | ref    | 11382520
N1_pool_6-8h_rep1       | alt    | 4633633
N1_pool_6-8h_rep1       | amb    | 81661815
N1_pool_6-8h_rep1       | err    | 28772
N1_pool_6-8h_rep1       | ref    | 16796516
N1_pool_6-8h_rep2       | alt    | 4441527
N1_pool_6-8h_rep2       | amb    | 84349447
N1_pool_6-8h_rep2       | err    | 29883
N1_pool_6-8h_rep2       | ref    | 16415263
N1_pool_heads_rep1      | alt    | 2716284
N1_pool_heads_rep1      | amb    | 52739557
N1_pool_heads_rep1      | err    | 15400
N1_pool_heads_rep1      | ref    | 8541445
N1_pool_heads_rep2      | alt    | 2262824
N1_pool_heads_rep2      | amb    | 43347018
N1_pool_heads_rep2      | err    | 12119
N1_pool_heads_rep2      | ref    | 6854975
N1sex_pool_6-8h_rep1    | alt    | 6353413
N1sex_pool_6-8h_rep1    | amb    | 86101697
N1sex_pool_6-8h_rep1    | err    | 33870
N1sex_pool_6-8h_rep1    | ref    | 16553408
N1sex_pool_6-8h_rep2    | alt    | 5458129
N1sex_pool_6-8h_rep2    | amb    | 76655742
N1sex_pool_6-8h_rep2    | err    | 32557
N1sex_pool_6-8h_rep2    | ref    | 14211708
