#!/bin/bash
path=/g/korbel/shared/projects/drosophila_balancer/data/raw/wgs
threads=6
bwa mem -t $threads TE/te.fa $path/000000000-AA7FA_WGS_CROSS_15s009796-5-1_Ghavi-helm_lane115s009796_* \
    | samtools view -F 4 -bT TE/te.fa - > mapped_to_TE/AA7FA_WGS_CROSS.bam
bwa mem -t $threads TE/te.fa $path/000000000-AA703_WGS_CROSS_15s009796-2-1_Ghavi-helm_lane115s009796_* \
    | samtools view -F 4 -bT TE/te.fa - > mapped_to_TE/AA703_WGS_CROSS.bam
bwa mem -t $threads TE/te.fa $path/000000000-AA7FA_WGS_CROSS_15s009796-5-1_Ghavi-helm_lane115s009796_* \
    | samtools view -F 4 -bT TE/te.fa - > mapped_to_TE/AA7FA_WGS_CROSS.bam
bwa mem -t $threads TE/te.fa $path/000000000-AABU7_WGS_VRG_15s009795-1-1_Ghavi-helm_lane115s009795_* \
    | samtools view -F 4 -bT TE/te.fa - > mapped_to_TE/AABU7_WGS_VRG.bam
bwa mem -t $threads TE/te.fa $path/000000000-AABWP_WGS_VRG_15s009795-2-1_Ghavi-helm_lane115s009795_* \
    | samtools view -F 4 -bT TE/te.fa - > mapped_to_TE/AABWP_WGS_VRG.bam
bwa mem -t $threads TE/te.fa $path/000000000-AACKU_VRG_WGS_15s009530-1-1_Ghavi-helm_lane115s009530_* \
    | samtools view -F 4 -bT TE/te.fa - > mapped_to_TE/AACKU_VRG_WGS.bam
bwa mem -t $threads TE/te.fa $path/000000000-AACLV_WGS_CROSS_15s009796-1-1_Ghavi-helm_lane115s009796_* \
    | samtools view -F 4 -bT TE/te.fa - > mapped_to_TE/AACLV_WGS_CROSS.bam
bwa mem -t $threads TE/te.fa $path/000000000-AAD51_DB_CROSS_WGS_15s009531-1-1_Ghavi-helm_lane115s009531_* \
    | samtools view -F 4 -bT TE/te.fa - > mapped_to_TE/AAD51_DB_CROSS.bam
bwa mem -t $threads TE/te.fa $path/000000000-AAD9C_WGS_CROSS_15s009796-4-1_Ghavi-helm_lane115s009796_* \
    | samtools view -F 4 -bT TE/te.fa - > mapped_to_TE/AAD9C_WGS_CROSS.bam
bwa mem -t $threads TE/te.fa $path/000000000-AAD9M_WGS_CROSS_15s009796-3-1_Ghavi-helm_lane115s009796_* \
    | samtools view -F 4 -bT TE/te.fa - > mapped_to_TE/AAD9M_WGS_CROSS.bam

for x in mapped_to_TE/*.bam
do 
    echo "$x --> ${x%.bam}.txt"
    samtools view $x | cut -f1 | uniq > ${x%.bam}.txt; 
done