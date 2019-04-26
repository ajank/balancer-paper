#!/bin/bash

dir="analysis/balancer_cap2/annotations"
mkdir -p $dir

#
#  Create temporary directory
#

tmpdir=`mktemp -d -p /scratch/jankowsk`

#
#  CAD4
#

cp -p /g/furlong/project/CAD/CAD4_corrected_names/CAD4_plus_vienna_minus_inactive_corrected_names_dm6.bed $dir

#
#  Insulator CRMs
#

indir="/g/furlong/project/59_ChIP-CRM-insulators/analysis/insulator_cluster_sites"
for f in `cd $indir; ls *.gff`; do
  echo $f
  sed -e "s/ID=/ID.orig=/" $indir/$f > $dir/$f
done

#
#  Zinzen et al. 2009 mesoderm CRMs
#

liftOver -gff /g/furlong/DataWarehouse/Data/datawarehouse/furlong_data/dmelanogaster/features/crms/Meso_CRM/dm3/MesoCRM_dm3_Nature_Zinzen2009.gff /g/furlong/genome/chain_files/dm3ToDm6.over.chain.gz $dir/MesoCRM_dm6_Nature_Zinzen2009.gff /dev/null

#
#  Thomas et al. 2011 DNase-seq data
#

cp -p /g/furlong/DataWarehouse/Data/datawarehouse/external_data/Thomas_et_al/dmelanogaster/sequencing/DNase-seq/whole_embryo/*/signal/dm6/*.bigwig $dir
cp -p /g/furlong/DataWarehouse/Data/datawarehouse/external_data/Thomas_et_al/dmelanogaster/sequencing/DNase-seq/whole_embryo/*/peaks/dm6/*.bed $dir

cat \
  $dir/DNase_HS_sites_stage9_HotSpot_peaks_FDR-1pc_liftedToDm6.bed \
  $dir/DNase_HS_sites_stage10_HotSpot_peaks_FDR-1pc_liftedToDm6.bed \
  $dir/DNase_HS_sites_stage11_HotSpot_peaks_FDR-1pc_liftedToDm6.bed \
  | bedtools sort | bedtools merge > $dir/DNase_HS_sites_stages9-11_HotSpot_peaks_FDR-1pc_liftedToDm6.bed

#
#  TADs and chromatin colors
#

cp -p ../33_Hi-C/data/hic_annotations/dm6/TADs_HiC_DB_6-8h_All_IS100k_dm6.bed.gz $dir
cp -p ../33_Hi-C/data/hic_annotations/dm6/TADs_HiC_DB_6-8h_All_HiCExplorer_dm6.bed.gz $dir
cp -p ../33_Hi-C/data/hic_annotations/dm6/TADs_HiC_DB_6-8h_VRG_HiCExplorer_dm6.bed.gz $dir
cp -p ../33_Hi-C/data/hic_annotations/dm6/TADs_HiC_DB_6-8h_BAL_HiCExplorer_dm6.bed.gz $dir
cp -p ../33_Hi-C/data/hic_annotations/dm6/TADs_Ramirez_2017_Kc167_dm6.bed.gz $dir
cp -p ../33_Hi-C/data/hic_annotations/dm6/TADs_Sexton_2012_embryo_16-18h_dm6.bed.gz $dir
cp -p ../33_Hi-C/data/hic_annotations/dm6/Filion_2010_colors_dm6.bed.gz $dir

#
#  Balancer assembly and viewpoint annotations
#

cp -p ../39_Balancer_Hi-C/analysis/assembly_dm6bal3.bed $dir
cp -p analysis/balancer_cap2/viewpoints.gff $dir

#
#  Balancer RNA-seq data
#

rnadir=/g/korbel/shared/projects/drosophila_balancer/analyses/ase/tracks/rna
chromsizes=/g/furlong/genome/D.melanogaster/Dm6/chrsizes/dm6.ucsc.chrom.sizes

for allele in All VRG BAL; do
  case "$allele" in
    "All") oa=ALL ;;
    "VRG") oa=REF ;;
    "BAL") oa=ALT ;;
  esac

  for strand in fwd rev; do
    for f in `cd $rnadir; ls N1*_pool_6-8h_rep*.$oa.$strand.bedGraph.gz`; do
      echo $f
      /g/furlong/jankowsk/ucsc/bedSort $rnadir/$f $tmpdir/$f.bedGraph
      /g/furlong/jankowsk/ucsc/bedGraphToBigWig $tmpdir/$f.bedGraph $chromsizes $tmpdir/$f.bigWig
    done

    /g/furlong/jankowsk/ucsc/bigWigMerge $tmpdir/*.bigWig $tmpdir/merged.bedGraph
    /g/furlong/jankowsk/ucsc/bedGraphToBigWig $tmpdir/merged.bedGraph $chromsizes $dir/RNA-seq_N1_6-8h_$allele.$strand.bigWig
    rm -f $tmpdir/*
  done
done

#
#  Sascha's SNVs (also MNVs and small indels)
#

rm -f $dir/wgs.freebayes-k.filter.norm.vcf.gz
cp -p /g/korbel/shared/projects/drosophila_balancer/data/variants/SNVs2/wgs.freebayes-k.filter.norm.vcf.gz $dir
tabix -f -p vcf $dir/wgs.freebayes-k.filter.norm.vcf.gz

(
  zcat $dir/wgs.freebayes-k.filter.norm.vcf.gz | grep '^#'
  zcat $dir/wgs.freebayes-k.filter.norm.vcf.gz | awk -F "\t" '{ if ($10~/0\/1/ && $11~/1\/1/) print }'
) | bgzip > $dir/wgs.freebayes-k.filter.norm.VRG.vcf.gz
tabix -f -p vcf $dir/wgs.freebayes-k.filter.norm.VRG.vcf.gz

(
  zcat $dir/wgs.freebayes-k.filter.norm.vcf.gz | grep '^#';
  zcat $dir/wgs.freebayes-k.filter.norm.vcf.gz | awk -F "\t" '{ if ($10~/0\/1/ && $11~/0\/0/) print }'
) | bgzip > $dir/wgs.freebayes-k.filter.norm.BAL.vcf.gz
tabix -f -p vcf $dir/wgs.freebayes-k.filter.norm.BAL.vcf.gz

(
  zcat $dir/wgs.freebayes-k.filter.norm.vcf.gz | grep '^#'
  zcat $dir/wgs.freebayes-k.filter.norm.vcf.gz | awk -F "\t" '{ if ($10~/1\/1/ && $11~/1\/1/) print }'
) | bgzip > $dir/wgs.freebayes-k.filter.norm.common.vcf.gz
tabix -f -p vcf $dir/wgs.freebayes-k.filter.norm.common.vcf.gz

#
#  Sascha's deletion and duplication tracks
#

delfile="/g/korbel/shared/projects/drosophila_balancer/analyses/SV_filter/FINAL/DEL.final.bed"
grep "0/1_1/1" $delfile | cut -f 1-3 > $dir/DEL.final.VRG.bed
grep "0/1_0/0" $delfile | cut -f 1-3 > $dir/DEL.final.BAL.bed
grep "1/1_1/1" $delfile | cut -f 1-3 > $dir/DEL.final.common.bed

dupfile="/g/korbel/shared/projects/drosophila_balancer/analyses/SV_filter/FINAL/DUP.final.bed"
grep "0/1_0/1" $dupfile | cut -f 1-3 > $dir/DUP.final.VRG.bed
grep "0/1_0/0" $dupfile | cut -f 1-3 > $dir/DUP.final.BAL.bed

#
#  Sascha's differential MEI calls
#

meifile="/g/korbel/shared/projects/drosophila_balancer/analyses/MEI_detection/own2/final.predicted_MEI.bed"
awk 'BEGIN { FS="\t"; OFS="\t"; print "track name=\"Differential MEI calls\" description=\"\" visibility=1 itemRgb=on" } \
  { print $1,$2,$3,".",1000,".",$2,$3,($4 == "balancer")? "55,126,184" : "77,175,74" }' $meifile > $dir/final.predicted_MEI.bed

#
#  Final cleanup
#

rmdir $tmpdir
