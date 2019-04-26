#!/bin/bash

cd analysis/balancer_cap2

rm -f igv_sessions.zip
(
  ls annotations/assembly_dm6bal3.bed
  ls annotations/CAD4_plus_vienna_minus_inactive_corrected_names_dm6.bed
  ls annotations/DEL.final.BAL.bed
  ls annotations/DEL.final.VRG.bed
  ls annotations/DNase_HS_sites_Stage10_merged_norm1X_bin_20bp.bigwig
  ls annotations/DNase_HS_sites_Stage11_merged_norm1X_bin_20bp.bigwig
  ls annotations/DNase_HS_sites_Stage9_merged_norm1X_bin_20bp.bigwig
  ls annotations/DNase_HS_sites_stages9-11_HotSpot_peaks_FDR-1pc_liftedToDm6.bed
  ls annotations/DpnII_dm6_balancer_affected_RS.bed
  ls annotations/DUP.final.BAL.bed
  ls annotations/DUP.final.VRG.bed
  ls annotations/Filion_2010_colors_dm6.bed.gz
  ls annotations/final.predicted_MEI.bed
  ls annotations/Insulator_CRMs.chip_seq.cell_lines.gff
  ls annotations/MesoCRM_dm6_Nature_Zinzen2009.gff
  ls annotations/RNA-seq_N1_6-8h_All.fwd.bigWig
  ls annotations/RNA-seq_N1_6-8h_All.rev.bigWig
  ls annotations/RNA-seq_N1_6-8h_BAL.fwd.bigWig
  ls annotations/RNA-seq_N1_6-8h_BAL.rev.bigWig
  ls annotations/RNA-seq_N1_6-8h_VRG.fwd.bigWig
  ls annotations/RNA-seq_N1_6-8h_VRG.rev.bigWig
  ls annotations/TADs_HiC_DB_6-8h_All_HiCExplorer_dm6.bed.gz
  ls annotations/TADs_HiC_DB_6-8h_BAL_HiCExplorer_dm6.bed.gz
  ls annotations/TADs_HiC_DB_6-8h_VRG_HiCExplorer_dm6.bed.gz
  ls annotations/TADs_Ramirez_2017_Kc167_dm6.bed.gz
  ls annotations/TADs_Sexton_2012_embryo_16-18h_dm6.bed.gz
  ls annotations/viewpoints.gff

  ls contacts_noARS_all/igv_sessions/*_igv_session.xml
  ls contacts_noARS_all/igv_sessions_differential/*_igv_session.xml
  ls contacts_noARS_all/igv_sessions_differential_CAD4_or_MesoCRM/*_igv_session.xml
  ls contacts_noARS_all/igv_sessions_differential_DHS_or_CAD4_or_MesoCRM/*_igv_session.xml
  ls contacts_all/tracks/*_All.RPGC.bedGraph.gz
  ls contacts_all/tracks/*_BAL.RPGC.bedGraph.gz
  ls contacts_all/tracks/*_VRG.RPGC.bedGraph.gz
  # ls contacts_min40bp/tracks/*_All_interactions.gff
  # ls contacts_min40bp/tracks/*_All_interactions_log.p.bedGraph.gz
  # ls contacts_min40bp/tracks/*_All_interactions_score.bedGraph.gz
  ls contacts_noARS_all/tracks/*_differential.gff
  # ls contacts_noARS_all/tracks/*_differential_BAL.gff
  # ls contacts_noARS_all/tracks/*_differential_VRG.gff
) | sort | zip -@ -q igv_sessions

cd ../..
