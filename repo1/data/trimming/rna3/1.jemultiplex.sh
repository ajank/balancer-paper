#!/bin/bash

#BSUB -J jemulti[1-1]
#BSUB -M 2000

_input=(dummy HTFMVBGXX_RNA_CROSS_NEW_16s002675-1-1_Ghavi-helm_lane116s002675 )
_barcode=(dummy CROSS_NEW.txt)

f1=../../raw/rna3/${_input[$LSB_JOBINDEX]}_1_sequence.txt.gz
f2=../../raw/rna3/${_input[$LSB_JOBINDEX]}_2_sequence.txt.gz

java -Xmx4g -jar /g/funcgen/gbcs-tools/jemultiplexer/jemultiplexer.jar \
    F1=$f1 F2=$f2 BF=../../raw/rna3/${_barcode[$LSB_JOBINDEX]} \
	UF1=${_barcode[$LSB_JOBINDEX]}.unassigned.1.fq.gz \
	UF2=${_barcode[$LSB_JOBINDEX]}.unassigned.2.fq.gz

