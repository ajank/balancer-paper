#!/bin/bash

#BSUB -J jemulti[1-3]
#BSUB -M 2000

_input=(dummy
        C840RACXX_Heads_1_15s015373-1-1_Ghavi-helm_lane615s015373
		C88LWACXX_Heads_2_15s015592-1-1_Ghavi-helm_lane415s015592
		C88LWACXX_Heads_3_15s015593-1-1_Ghavi-helm_lane515s015593 )
_barcode=(dummy
           Heads_1.txt
           Heads_2.txt
		   Heads_3.txt )

f1=../../raw/rna2/${_input[$LSB_JOBINDEX]}_1_sequence.txt.gz
f2=../../raw/rna2/${_input[$LSB_JOBINDEX]}_2_sequence.txt.gz

java -Xmx4g -jar /g/funcgen/gbcs-tools/jemultiplexer/jemultiplexer.jar \
    F1=$f1 F2=$f2 BF=../../raw/rna2/${_barcode[$LSB_JOBINDEX]} \
	UF1=${_barcode[$LSB_JOBINDEX]}.unassigned.1.fq.gz \
	UF2=${_barcode[$LSB_JOBINDEX]}.unassigned.2.fq.gz

