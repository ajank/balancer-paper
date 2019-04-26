#!/bin/bash

deseq=../ase/deseq/DESeq.N1_6-8h.standardFormat.txt

awk '$7>=0.05           {print $1}'     $deseq > FBgn.control.txt
awk '$7 <0.05           {print $1}'     $deseq > FBgn.all_ase.txt
awk '$7 <0.05 && $3 > 0 {print $1}'     $deseq > FBgn.all_up-regulated.txt
awk '$7 <0.05 && $3 < 0 {print $1}'     $deseq > FBgn.all_down-regulated.txt
awk '$7 <0.05 && ($3<-0.585||$3>0.585)' $deseq > FBgn.ase_with_fc1.5.txt
