#! /usr/bin/bash

fq1=$1
fq2=$2
fname=$3

mkdir cleandata jellyfish spades corona quast_corona quast_scaffolds

fastp --thread=4 --detect_adapter_for_pe --cut_front --cut_tail --cut_right  --cut_front_window_size=1 \
--cut_front_mean_quality=3 --cut_tail_window_size=1 --cut_tail_mean_quality=3 --cut_right_window_size=4 \
--cut_right_mean_quality=15 --length_required=36 --trim_front1=10 --trim_front2=10 --compression=4 \
-i ${fq1} \
-o ./cleandata/${fname}_r1.clean.fq.gz \
-I ${fq2} \
-O ./cleandata/${fname}_r2.clean.fq.gz \
--json ${fname}.fastp.json \
--html ${fname}.fastp.html \
--report_title "${fname} QC Reporter"

jellyfish count -C -m 21 -s 1000M -t 4 -o ./jellyfish/${fname}.jf -F 2 <(zcat ./cleandata/${fname}_r1.clean.fq.gz) <(zcat ./cleandata/${fname}_r2.clean.fq.gz)

jellyfish histo -t 8 ./jellyfish/${fname}.jf > ./jellyfish/${fname}.histo

spades.py --pe1-1 ./cleandata/${fname}_r1.clean.fq.gz \
--pe1-2 ./cleandata/${fname}_r2.clean.fq.gz \
--isolate -t 4 -o ./spades

coronaspades.py -1 ./cleandata/${fname}_r1.clean.fq.gz \
-2 ./cleandata/${fname}_r2.clean.fq.gz \
-t 4 -o ./corona

sed -n '1p' ./spades/scaffolds.fasta | sed 's/>//' > ./spades/id.list
seqtk subseq ./spades/scaffolds.fasta ./spades/id.list > ${fname}.fa
