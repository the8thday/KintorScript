#! /usr/bin/env bash

fname=$1


quast.py ./corona/scaffolds.fasta -t 4 -R /mnt/d/Brazil_seq/covid19.fasta \
-g /mnt/d/Brazil_seq/GCF_009858895.2_ASM985889v3_genomic.gff \
-l ${fname} -o quast_corona
quast.py ./spades/scaffolds.fasta -t 4 -R /mnt/d/Brazil_seq/covid19.fasta \
-l ${fname} -o quast_scaffolds

