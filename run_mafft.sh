#!/usr/bin/env bash

mafft --thread 2 --auto cov100_spades.fna > cov100_spades_aln.fasta

# trimal -in cov100_spades_aln.fasta -out cov100_spades_aln.trim.fasta -htmlout cov100_spades_aln_trimAl.html -gt 0.8 -cons 60
trimal -in cov100_spades_aln.fasta -out cov100_spades_aln.trim.fasta -automated1 -htmlout cov100_trimal.html

iqtree -s cov100_spades_aln.trim.fasta -B 1000 --bnni -T AUTO -ntmax 4 -m MFP --prefix cov100_spades



# call mutation
perl align.pl --clean F \
--refile GCF_009858895.2_ASM985889v3_genomic.fna \
--multi /mnt/d/Brazil_seq/h348_1611/CovMutation/cov100/cov100.fasta \
--out ALIGN_out.tsv

perl /mnt/d/script/CorGAT-master/annotate.pl --in ALIGN_out.tsv

