# 系统发育树

# 常用比对工具
FastTree, Clustal, Prank, Mafft, MUSCLE
# 常用进化树工具
Phylip, MEGA, Phyml, RaxML, Iqtree, PAUP

# 整合性工具
PhyloSuite, Newick Utilities
# 直系同源基因
orthofinder

# 系统发育树的一般过程
计算推断系统发育树的主要任务是求出最优树的拓扑结构和估计分支长度；评估的目的是对已经得出的系统发育树的置信度进行评估，常用的方法有自举检验法
# 多重序列比对
mafft --thread 16 --auto input.fasta > output.fasta
mafft --maxiterate 1000 --localpair in > out
# 获得精确的多序列比对后，要过滤掉一些低质量以及高变异度的序列区域, trimAI/Gblocks
# 序列比对完后，用于建树的序列位点必须保证具有良好的同源性。所以需要删除序列分歧很大的区域和gap区域。
trimal -in example1 -out output1 -htmlout output1.html -automated1
Gblocks 80-AcoeOut_aln.fasta -t=d -b5=n
# IQ-tree软件构建系统发育树
iqtree -s example.phy -m MFP -b 100 -T AUTO
#muscle 生成phy格式
#muscle -in mouse_J.pro  -physout seqs.phy
# 评估构建好的进化树

# 可视化网站、工具
iTOL, FigTree, ggtree

# 数据库
https://www.biostars.org/p/427840/
https://www.covid19dataportal.org/
https://www.ncbi.nlm.nih.gov/sars-cov-2/
https://bigd.big.ac.cn/ncov/
https://www.gisaid.org
https://covid-19.ensembl.org/index.html
https://virological.org/
https://nextstrain.org/sars-cov-2/

http://bis.zju.edu.cn/overcovid/ # 整合covid工具的网站


# ------------------------ shell pipeline --------------------------

#!/usr/bin/bash
# 单样本运行

set -e
set -u
set -o pipefail

infa=$1
fname=$2


mafft --thread 16 --auto ${infa} > ${fname}_aln.fasta
trimal -in ${fname}_aln.fasta -out ${fname}.trim.fasta -htmlout ${fname}_trimAl.html -gt 0.9 -cons 60
# 最大似然法
iqtree -s ${fname}.trim.fasta -B 1000 --bnni -T AUTO  --threads-max 16 \
--seqtype DNA \
--mem 16G \
-m MFP --prefix ${fname}
# 系统发育树
# NJ
# Bayesian
# build phylogenic tree
iTOL


# ------------------------ covid19 WGS --------------------------
https://www.genomedetective.com/app/typingtool/virus/
https://pangolin.cog-uk.io/
http://giorgilab.unibo.it/coronannotator/

fq1=$1
fq2=$2
outdir=$3
fname=$4

mkdir cleandata

# QC, 是否需要去除宿主基因组
fastp --thread=4 --detect_adapter_for_pe --cut_front --cut_tail --cut_right  --cut_front_window_size=1 \
--cut_front_mean_quality=3 --cut_tail_window_size=1 --cut_tail_mean_quality=3 --cut_right_window_size=4 \
--cut_right_mean_quality=15 --length_required=36 --trim_front1=10 --trim_front2=10 --compression=4 \
-i 210317258935-2-1-1_S2_L001_R1_001_FQ39428.fastq.gz \
-o ./cleandata/210317258935_r1.clean.fq.gz \
-I 210317258935-2-1-1_S2_L001_R2_001_FQ39428.fastq.gz \
-O ./cleandata/210317258935_r2.clean.fq.gz \
--json 210317258935.fastp.json \
--html 210317258935.fastp.html --report_title "210317258935 QC Reporter"

# KmerGenie可以进行k-mer分析及基因组大小评估；jellyfish/KMC && GenomeScope评估基因组
# kmergenie fastq_list.txt -k 121 -l 15 -t 8 -o ./kmergenie/CRR197319
# jellyfish count -C -m 21 -s 20G -t 8 *.fastq -o reads.jf
# kmc -k21 -t10 -m32 -ci1 -cs10000 @FILES reads tmp/
# kmc_tools transform reads histogram reads.histo -cx10000
jellyfish count -C -m 21 -s 1000000000 -t 4 -o ./jellyfish/reads.jf <(zcat ./cleandata/210317258935_r1.clean.fq.gz) <(zcat ./cleandata/210317258935_r2.clean.fq.gz)
jellyfish stats reads.jf
jellyfish histo -t 8 reads.jf > reads.histo
Rscript /path/genomescope2.0/genomescope.R -i reads.histo -o output_p3 -k 21 -p 3 -n NAME_PREFIX
jellyfish dump -c -t -o mer_counts.txt reads.jf

# Spades 组装 or megahit
# 对于病毒分离株推荐使用--isolate参数，不过从体液中提取的数据似乎应该不算
nohup spades.py -1 ./cleandata/210317258935_r1.clean.fq.gz \
-2 ./cleandata/210317258935_r2.clean.fq.gz \
--careful \
# --isolate \
-t 4 -o ./spades &

coronaspades.py -1 ./cleandata/210317258935_r1.clean.fq.gz \
-2 ./cleandata/210317258935_r2.clean.fq.gz \
-t 4 -o ./corona
# megahit -t 2 -m 0.95 --min-contig-len 500 -1 a1_.fq -2  a2_.fq -o  a_Contig.out
# megahit -t 2 -m 0.95 --min-contig-len 500 -1 b1_.fq -2  b2_.fq -o  b_Contig.out

# 组装结果的评估，QUAST
quast.py ./corona/scaffolds.fasta -t 4 -R /mnt/d/Brazil_seq/covid19.fasta \
-g /mnt/d/Brazil_seq/GCF_009858895.2_ASM985889v3_genomic.gff \
-l 210317258935 -o quast_corona
quast.py ./spades/scaffolds.fasta -t 4 -R /mnt/d/Brazil_seq/covid19.fasta \
-l 210317258935 -o quast_scaffolds

# covid19 genotype
pangolin

# align工具选择
blast
bwa mem

# variant caller
lofreq

# annotation of mutation
https://www.biorxiv.org/content/10.1101/2020.05.31.124966v3.full #coronapp http://giorgilab.unibo.it/coronannotator/
http://corgat.cloud.ba.infn.it/galaxy # https://corgat.readthedocs.io/en/latest/
https://github.com/rcs333/VAPiD
https://virological.org/t/global-framework-for-sars-cov-2-data-analysis-application-to-intrahost-variation-part-1/623  SnpEff specifically its 4.5covid19 version
annovar 


# ------------------------ covid19 多序列突变分析 from fasta --------------------------
# 首先从gisaid下载所有的complete序列，或可进一步对这些序列进行筛选

# multiple sequence alignment (MSA) is carried out by using Clustal Omega
# 对比对后的序列和reference通过mummer进行比对, mummer需要多重比对吗？
nucmer --maxgap=90 --mincluster=65 --forward --prefix=ref_qry GCF_009858895.2_ASM985889v3_genomic.fna gisaid_hcov-19_2021_03_04_amazon.fasta
delta-filter -m ref_qry.delta > filter.delta
show-coords -r -c -l filter.delta > prefix.coords
show-aligns filter.delta NC_045512.2 MN908947.3 | less # 查看具体一条的信息
show-snps -r -T -l filter.delta > filter.snps
delta2vcf < filter.delta > file.vcf # 结果没有含有样本信息，amazing
# show-tilling
dnadiff -d ref_qry.delta -p prefix #整合的工具

# 单个fasta比对
nucmer --maxgap=90 --mincluster=65 --prefix=ref_qry GCF_009858895.2_ASM985889v3_genomic.fna MN908947.3.fasta
delta-filter -1 ref_qry.delta | delta2vcf > filter.vcf
show-snps -Clr ref_qry.delta
show-snps -r -T -l filter.delta > MN908947.3.snps
annovar 


# annotation of mutation
# 其的输入应为：
perl /cluster/home/liucong/Software/CorGAT-master/align.pl --clean F \
--refile GCF_009858895.2_ASM985889v3_genomic.fna \
--multi ~/work/phylogenetic/covid19_china_brazil_sf.trim.fasta --out ALIGN_out.tsv
perl /cluster/home/liucong/Software/CorGAT-master/annotate.pl --in ALIGN_out.tsv
# 如果需要对应单独样本
perl align.pl --clean F --refile GCF_009858895.2_ASM985889v3_genomic.fna \
--multi ../covid19/gisaid_hcov-19_2021_03_04_amazon.fasta --out ALIGN_out.tsv

mv align.tmp single_input
for i in `ls ./single_input/*.fasta`
do
fname=`basename ${i} | cut -f 1 -d "."`
echo ${fname}
perl align.pl --multi ${i} --out ./result/${fname}.ALIGN_out.tsv

perl annotate.pl --in ./result/${fname}.ALIGN_out.tsv --out ./result/${fname}.CorGAT_out.tsv
done


