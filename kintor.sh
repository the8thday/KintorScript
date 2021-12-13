# GDSC数据库 https://www.cancerrxgene.org/ 是肿瘤细胞药物敏感性和肿瘤治疗基因组数据最大公共资源

# 主要关注的几株新冠
1，P.1 variants在整个巴西人群中的比例, P.1 has also been called 'B.1.1.28.1'/'B.1.1.248'
2，巴西有哪些实验室/公司在做这方面的研究, 病毒突变株的测序
3，南非的B.1.351突变株, 和P.1共同拥有N501Y突变

You can log in at https://gisaid.org
using your user ID: Ericaren@Kintor
and your initial password: Kintor9939


# annswer for tow Questions above
Q1: https://nymag.com/intelligencer/article/what-we-know-about-the-p1-variant-of-the-coronavirus.html  42%
Three-quarters attack rate of SARS-CoV-2 in the Brazilian Amazon during a largely unmitigated epidemic  66%
Resurgence of COVID-19 in Manaus, Brazil, despite high seroprevalence
https://virological.org/t/genomic-characterisation-of-an-emergent-sars-cov-2-lineage-in-manaus-preliminary-findings/586



# 统计fastq inedx
BEGIN { FS = ":"; }
((NR % 4) == 1) { barcodes[$10]++; }
END {
    for (bc in barcodes) {
        print bc": "barcodes[bc]"";
    }
}


# TF database
http://bioinfo.life.hust.edu.cn/AnimalTFDB
https://maayanlab.cloud/Harmonizome/dataset/ENCODE+Transcription+Factor+Targets
http://jaspar.genereg.net/matrix/MA0007.2/
https://www.grnpedia.org/trrust/result.php?gene=AR&species=human&confirm=1
https://bigd.big.ac.cn/databasecommons/
http://bioinfo.life.hust.edu.cn/HumanTFDB#!/seach_result?query=AR
http://cistrome.org/db/#/
http://bioinfo.life.hust.edu.cn/hTFtarget#!/
http://gtrd.biouml.org/#!






grep -wEA1 --no-group-separator 'NODE_1_length_29866_cov_446.601235' scaffolds.fasta
awk 'NR==1{printf $0"\t";next}{printf /^>/ ? "\n"$0"\t" : $0}' multi_fasta.fa | awk -F"\t" 'BEGIN{while((getline k < "IDs.txt")>0)i[k]=1}{gsub("^>","",$0); if(i[$1]){print ">"$1"\n"$2}}'

# sort multi fasta by length
awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}'  input.fasta  |\
    awk -F '\t' '{printf("%d\t%s\n",length($2),$0);}' |\
    sort -k1,1n | cut -f 2- | tr "\t" "\n"

# 以下为操作步骤
/mnt/d/script/run_spades.sh
# link 文件
head ../id.txt | while read line; do mkdir ${line}; find /mnt/d/Brazil_seq/h348_1611/ -name "${line}*.fastq.gz" -exec ln -s {} ${line} \; ; done
# 生成run spades的脚本
find . -name "*fastq.gz" | sed 'N;s/\n/\t/' | awk '{split($1,a,"/"); split($2,b,"/"); print "cd "a[2]";bash /mnt/d/script/run_spades.sh "a[3]" "b[3]" "a[2]";cd ../"}' > run.sh

rm -r quast_scaffolds quast_corona cleandata jellyfish spades corona 210317258935.fastp.*
# 找到corona文件夹下的文件并对第一行改名
find ./*/corona -name "gene_clusters.fasta" -type f | grep -v 'K' | awk '{split($1,a,"/"); print "seqkit sort -l -r --quiet "$1" > ./"a[2]"/corona/gene_clusters_sort.fasta;sed -i ""'\''""1c >"a[2]"'\''"" ./"a[2]"/corona/gene_clusters_sort.fasta"}'
find ./*/corona -name "gene_clusters_sort.fasta" | grep -v 'K' | awk '{split($1,a,"/"); print "cd "a[2]";sed -n ""'\''""1p""'\''"" ./corona/gene_clusters_sort.fasta | sed ""'\''""s/>//""'\''"" > ./corona/id.list;cd ../"}'
find ./*/corona -name "gene_clusters_sort.fasta" | grep -v 'K' | awk '{split($1,a,"/"); print "cd "a[2]"/corona;seqtk subseq gene_clusters_sort.fasta id.list > "a[2]".fna;cd ../../"}'


# 提取fasta的第一个fasta序列
seqtk subseq scaffolds.fasta id.list
seqkit sort -l -r --quiet scaffolds.fasta
# 提取spades第一个fasta, 先提取第一行的名字，似乎结果中将最长的放在了第一位
find ./*/spades -name "scaffolds.fasta" | grep -v 'K' | awk '{split($1,a,"/"); print "cd "a[2]";sed -n ""'\''""1p""'\''"" ./spades/scaffolds.fasta | sed ""'\''""s/>//""'\''"" > ./spades/id.list;cd ../"}'
find ./*/spades -name "scaffolds.fasta" | grep -v 'K' | awk '{split($1,a,"/"); print "cd "a[2]"/spades;seqtk subseq scaffolds.fasta id.list > "a[2]".fna;cd ../../"}'
find ./*/spades -name "*.fna" -type f | awk '{split($1,a,"/"); print "sed -i ""'\''""1c >"a[2]"'\''"" "$1}'


# run rna seq
bash rna_pipe.sh /data/home/liuc/work/rna_test/5_combined_R1.fastq.gz /data/home/liuc/work/rna_test/5_combined_R2.fastq.gz /data/home/liuc/work/rna_test/ U5
nohup bash rna_pipe.sh /data/home/liuc/work/rna_test/U6/6_combined_R1.fastq.gz \
/data/home/liuc/work/rna_test/U6/6_combined_R2.fastq.gz /data/home/liuc/work/rna_test/U6 U6 &

# use refgenie
export REFGENIE='genome_config.yaml'
refgenie init -c $REFGENIE
refgenie listr
refgenie pull hg38/salmon_sa_index
refgenie pull hg38/salmon_partial_sa_index


nohup salmon quant -l MU -i ~/Software/alias/hg38/salmon_partial_sa_index/default \
--validateMappings \
-p 12 \
-1 ./cleandata/U5_r1.clean.fq.gz -2 ./cleandata/U5_r2.clean.fq.gz -o ./salmon/U5 &

nohup salmon quant --gcBias -l MU -p 12 \
--validateMappings \
-i /data/home/liuc/Reference/RNA/salmon/salmon_index_all/salmon_index_all -g /data/home/liuc/Reference/RNA/salmon/GRCh38.primary_assembly.genome.fa \
-1 ./cleandata/U6_r1.clean.fq.gz -2 ./cleandata/U6_r2.clean.fq.gz -o ./salmon_all &


featureCounts -T 8 -p \
-t exon -g gene_id \
-a  -o ./featureCounts_all.txt ./*/bam/*.bam \
1>all.featureCounts.sh.out 2>all.featureCounts.sh.err




