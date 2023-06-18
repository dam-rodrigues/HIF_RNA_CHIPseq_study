#connecting to the server
ssh -p 12034 acmo02@cluster.di.fct.unl.pt

#password
kWwSy?EQ3jrAh@VG

#enter the clusters for x hours
oarsub -l walltime=3:00 -I

#downloading our data (HCT-116 -> colorectal cancer cell lines)
#used "--split-files because we have paired data"
#3 conditions of Normoxia and 3 of Hypoxia
docker run -v $PWD/:/data/ -w=/data/ --rm ncbi/sra-tools:2.11.1 fasterq-dump --split-files SRR18643793
docker run -v $PWD/:/data/ -w=/data/ --rm ncbi/sra-tools:2.11.1 fasterq-dump --split-files SRR18643794
docker run -v $PWD/:/data/ -w=/data/ --rm ncbi/sra-tools:2.11.1 fasterq-dump --split-files SRR18643795
docker run -v $PWD/:/data/ -w=/data/ --rm ncbi/sra-tools:2.11.1 fasterq-dump --split-files SRR18643796
docker run -v $PWD/:/data/ -w=/data/ --rm ncbi/sra-tools:2.11.1 fasterq-dump --split-files SRR18643797
docker run -v $PWD/:/data/ -w=/data/ --rm ncbi/sra-tools:2.11.1 fasterq-dump --split-files SRR18643798

#after renaming the files, we aligned it with Kallisto to a reference human genome
docker run -v $PWD/:/data/ -w=/data/ --rm zlskidmore/kallisto:0.46.0 kallisto quant -i gencode.v39.transcriptome.idx -o Normoxia_rep3 HCT-116_Normoxia_rep3_1.fastq HCT-116_Normoxia_rep3_2.fastq
docker run -v $PWD/:/data/ -w=/data/ --rm zlskidmore/kallisto:0.46.0 kallisto quant -i gencode.v39.transcriptome.idx -o Normoxia_rep2 HCT-116_Normoxia_rep2_1.fastq HCT-116_Normoxia_rep2_2.fastq
docker run -v $PWD/:/data/ -w=/data/ --rm zlskidmore/kallisto:0.46.0 kallisto quant -i gencode.v39.transcriptome.idx -o Normoxia_rep1 HCT-116_Normoxia_rep1_1.fastq HCT-116_Normoxia_rep1_2.fastq
docker run -v $PWD/:/data/ -w=/data/ --rm zlskidmore/kallisto:0.46.0 kallisto quant -i gencode.v39.transcriptome.idx -o Hypoxia_rep3 HCT-116_Hypoxia_rep3_1.fastq HCT-116_Hypoxia_rep3_2.fastq
docker run -v $PWD/:/data/ -w=/data/ --rm zlskidmore/kallisto:0.46.0 kallisto quant -i gencode.v39.transcriptome.idx -o Hypoxia_rep2 HCT-116_Hypoxia_rep2_1.fastq HCT-116_Hypoxia_rep2_2.fastq
docker run -v $PWD/:/data/ -w=/data/ --rm zlskidmore/kallisto:0.46.0 kallisto quant -i gencode.v39.transcriptome.idx -o Hypoxia_rep1 HCT-116_Hypoxia_rep1_1.fastq HCT-116_Hypoxia_rep1_2.fastq