# Get all necessary files: raw data and the genome
# we are using our own ChIPseq data (we could not find the input data) https://www.ncbi.nlm.nih.gov/sra?term=SRX14746757
        
#get our reps
docker run -v $PWD/:/data/ -w=/data/ --rm ncbi/sra-tools:2.11.1 fasterq-dump SRR18643593

#copy the reference from the shared folder
cp -r /mnt/share/acmo/AppCOmics_data/Genomes/GRCh38_noalt_as .

#align to the reference
docker run -v $PWD:$PWD -w=$PWD --rm staphb/bowtie2:latest bowtie2 -k 1 --threads 10 -x ./GRCh38_noalt_as/GRCh38_noalt_as -1 SRR18643592_1.fastq,SRR18643593_1.fastq -2 SRR18643592_2.fastq,SRR18643593_2.fastq -S HIF1A_ChIPseq.sam

#results
"""
27891617 reads; of these:
27891617 (100.00%) were paired; of these:
2418042 (8.67%) aligned concordantly 0 times
25176194 (90.26%) aligned concordantly exactly 1 time
297381 (1.07%) aligned concordantly >1 times
----
2418042 pairs aligned concordantly 0 times; of these:
1081834 (44.74%) aligned discordantly 1 time
----
1336208 pairs aligned 0 times concordantly or discordantly; of these:
2672416 mates make up the pairs; of these:
1690073 (63.24%) aligned 0 times
935827 (35.02%) aligned exactly 1 time
46516 (1.74%) aligned >1 times
"""
 
#convert the output file in bam format and sort (i.e. sort the reads by genomic position) using sambamba tool
docker run -v $PWD:$PWD -w=$PWD --rm miguelpmachado/sambamba:0.7.1-01 sambamba view -t 10 -S -f bam HIF1A_ChIPseq.sam -o temp.bam
docker run -v $PWD:$PWD -w=$PWD --rm miguelpmachado/sambamba:0.7.1-01 sambamba sort -t 10 -o HIF1A_ChIPseq.bam temp.bam


#identification of enriched ChIPseq regions (i.e. peaks) comparing the ChIP and input samples with MACS2 (removed the -c command because we have no input)
docker run -v $PWD:$PWD -w=$PWD --rm resolwebio/chipseq:5.1.3 macs2 callpeak -t HIF1A_ChIPseq.bam -f BAM -g 2.7e9 -q 0.05 -n HIF1A --outdir macs2

#Filter out blacklisted regions from ChIPseq results
#we need to intersect the bed file we generated with the blacklist (summits)
docker run -v $PWD:$PWD -w=$PWD --rm biocontainers/bedtools:v2.27.1dfsg-4-deb_cv1 bedtools intersect -v -a HIF1A_summits.bed -b hg38.blacklist.bed > Intersect.bed

#get the file genomeSize.txt (text file containing the size of each human chromosome) - **used in the first step**
cp -r /mnt/share/acmo/AppCOmics_data/Genomes/genomeSize.txt .
    
#get the file GCA_000001405.15_GRCh38_no_alt_analysis_set.fna (entire human genome sequence) - **used in the second step**
cp -r /mnt/share/acmo/AppCOmics_data/Genomes/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna .
cp -r /mnt/share/acmo/AppCOmics_data/Genomes/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fai .

#get the file gencode.v26.annotation_promoters.fasta (all gene promoters as background sequences) - **used in the third step**
cp -r /mnt/share/acmo/AppCOmics_data/Genomes/gencode.v26.annotation_promoters.fasta .

# **First:** extend the genomic coordinates of the peak summit by 200bp upstream and downstream using the function “slop” from BedTools. The output of the function should be a bed file slop https://bedtools.readthedocs.io/en/latest/content/tools/slop.html?highlight=slop
docker run -v $PWD:$PWD -w=$PWD --rm biocontainers/bedtools:v2.27.1dfsg-4-deb_cv1 bedtools slop -i Intersect.bed -g genomeSize.txt -b 200 > extension.bed

#**Second:** you get the DNA sequences for each peak summit using the function “getFasta” from Bedtools. The output of the function should be a fasta file (extension .fasta) https://bedtools.readthedocs.io/en/latest/content/tools/getfasta.html?highlight=getfasta
docker run -v $PWD:$PWD -w=$PWD --rm biocontainers/bedtools:v2.27.1dfsg-4-deb_cv1 bedtools getfasta -fi GCA_000001405.15_GRCh38_no_alt_analysis_set.fna -bed extension.bed > extension.fasta

#**Third:** perform the motif enrichment analysis using the “findMotifs” function from HOMER tool (Docker image nfcore/chipseq:latest)  http://homer.ucsd.edu/homer/motif/fasta.html
docker run -v $PWD:$PWD -w=$PWD --rm nfcore/chipseq:latest findMotifs.pl extension.fasta fasta homer_results.fasta -fasta gencode.v26.annotation_promoters.fasta