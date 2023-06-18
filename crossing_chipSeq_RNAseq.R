# 1st Part: get the annotations for our Chip_seq peaks ####

BiocManager::install("ChIPpeakAnno")
library("ChIPpeakAnno")
library("GenomicRanges")

#read in the bed file
df1<-read.table("HIF1A_filtered_peaks.bed", header=FALSE)

#convert the peaks to a GRanges object
gr1 <- GRanges(seqnames=df1$V1, ranges=IRanges(start=df1$V2, end=df1$V3))

##annotation
BiocManager::install("EnsDb.Hsapiens.v86")
library("EnsDb.Hsapiens.v86")

## create annotation file from EnsDb (Ensembl) or TxDb (transcript annotation) packages
annoData <- toGRanges(EnsDb.Hsapiens.v86, feature="gene")
annoData[1:2]

# annotate the peak GRanges object with peaks mapped to gene with a -2000 and 500 bp window around the TSS
anno.gr1 <- annotatePeakInBatch(gr1, 
                                AnnotationData=annoData, 
                                output="nearestBiDirectionalPromoters",
                                bindingRegion=c(-2000, 500))

#trim out of bound ranges
anno.gr1 <- trim(anno.gr1)

#annotate with Entrez IDs
# BiocManager::install("org.Hs.eg.db")
library("org.Hs.eg.db")
anno.gr1 <- addGeneIDs(anno.gr1,"org.Hs.eg.db",IDs2Add = "entrez_id")

# list annotated peaks
head(anno.gr1)

# get the gene names from your annotated data
annotated_genes <- anno.gr1$gene_name

# 2nd Part: get the gene names from differentiated expressed genes (RNAseq) ####
library(readxl)
DEGs <- read_excel("DEGs.xlsx")

diffExp <- as.list(DEGs$...1)
upReg <- as.list(read.table("upReg.txt"))
upReg <- upReg[["V1"]]
downReg <- as.list(read.table("downReg.txt"))
downReg <- downReg[["V1"]]

# 3rd Part: intersect results #####

annotated_upReg <- intersect(annotated_genes, upReg)

library(VennDiagram)

venn.diagram(list("Annotated"=annotated_genes, "DEGs"=diffExp), fill = c("yellow","cyan"), cex = 1, filename="Crossing_chip_Rna/Venn_diagram_annotated_DEGs.png")

venn.diagram(list("Annotated"=annotated_genes, "Upregulated"=upReg, "Downregulated"=downReg), fill = c("yellow","red","cyan"), cex = 1.5,filename="Crossing_chip_Rna/Venn_diagram_annotated_up_down.png")

venn.diagram(
  list("Annotated" = annotated_genes, "Upregulated" = upReg, "Downregulated" = downReg),
  fill = c("yellow", "red", "cyan"),
  cex = 1.5,
  filename = "Crossing_chip_Rna/Venn_diagram_annotated_up_down.png",
  margin = 0.1  # Adjust the margin value as per your requirement
)

