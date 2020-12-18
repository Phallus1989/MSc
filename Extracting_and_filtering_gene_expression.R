#Gene expression count extraction

#William Taylor 18/12/2020
#wil_tr@live.co.uk




#Including ENTREZID for Set Gene ontology and Kegg analysis
setwd("/C:")
library(tximportData)
library(tximport)
library(GenomicFeatures)
library(edgeR)
library(limma)
library(readr)
library(rjson)
library(csaw)
library(data.table)
library(Glimma)
library(textshape)
library(org.Hs.eg.db)

#making TxDb object with genome
gffFile <- system.file("extdata","GTF_files","GCF_000001405.38_GRCh38.p12_genomic.gff",package="GenomicFeatures")
txdb <- makeTxDbFromGFF(file=gffFile,
                        dataSource="gff file",
                        organism="Homo sapiens")

#using tximport to take quant.sf files from salmon 
dir <- system.file("extdata", package = "tximportData")
list.files(dir)

k = keys(txdb,keytype = "GENEID")
df <- select(txdb, keys = k,  columns = "TXNAME", keytype = "GENEID")
tx2gene <- df[, 2:1]
dir <- system.file("extdata", package = "tximportData")
list.files(dir)
k = keys(txdb, keytype = "GENEID")
df <- select(txdb, keys = k,  columns = "TXNAME", keytype = "GENEID")
tx2gene <- df[, 2:1]
files <- list.files( pattern = "quant.sf",full.names = TRUE)

#make sure the sample names match up, and if they dont, link them to the appropriate ones in metadata
write.table(files,"samplematchingInfo_RNA-unStranded.txt",sep = "\t",col.names = TRUE, quote = F )
all(file.exists(files))

#calculating the gene counts file and exporting count values
txi_counts = tximport::tximport(files, type = "salmon",txIn = TRUE, tx2gene = tx2gene,dropInfReps = TRUE)
write.csv(txi_counts$counts, "counts.csv")

count <- rownames_to_column(counts)
list <- count[1]
list<- lapply(list, as.character)
hs <- org.Hs.eg.db
counter_list <- mapIds(hs, list$rowname, "ENTREZID", "SYMBOL")
list2 <- as.data.frame(counter_list)
list2 <- rownames_to_column(list2)
colnames(list2) <- c("SYMBOL", "ENTREZID")
x <- DGEList(counts = counts, lib.size = colSums(counts), norm.factors =calcNormFactors(counts, method = "TMM"), samples = des, group = des$Tissue, genes = list2, remove.zeros = FALSE)
## Name dge column names 
samplenames <- c("RT1T", "RT10T", "RT11T", "RT12T", "RT13T", "RT14T", "RT15T", "RT16T", "RT17T", "RT18T", "RT19T", "RT2T", "RT20T", "RT1N", "RT2N", "RT3N", "RT4N", "RT5N", "RT6N", "RT7N", "RT8N", "RT9N", "RT3T", "RT10N", "RT11N", "RT12N", "RT13N", "RT14N", "RT15N", "RT16N", "RT17N", "RT18N", "RT19N", "RT4T", "RT21T", "RT22T", "RT23T", "RT24T", "RT25T", "RT26T", "RT27T", "RT28T", "RT29T", "RT30T", "RT5T", "RT31T", "RT32T", "RT33T", "RT34T", "RT35T", "RT36T", "RT37T", "RT38T", "RT39T", "RT40T", "RT6T", "RT21N", "RT22N", "RT23N", "RT24N", "RT25N", "RT26N", "RT27N", "RT28N", "RT29N", "RT30N", "RT7T", "RT31N", "RT32N", "RT33N", "RT34N", "RT35N", "RT36N", "RT37N", "RT38N", "RT39N", "RT40N", "RT8T", "RT20N", "RT9T")
colnames(x)<-samplenames
summary(x)
x
o <- order(rowSums(x$counts), decreasing = TRUE)
x <- x[o,]
dim(x)
#38461    80
dev.off()
#exploration
plotMDS(cpm(x, log = TRUE), main = "MDS plot of logFC between samples")


cpm <- cpm(x)
lcpm <- cpm(x, log = TRUE)
write.csv(lcpm, "lcpm_counts.csv")
write.csv(cpm, "cpm_counts.csv")
L <- mean(x$samples$lib.size) * 1e-6
M <- median(x$samples$lib.size) * 1e-6
c(L, M)
# 36.14167 35.06695
#check number of 0 count genes across samples
table(rowSums(x$counts==0)==80)
#FALSE  TRUE 
#36622  1839  
summary(lcpm)
keep <- filterByExpr(x, group = des$Tissue, min.count = 100)
x <- x[keep, keep.lib.sizes = FALSE]
dim(x)
# 13613       80
lcpm <- cpm(x, log = TRUE)
write.csv(x$counts, "genecounts.csv")
