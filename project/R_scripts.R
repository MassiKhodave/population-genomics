
# Running tximport to combine files

library(tximportData)
library(tximport)

#locate the directory containing the files. 
dir <- "/data/project_data/brachy/transcripts_quant"
list.files(dir)

# read in table with sample ids
samples <- read.table("/data/project_data/brachy/transcripts_quant/brachy_sampleNO_HJ95.txt", header=TRUE)

# now point to quant files
all_files <- file.path(dir, samples$Index, "quant.sf")
names(all_files) <- samples$Index

# to be able to run without tx2gene
txi <- tximport(all_files, type = "salmon", txOut=TRUE)  
names(txi)

head(txi$counts)

countsMatrix <- txi$counts
dim(countsMatrix)


# To write out
write.table(countsMatrix, file = "/data/project_data/brachy/brachy_countsMatrix_NO_HJ95.txt", col.names = T, row.names = T	, quote = F)

#######################################################################################
# Running DESeq2

```
library(DESeq2)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(ggpubr)
library(wesanderson)
```

# reading the count table 
countsTable <- read.table(
    "C:\\Users\\jamsh\\OneDrive\\Desktop\\Population Genomics\\project\\brachy_countsMatrix_NO_HJ95.txt", header= TRUE, row.names=1)

# first few lines of the table, name of the columns, and the dimensions
head(countsTable)
names(countsTable)
dim(countsTable)

# rounding the counts to integers
countsTableRound <- round(countsTable)

# read the description text file
conds<-read.table("C:\\Users\\jamsh\\OneDrive\\Desktop\\Population Genomics\\project\\brachy_sampleNO_HJ95.txt", header=TRUE ,stringsAsFactors=TRUE, row.names=1)

head(conds)

# let's see how many reads we have from each sample
colSums(countsTableRound)
mean(colSums(countsTableRound))

png("counts_barplots.png",width = 30,height = 40)
dev.new()
par(mar=c(7,7,5,5))
barplot(colSums(countsTableRound), names.arg=colnames(countsTableRound), cex.names = 1, las=3, col=c("darkorange","darkorange","darkorange","darkseagreen", "darkseagreen","pink","pink","pink"))
abline(h=mean(colSums(countsTableRound)), col="blue", lwd=2)
dev.off()

#the average number of counts per gene
rowSums(countsTableRound)
mean(rowSums(countsTableRound)) #11930.81
median(rowSums(countsTableRound)) #4018.946

apply(countsTableRound,2,mean) #2 in the apply function does the action across columns
apply(countsTableRound,1,mean) #1 in the apply function does the action across rows
dev.new()
hist(apply(countsTableRound,1,mean),xlim=c(0,1000), breaks=10000)


## PCA

# create a DESeq object and define the experimental design here with the tilde 
dds <-DESeqDataSetFromMatrix(countsTableRound,
                             colData = conds, 
                             design =~ Genotype)
dim(dds)
# filter out genes with too few reads
#- keep with average >10 reads per sample
dds <-dds[rowSums(counts(dds)) >80]
dim(dds)
# Run the DESeq model to test for differential gene expression
dds <- DESeq(dds)

# List the results you've generated
resultsNames(dds)

# Let's start with a PCA to visualize global gene expression patterns
# We transform the output of the DESeq through variance stabilizing method for the visuzalization (and later clustering) purposes 
vsd <- vst(dds, blind=FALSE)

dev.new()
data <- plotPCA(vsd, intgroup=c("Genotype"), returnData=TRUE)
percentVar <- round(100 * attr(data,"percentVar"))

ggplot(data, aes(PC1,PC2, color=Genotype)) +
    geom_point(size=4, alpha=0.85) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) +
    theme_minimal()
dev.off()


#######################
#  TEST FOR EFFECT OF GENOTYPE
#######################

dds <- DESeqDataSetFromMatrix(countData = countsTableRound, 
                              colData = conds, 
                              design = ~ Genotype)

dds <- DESeq(dds, test="Wald")
# List the results you've generated
resultsNames(dds)

# Contrasting HDAC19 vs. WT
resHW <- results(dds, contrast = c("Genotype","HDAC19","Wild_Type"), alpha=0.05)
summary(resHW)
resHW <- resHW[order(resHW$padj),]
resHW <- resHW[!is.na(resHW$padj),]
degsHW <- row.names(resHW[resHW$padj < 0.05,]) 

LFCHW <- resHW[resHW$padj < 0.05,]$log2FoldChange
degsHW_down <- degsHW[LFCHW<0]
degsHW_up <- degsHW[LFCHW>0]

summary(degsHW_down) #248
summary(degsHW_up)   #126


# Contrasting GF14H vs. WT
resGW <- results(dds, contrast = c("Genotype","GF14H","Wild_Type"), alpha=0.05)
summary(resGW)
resGW <- resGW[order(resGW$padj),]
resGW <- resGW[!is.na(resGW$padj),]
degsGW <- row.names(resGW[resGW$padj < 0.05,]) 
LFCGW <- resGW[resGW$padj < 0.05,]$log2FoldChange
degsGW_down <- degsGW[LFCGW<0]
degsGW_up <- degsGW[LFCGW>0]

summary(degsGW_down) #76
summary(degsGW_up)   #98

###############################
# Volcano plot
###############################
#BiocManager::install('EnhancedVolcano')
library(EnhancedVolcano)
# plot adding up all layers we have seen so far

dev.new()
EnhancedVolcano(resGW,
                lab = rownames(resGW),
                x = 'log2FoldChange',
                y = 'padj',
                ylab='-Log adjusted P-value',
                title = 'GF14h vs. WT',
                pCutoff = 0.05,
                FCcutoff = 1.5,
                pointSize = 2.0,
                labSize = 3.5,
                shape = c(1, 4, 23, 25),
                colAlpha = 1,
                drawConnectors = TRUE,
                boxedLabels = TRUE)

dev.new()
EnhancedVolcano(resHW,
                lab = rownames(resHW),
                x = 'log2FoldChange',
                y = 'padj',
                ylab='-Log adjusted P-value',
                title = 'HDAC19 vs. WT',
                pCutoff = 0.05,
                FCcutoff = 1.5,
                pointSize = 2.0,
                labSize = 3.5,
                shape = c(1, 4, 23, 25),
                colAlpha = 1,
                drawConnectors = TRUE,
                boxedLabels = TRUE)



### Venn diagram

library(eulerr)

# Up-regulated genes

# number of common DEGs
nOV <- length(intersect(degsHW_up, degsGW_up))
# number of DEGs for HDAC19 vs WT
nHW <- length(degsHW_up)
# number of DEGs for GF14H vs WT
nGW <- length(degsGW_up)

# creating the Venn diagram using "euler"
fit1 <- euler(c("HDAC19 vs. WT" = nHW, "GF14H vs. WT" = nGW, "HDAC19 vs. WT&GF14H vs. WT" = nOV))

dev.new()
par(mar=c(15,15,10,10))
plot(fit1,  lty = 1:3, labels = list(cex=0), quantities = list(cex=5),fills = c("darkolivegreen3", "darkorchid2"))


# down-regulated genes

# number of common DEGs
nOV <- length(intersect(degsHW_down, degsGW_down))
# number of DEGs for HDAC19 vs WT
nHW <- length(degsHW_down)
# number of DEGs for GF14H vs WT
nGW <- length(degsGW_down)

# creating the Venn diagram using "euler"
fit1 <- euler(c("HDAC19 vs. WT" = nHW, "GF14H vs. WT" = nGW, "HDAC19 vs. WT&GF14H vs. WT" = nOV))

dev.new()
par(mar=c(10,10,18,18))
plot(fit1,  lty = 1:3, labels = list(cex=0), quantities = list(cex=5),fills = c("darkolivegreen3", "darkorchid2"))


### Plot Individual genes ### 

# Counts of specific top interaction gene! (important validatition that the normalization, model is working)
d <- plotCounts(dds, gene="Bradi1g08340", intgroup = (c("Genotype")), returnData=TRUE)
d

p <-ggplot(d, aes(x=Genotype, y=count, color=Genotype)) + 
    theme_minimal() + theme(text = element_text(size=20), panel.grid.major=element_line(colour="grey"))
p <- p + geom_point(position=position_jitter(w=0.2,h=0), size=6)
p <- p + stat_summary(fun = mean, geom = "line") + ggtitle("Bradi1g08340")
#p <- p + stat_summary(fun = mean, geom = "point")
dev.new()
p


# Heatmap of top 20 genes sorted by pvalue

library(pheatmap)

topgenes <- head(intersect(degsHW, degsGW),20)
mat <- assay(vsd)[topgenes,]
# make the VSD matrix zero-mean
mat <- mat - rowMeans(mat)
df <- as.data.frame(colData(dds)[c("Genotype")])

dev.new()
pheatmap(mat, annotation_col=df)


# plotting for one single gene (VRN1)
dev.new()
plotCounts(dds, gene=c("Bradi4g13970.1"), intgroup=c("Genotype"))

dev.new()
#plotCounts(dds, gene=head(rownames(resGW),1), intgroup=c("Genotype"))
plotCounts(dds, gene=topgenes[3], intgroup=c("Genotype"))



##################################################
# ENRICHMENT ANALYSIS
##################################################

annot.info <- read.csv("C:\\Users\\jamsh\\OneDrive\\Desktop\\Population Genomics\\project\\Bdistachyon_556_v3.2.annotation_info.txt", header=TRUE, sep="\t")

HW_GO = c()
for(i in 1:dim(resHW)[1]) { 
    goterm <- annot.info[annot.info$transcriptName==
                          row.names(resHW)[i],]$GO
    HW_GO <- append(HW_GO, gsub(" ",";",goterm))
}

HW_genes_wGO = row.names(resHW)[HW_GO!=""]
HW_GO_wGO = HW_GO[HW_GO!=""]

write.table(data.frame(HW_genes_wGO, HW_GO_wGO), file='C:\\Users\\jamsh\\OneDrive\\Desktop\\Population Genomics\\project\\GO\\HW_gene_GO.tsv', sep='\t',col.names=F, row.names=F,quote=F)

HW_LFC_wGO = resHW$log2FoldChange[HW_GO!=""]
write.table(data.frame(HW_genes_wGO, HW_LFC_wGO), file='C:\\Users\\jamsh\\OneDrive\\Desktop\\Population Genomics\\project\\GO\\HW_gene_LFC.tsv', sep='\t',col.names=F, row.names=F,quote=F)

HW_up_deflines = c()
for(i in 1:length(degsHW_up)) { 
    defline <- annot.info[annot.info$transcriptName==
                              degsHW_up[i],]$best_rice_defline
    HW_up_deflines <- append(HW_up_deflines, defline)
}


#######################################
# GO ANALYSIS
#######################################
library(dplyr)

# ------------- GF14h

# The names of genes of interest (DEGs), as character values
GF_geneID <- data.frame(row.names(resGW[,]))
GF_LFC <- data.frame(resGW[,2])
GF_geneID <- cbind(GF_geneID,GF_LFC)
names(GF_geneID) <- c("gene","logFoldChange")
write.csv(GF_geneID,file="GF14H_logFoldChange.csv",row.names=F)


# ----------- Gene Universe
# The gene Universe aka background set for GO analysis, which will contain the names of all 21897 transcripts
# Only need to make the geneUniverse.tab file once since it will be the same for both GF14h and HDAC19
transcriptName <- GF_geneID %>%
    select(1)
names(transcriptName)[1] <- "transcriptName"

# Read in the annotation file (with GO terms)
annotations <- read.delim("Bdistachyon_556_v3.2.annotation_info.txt")

# Selecting only the columns we need,transcript names and GO categories
GoMap <- annotations %>% 
    select(transcriptName,GO)

transcriptName <- merge(transcriptName,GoMap,all.x=T)


a=""
for (x in c(0:length(transcriptName[,2]))) {
    transcriptName[x,2] <- gsub(" ",";",transcriptName[x,2])
    if(identical(a,transcriptName[x,2])){
        transcriptName[x,2] <- gsub("","unknown",transcriptName[x,2])
    }
}

write.table(transcriptName, file="geneUniverse.tab", sep="\t", row.names = F,col.names = F)

# ------------- HDAC19
# The names of genes of interest (DEGs), as character values
HD_geneID <- data.frame(row.names(resHW[,]))
HD_LFC <- data.frame(resHW[,2])
HD_geneID <- cbind(HD_geneID,HD_LFC)
names(HD_geneID) <- c("gene","logFoldChange")
write.csv(HD_geneID,file="HDAC19_logFoldChange.csv",row.names=F)

