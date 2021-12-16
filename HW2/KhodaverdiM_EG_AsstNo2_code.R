library(DESeq2)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(ggpubr)
library(wesanderson)

# reading the count table for a generation (either F1, which we already did in class, or F3)
countsTable <- read.table(
  "C:\\Users\\jamsh\\OneDrive\\Desktop\\Population Genomics\\DE_counts_F3.txt", header= TRUE, row.names=1)

# first few lines of the table, name of the columns, and the dimensions
head(countsTable)
names(countsTable)
dim(countsTable)

# rounding the counts to integers
countsTableRound <- round(countsTable)

# read the description text file
conds<-read.delim("C:\\Users\\jamsh\\OneDrive\\Desktop\\Population Genomics\\RT_tonsa_F3_Samples.txt", header=TRUE ,stringsAsFactors=TRUE, row.names=1)

head(conds)

# let's see how many reads we have from each sample
colSums(countsTableRound)
mean(colSums(countsTableRound))
barplot(colSums(countsTableRound), names.arg=colnames(countsTableRound), cex.names = 0.5, las=3, ylim=c(0,2))
abline(h=mean(colSums(countsTableRound)), col="blue", lwd=2)


#the average number of counts per gene
rowSums(countsTableRound)
mean(rowSums(countsTableRound)) #11930.81
median(rowSums(countsTableRound)) #2226

apply(countsTableRound,2,mean) #2 in the apply function does the action across columns
apply(countsTableRound,1,mean) #1 in the apply function does the action across rows
hist(apply(countsTableRound,1,mean,xlim=c(0,1000), ylim=c(0,10), breaks=10000))


#### create a DESeq object and define the experimental design here with the tilde 

dds <-DESeqDataSetFromMatrix(countsTableRound, colData = conds, 
                             design =~ line+ environment+line:environment)
dim(dds)
# filter out genes with too few reads
#- keep with average >10 reads per sample
dds <-dds[rowSums(counts(dds)) >160]
dim(dds)
# Run the DESeq model to test for differential gene expression
dds <- DESeq(dds)

# List the results you've generated
resultsNames(dds)
# [1] "Intercept"                  "line_combined_vs_ambient"   "environment_HH_vs_AA"      
# [4] "linecombined.environmentHH"
# Let's start with a PCA to visualize global gene expression patterns
vsd <- vst(dds, blind=FALSE)

data <- plotPCA(vsd, intgroup=c("line","environment"), returnData=TRUE)
percentVar <- round(100 * attr(data,"percentVar"))

ggplot(data, aes(PC1,PC2, color=environment, shape=line)) +
  geom_point(size=4, alpha=0.85) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  theme_minimal()
dev.off()
# What patterns do we see? Clustering by groups, line and environment. 
# What gene expression results do we expect for each factor, main effects and/or interactions?
###############################################################


###### Continue on dds model from before the PCA
###### Order and summarize the results from specific contrasts

resInteraction <- results(dds, alpha=0.05)
resInteraction <- resInteraction[order(resInteraction$padj),]
head(resInteraction) 
#######################
############################################## TEST FOR EFFECT OF ENVIRONMENT
#######################

dds <- DESeqDataSetFromMatrix(countData = countsTableRound, colData = conds, 
                              design = ~ line + environment)

dds <- DESeq(dds, test="LRT", reduced=~line)
# List the results you've generated
resultsNames(dds)

# Order and list and summarize results from specific contrasts
resEnv <- results(dds, alpha = 0.05)
resEnv <- resEnv[order(resEnv$padj),]
head(resEnv)

summary(resEnv)

resEnv <- resEnv[!is.na(resEnv$padj),]

degsEnv <- row.names(resEnv[resEnv$padj < 0.05,]) 
#######################
##############################################  TEST FOR EFFECT OF LINE
#######################

dds <- DESeqDataSetFromMatrix(countData = countsTableRound, colData = conds, 
                              design = ~ environment + line)

dds <- DESeq(dds, test="LRT", reduced=~environment)
resultsNames(dds)

resLine <- results(dds, alpha = 0.05)
resLine <- resLine[order(resLine$padj),]
head(resLine)


summary(resLine)


resLine <- resLine[!is.na(resLine$padj),]

degsline <- row.names(resLine[resLine$padj < 0.05,])

#######################
##############################################  TEST FOR INTERACTION
#######################

dds <- DESeqDataSetFromMatrix(countData = countsTableRound, colData = conds, 
                              design = ~ environment + line + environment:line)

dds <- DESeq(dds, test="LRT", reduced=~environment + line)
resultsNames(dds)

resInt <- results(dds, alpha = 0.05)
resInt <- resInt[order(resInt$padj),]
head(resInt)

summary(resInt)


resInt <- resInt[!is.na(resInt$padj),]

degsInt <- row.names(resInt[resInt$padj < 0.05,])


### Plot Individual genes ### 

# Counts of specific top interaction gene! (important validatition that the normalization, model is working)
d <-plotCounts(dds, gene="TRINITY_DN138549_c1_g2", intgroup = (c("line","environment")), returnData=TRUE)
d

p <-ggplot(d, aes(x=environment, y=count, color=line, shape=line, group=line)) + 
  theme_minimal() + theme(text = element_text(size=20), panel.grid.major=element_line(colour="grey"))
p <- p + geom_point(position=position_jitter(w=0.2,h=0), size=3)
p <- p + stat_summary(fun = mean, geom = "line")
p <- p + stat_summary(fun = mean, geom = "point", size=5, alpha=0.7) 
p

#######################
############################################## PLOT OVERLAPPING DEGS IN VENN DIAGRAM
#######################

library(eulerr)

# Total
length(degsEnv)  # 448
length(degsline)  # 226
length(degsInt)  # 3854

# Intersections
length(intersect(degsEnv,degsline))  # 37
length(intersect(degsEnv,degsInt))  # 44
length(intersect(degsInt,degsline))  # 34

intEL <- intersect(degsEnv,degsline)
length(intersect(degsInt,intEL)) # 7

# Number unique
448-44-37-7 # 360
226-37-34-7 # 148
3854-44-34-7 # 3769


fit1 <- euler(c("Env" = 360, "Line" = 148, "Interaction" = 3769, "Env&Line" = 37, "Env&Interaction" = 44, "Line&Interaction" = 34, "Env&Line&Interaction" = 7))

plot(fit1,  lty = 1:3, quantities = TRUE)

plot(fit1, quantities = TRUE, fill = "transparent",
     lty = 1:3,
     labels = list(font = 4))

# Heatmap of top 20 genes sorted by pvalue

library(pheatmap)

# Interaction

topgenes <- head(row.names(resInt),20)
mat <- assay(vsd)[topgenes,]
mat <- mat - rowMeans(mat)
df <- as.data.frame(colData(dds)[,c("line","environment")])
pheatmap(mat, annotation_col=df)

# By line

topgenes <- head(rownames(resLine),20)
mat <- assay(vsd)[topgenes,]
mat <- mat - rowMeans(mat)
df <- as.data.frame(colData(dds)[,c("line","environment")])
pheatmap(mat, annotation_col=df)


#################
# VENN DIAGRAM OF F1 vs. F3 LINE DEGs

# I already saved DEGs for F1 and F3 obtained above 
degsLineF1 <- read.csv("C:\\Users\\jamsh\\OneDrive\\Desktop\\Population Genomics\\F1\\degs_line.csv", header=FALSE)
degsLineF3 <- read.csv("C:\\Users\\jamsh\\OneDrive\\Desktop\\Population Genomics\\F3\\degs_line.csv", header=FALSE)

# intersection between the DEG2
OV <- dim(intersect(degsLineF1,degsLineF3))[1]
# number of DEGs for F1
SF1 <- dim(degsLineF1)[1]
# number of DEGs for F3
SF3 <- dim(degsLineF3)[1]

# creating the Venn diagram using "euler"
fit1 <- euler(c("Gen.1" = SF1, "Gen.3" = SF3, "Gen.1&Gen.3" = OV))
# plotting and saving
png(filename="C:\\Users\\jamsh\\OneDrive\\Desktop\\Population Genomics\\venn.png", width=450, height=350)
plot(fit1,  lty = 1:3, labels = list(cex=2), quantities = list(cex=1.5))
dev.off()
