setwd("C:/Users/jamsh/OneDrive/Desktop/Population Genomics/Chr10_sweeps/")

library(ggplot2)
library(gridExtra)

# Get the list of admixed individuals:
Admixed <- read.table("Admixed.Inds",header=F)

# Read in list of positions
snps <- read.table("Chr10.kept.sites",sep="\t", header=T)


# loading Plink results that were generated for climate traits
cGDDfreeze <- read.table("plink2.mean_cGDDfreeze.glm.linear",skip=1,sep="\t",header=F)
names(cGDDfreeze) = c("CHROM",    "POS",  "ID",   "REF",  "ALT",  "A1",   "TEST", "OBS_CT",   "BETA", "SE",   "T_STAT",   "P")
cGDDfreeze2 <- cGDDfreeze[which(cGDDfreeze$TEST=="ADD"),]


# Define association outliers as the upper 0.1% of p-values
cGDDfreeze2 <- cbind(snps, cGDDfreeze2[,-c(1:2)])
cGDDfreeze2$outlier = ifelse(cGDDfreeze2$P<quantile(cGDDfreeze2$P,0.01),2,1)

p1 <- ggplot(cGDDfreeze2,aes(x=POS,y=-1*log10(P))) +
  geom_point(size=cGDDfreeze2$outlier, color=cGDDfreeze2$outlier) + 
  xlab("Position (bp) along chromosome") +
  ylab("-log10 P-value") +
  ggtitle("Growing Degree Days")

#######################################
# repeat above code for the next climate variable: med_DD0
DD0 <- read.table("plink2.med_DD0.glm.linear",skip=1,sep="\t",header=F)
names(DD0) = c("CHROM",  "POS",  "ID",   "REF",  "ALT",  "A1",   "TEST", "OBS_CT",   "BETA", "SE",   "T_STAT",   "P")
DD02 <- DD0[which(DD0$TEST=="ADD"),]
DD02 <- cbind(snps, DD02[,-c(1,2)])
DD02$outlier = ifelse(DD02$P<quantile(DD02$P,0.01),2,1)

p2 <- ggplot(DD02,aes(x=POS,y=-1*log10(P))) +
  geom_point(size=DD02$outlier, color=DD02$outlier) + 
  xlab("Position (bp) along chromosome") +
  ylab("-log10 P-value") +
  ggtitle("Chilling Degree Days")

#######################################
# repeat above code for the next climate variable: mean_finalFreeze
finalFreeze <- read.table("plink2.mean_finalFreeze.glm.linear",skip=1,sep="\t",header=F)
names(finalFreeze) = c("CHROM",    "POS",  "ID",   "REF",  "ALT",  "A1",   "TEST", "OBS_CT",   "BETA", "SE",   "T_STAT",   "P")
finalFreeze2 <- finalFreeze[which(finalFreeze$TEST=="ADD"),]
finalFreeze2 <- cbind(snps, finalFreeze[,-c(1,2)])
finalFreeze2$outlier = ifelse(finalFreeze2$P<quantile(finalFreeze2$P,0.01),2,1)

p3 <- ggplot(finalFreeze2,aes(x=POS,y=-1*log10(P))) +
  geom_point(size=finalFreeze2$outlier, color=finalFreeze2$outlier) + 
  xlab("Position (bp) along chromosome") +
  ylab("-log10 P-value") +
  ggtitle("Last Freezing")

# png(heigh=50,width=30, 'snps_plink_clim.png')
dev.new()
grid.arrange(p1, p2, p3, p4, nrow = 4)



library(GenomicRanges)
library(GenomicFeatures)

CHR = "Chr10"

budflushGR <- GRanges(CHR,IRanges(budflush2$POS-2.5e4,budflush2$POS+2.5e4),POS=budflush2$POS, P=budflush2$P, outlier=budflush2$outlier)

# subsetting reduced set of windows to capturethe outlier regions
budflushGRout <- unlist(reduce(split(budflushGR, ~outlier)))
budflushGRout$outlier <- names(budflushGRout)
budflushGRCand <- subset(budflushGRout, outlier==2)

budflushGRCand # Print the candidate regions


cGDDfreezeGR <- GRanges(CHR,IRanges(cGDDfreeze2$POS-2.5e4,cGDDfreeze2$POS+2.5e4),POS=cGDDfreeze2$POS, P=cGDDfreeze2$P, outlier=cGDDfreeze2$outlier)
cGDDfreezeGRout <- unlist(reduce(split(cGDDfreezeGR, ~outlier)))
cGDDfreezeGRout$outlier <- names(cGDDfreezeGRout)
cGDDfreezeGRCand <- subset(cGDDfreezeGRout, outlier==2)


DD0GR <- GRanges(CHR,IRanges(DD02$POS-2.5e4,DD02$POS+2.5e4),POS=DD02$POS, P=DD02$P, outlier=DD02$outlier)
DD0GRout <- unlist(reduce(split(DD0GR, ~outlier)))
DD0GRout$outlier <- names(DD0GRout)
DD0GRCand <- subset(DD0GRout, outlier==2)


finalFreezeGR <- GRanges(CHR,IRanges(finalFreeze2$POS-2.5e4,finalFreeze2$POS+2.5e4),POS=finalFreeze2$POS, P=finalFreeze2$P, outlier=finalFreeze2$outlier)
finalFreezeGRout <- unlist(reduce(split(finalFreezeGR, ~outlier)))
finalFreezeGRout$outlier <- names(finalFreezeGRout)
finalFreezeGRCand <- subset(finalFreezeGRout, outlier==2)


overlap_BF_cGDDfReeze <- subsetByOverlaps(budflushGRCand, cGDDfreezeGRCand)
length(overlap_BF_cGDDfReeze)

overlap_BF_DD0 <- subsetByOverlaps(budflushGRCand, DD0GRCand)
length(overlap_BF_DD0)

overlap_BF_fFReeze <- subsetByOverlaps(budflushGRCand, finalFreezeGRCand)
length(overlap_BF_fFReeze)

write.table(candGenes$geneID, paste0("Annotation/candGenes",CHR,".txt"), quote=F, col.names=F, row.names=F, sep=",")


#################
# extracting genes
# Import the GFF annotation file and make a transcript database
txdb <- makeTxDbFromGFF("Ptrichocarpa_533_v4.1.gene.gff3.gz", format="gff3")

# How many chromosomes are present?
head(seqlevels(txdb))

# Subset the database for just your chromosome of interest
seqlevels(txdb) <- CHR # subset for just your chromosome

# Reduce the transcript database to just the non-redundant gene names, instead of multiple entries for all the variant transcript types per gene
genes <- unlist(reduce(transcriptsBy(txdb, by="gene"))) 
genes$geneID <- names(genes)

candGenes_BF_fFReeze <- subsetByOverlaps(genes, overlap_BF_fFReeze)
candGenes_BF_DD0 <- subsetByOverlaps(genes, overlap_BF_DD0)
candGenes_BF_cGDDfFreeze <- subsetByOverlaps(genes, overlap_BF_cGDDfReeze)


write.table(candGenes_BF_cGDDfFreeze$geneID, paste0("candGenes_BF_cGDDfFreeze",CHR,".txt"), quote=F, col.names=F, row.names=F, sep=",")
#
write.table(candGenes_BF_DD0$geneID, paste0("candGenes_BF_DD0",CHR,".txt"), quote=F, col.names=F, row.names=F, sep=",")
#
write.table(candGenes_BF_fFReeze$geneID, paste0("candGenes_BF_fFreeze",CHR,".txt"), quote=F, col.names=F, row.names=F, sep=",")

