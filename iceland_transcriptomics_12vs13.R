library(DESeq2)
library(RColorBrewer)
library(pheatmap)
library(tidyverse)
library(gplots)
library(plyr)
library(affy)
library(glmpca)
library(knitr)
library(gprofiler2)


# setwd("../Desktop/iceland_metat/")
data <- read.csv("iceland_genes_new.tsv", sep = "\t", header = T)
metadata <-  read.csv("metadata2.txt", sep = "\t")

dat1 <- ddply(data,"Gene.name",numcolwise(sum))


rownames(dat1) <- dat1[,1]
dat1 <- dat1[,-1]
# "10", "10", "10", "11", "11", "11", "12", "12", "12", "13", "13", "13", "14", "14", "14", "Neg", "Neg", "Neg

# condition =  factor(c("Cultive", "10", "10", "10", "11", "11", "11", "12", "12", "12", "13", "13", "13", "14", "14", "14", "Neg", "Neg", "Neg"))
# coldata <- data.frame(row.names=colnames(dat1), condition)


rownames(metadata)
colnames(dat1)

all(rownames(metadata) == colnames(dat1))


###### ##################################### sample comparison_12vs13



head(dat1)
# "Cultive", "10", "10", "10", "11", "11", "11", "12", "12", "12", "13", "13", "13", "14", "14", "14", "Neg", "Neg", "Neg


### subsetting count data and metadata for the sample sites 
dat1_12vs13 <- dat1 %>% select(starts_with(c("X12","X13")))

metadata_12vs13 <- metadata[metadata$Site %in% c(12, 13), ]

## Deseq Matrix
dds_12vs13 <- DESeqDataSetFromMatrix(countData= dat1_12vs13,
                                     colData = metadata_12vs13,
                                     design = ~ 1)

# removing NA rows

keep_12vs13 <- rowSums(counts(dds_12vs13)) > 1
dds_12vs13 <- dds_12vs13[keep_12vs13,]
nrow(dds_12vs13)

dds_12vs13 <- estimateSizeFactors(dds_12vs13)



# Count variance is related to mean

## Computing mean and variance
norm.counts_12vs13 <- counts(dds_12vs13, normalized=TRUE)
mean.counts_12vs13 <- rowMeans(norm.counts_12vs13)
variance.counts_12vs13 <- apply(norm.counts_12vs13, 1, var)

## sum(mean.counts==0) # Number of completely undetected genes

norm.counts.stats_12vs13 <- data.frame(
  min=apply(norm.counts_12vs13, 2, min),
  mean=apply(norm.counts_12vs13, 2, mean),
  median=apply(norm.counts_12vs13, 2, median),
  max=apply(norm.counts_12vs13, 2, max),
  zeros=apply(norm.counts_12vs13==0, 2, sum),
  percent.zeros=100*apply(norm.counts_12vs13==0, 2, sum)/nrow(norm.counts_12vs13),
  perc05=apply(norm.counts_12vs13, 2, quantile, 0.05),
  perc10=apply(norm.counts_12vs13, 2, quantile, 0.10),
  perc90=apply(norm.counts_12vs13, 2, quantile, 0.90),
  perc95=apply(norm.counts_12vs13, 2, quantile, 0.95)
)

kable(norm.counts.stats_12vs13)


## Mean and variance relationship
mean.var.col_12vs13 <- densCols(x=log2(mean.counts_12vs13), y=log2(variance.counts_12vs13))

png("Mean_and variance_relationship_12vs13.png", width = 8, height = 4, units = 'in', res = 300)

plot(x=log2(mean.counts_12vs13), y=log2(variance.counts_12vs13), pch=16, cex=0.5, 
     col=mean.var.col_12vs13, main="Mean-variance relationship",
     xlab="Mean log2(normalized counts) per gene",
     ylab="Variance of log2(normalized counts)",
     panel.first = grid())
abline(a=0, b=1, col="brown")

dev.off()
## Modelling read counts through a negative binomial

## Performing estimation of dispersion parameter
dds.disp_12vs13 <- estimateDispersions(dds_12vs13)

## A diagnostic plot which
## shows the mean of normalized counts (x axis)
## and dispersion estimate for each genes
png("Dispersion_12vs13.png", width = 8, height = 4, units = 'in', res = 300)
plotDispEsts(dds.disp_12vs13)
dev.off()

##### Performing differential expression call

alpha <- 0.0001
wald.test_12vs13 <- nbinomWaldTest(dds.disp_12vs13)
res.DESeq2_12vs13 <- results(wald.test_12vs13, alpha=alpha, pAdjustMethod="BH")

res.DESeq2_12vs13 <- res.DESeq2_12vs13[order(res.DESeq2_12vs13$padj),]

## Draw an histogram of the p-values
png("histogram_12vs13.png", width = 8, height = 4, units = 'in', res = 300)
hist(res.DESeq2_12vs13$padj, breaks=20, col="grey", main="DESeq2 p-value distribution", xlab="DESeq2 P-value", ylab="Number of genes")
dev.off()
###########Volcano plot

alpha <- 0.01 # Threshold on the adjusted p-value
cols_12vs13 <- densCols(res.DESeq2_12vs13$log2FoldChange, -log10(res.DESeq2_12vs13$pvalue))

png("volcano_12vs13.png", width = 8, height = 4, units = 'in', res = 300)
plot(res.DESeq2_12vs13$log2FoldChange, -log10(res.DESeq2_12vs13$padj), col=cols_12vs13, panel.first=grid(),
     main="Volcano plot", xlab="Effect size: log2(fold-change)", ylab="-log10(adjusted p-value)",
     pch=20, cex=0.6)
abline(v=0)
abline(v=c(-1,1), col="brown")
abline(h=-log10(alpha), col="brown")
gn.selected_12vs13 <- abs(res.DESeq2_12vs13$log2FoldChange) > 2 & res.DESeq2_12vs13$padj < alpha 
text(res.DESeq2_12vs13$log2FoldChange[gn.selected_12vs13],
     -log10(res.DESeq2_12vs13$padj)[gn.selected_12vs13],
     lab=rownames(res.DESeq2_12vs13)[gn.selected_12vs13 ], cex=0.4)
dev.off()

### Check the expression levels of the most differentially expressed gene

gn.most.sign_12vs13 <- rownames(res.DESeq2_12vs13)[1]
gn.most.diff.val_12vs13 <- counts(dds_12vs13, normalized=T)[gn.most.sign_12vs13,]

png("DFgenes_12vs13.png", width = 8, height = 4, units = 'in', res = 300)
barplot(gn.most.diff.val_12vs13, col=metadata_12vs13, main=gn.most.sign_12vs13, las=2, cex.names=0.5)
dev.off()

### MA plot

## Draw a MA plot.
## Genes with adjusted p-values below 1% are shown
png("MA_12vs13.png", width = 8, height = 8, units = 'in', res = 300)

plotMA(res.DESeq2_12vs13, colNonSig = "blue")
abline(h=c(-1:1), col="red")
dev.off()


# functional enrichment analysis using the list of induced genes. 
# This step will be performed using the gProfileR R library.


epsilon <- 1 # pseudo-count to avoid problems with log(0)

res.DESeq2.df_12vs13 <- na.omit(data.frame(res.DESeq2_12vs13))
induced.sign_12vs13 <- rownames(res.DESeq2.df_12vs13)[res.DESeq2.df_12vs13$log2FoldChange >= 2 &  res.DESeq2.df_12vs13$padj < alpha]
# head(induced.sign)
# names(term.induced)

term.induced_12vs13 <- gprofiler(query=induced.sign_12vs13, organism="scerevisiae")
term.induced_12vs13 <- term.induced_12vs13[order(term.induced_12vs13$p.value),]
# term.induced$p.value
kable(term.induced_12vs13[1:10,c("term.name",
                                 "term.size",
                                 "query.size",
                                 "overlap.size",
                                 "recall",
                                 "precision",
                                 "p.value", 
                                 "intersection")], 
      format.args=c(engeneer=TRUE, digits=3), caption="**Table: functional analysis wit gProfileR. ** ")

## now, using the list of repressed genes.

res.DESeq2.df_12vs13 <- na.omit(data.frame(res.DESeq2_12vs13))
repressed.sign_12vs13 <- rownames(res.DESeq2.df_12vs13)[res.DESeq2.df_12vs13$log2FoldChange <= -2 &  res.DESeq2.df_12vs13$padj < alpha]
head(repressed.sign_12vs13)


term.repressed_12vs13 <- gprofiler(query=repressed.sign_12vs13, organism="scerevisiae")
term.repressed_12vs13 <- term.repressed_12vs13[order(term.repressed_12vs13$p.value),]
kable(head(term.induced_12vs13[,c("p.value", "term.name","intersection")], 10))



#design(dds_12vs13) <- ~ Preservation + Site

dds_12vs13 <- DESeq(dds_12vs13)

res_12vs13 <-results(dds_12vs13)


summary(res_12vs13)


resBigFC_12vs13 <- results(dds_12vs13, lfcThreshold=1, altHypothesis="greaterAbs")
resSort_12vs13 <- res_12vs13[order(res_12vs13$pvalue),]

head(resSort_12vs13)


## visualizing results of the DF genes

# Expression Heatmap

normlzd_dds_12vs13 <- counts(dds_12vs13, normalized=T)

res_sig_12vs13 <- data.frame(normlzd_dds_12vs13[resSort_12vs13$pvalue, ])

# write.csv(res_sig_12vs13,file="res_sig_12vs13")

log_dds_12vs13<-rlog(dds_12vs13)

# Densities of log2(counts). Each curve corresponds to one sample.


##### Scatter plots
epsilon <- 1 # pseudo-count to avoid problems with log(0)


nb.pairs <-6

## Define a function to draw a scatter plot for a pair of variables (samples) with density colors
plotFun <- function(x,y){ 
  dns <- densCols(x,y); 
  points(x,y, col=dns, pch=".", panel.first=grid());  
  #  abline(a=0, b=1, col="brown")
}

## Plot the scatter plot for a few pairs of variables selected at random
set.seed(123) # forces the random number generator to produce fixed results. Should generally not be used, except for the sake of demonstration with a particular selection. 
png("scatter_12vs13.png", width = 8, height = 8, units = 'in', res = 300)

pairs(log2(dat1_12vs13[,sample(ncol(dat1_12vs13), nb.pairs)] + epsilon), 
      panel=plotFun, lower.panel = NULL)
dev.off()

































# Top Differential Expressed Genes sorted by p-value

resSort_12vs13 <- res.DESeq2_12vs13[order(res.DESeq2_12vs13$pvalue),]
head(resSort_12vs13)
sum(res.DESeq2_12vs13$pvalue < 0.05, na.rm=TRUE)
sum(!is.na(res.DESeq2_12vs13$pvalue))
resSig_12vs13 <- subset(res.DESeq2_12vs13, padj < 0.5)

# No significant genes significantly different between Freeze x RNA Later

head(resSig_12vs13[ order(resSig_12vs13$log2FoldChange), ])
head(resSig_12vs13[ order(resSig_12vs13$log2FoldChange, decreasing = TRUE), ])

# visualizing results of the DF genes with plots  ------------------------------- No DF GENES THEREFORE, NO PLOTS

topGene_12vs13 <- rownames(res.DESeq2_12vs13)[which.min(res.DESeq2_12vs13$padj)]

geneCounts_12vs13 <- plotCounts(dds_12vs13, gene = topGene_12vs13, intgroup = c("Site"),
                                returnData = TRUE)
png("DF_genecounts_12vs13.png", width = 8, height = 4, units = 'in', res = 300)
ggplot(geneCounts_12vs13, aes(x = "site", y = count, color = "site")) +
  scale_y_log10() +  geom_beeswarm(cex = 3)
dev.off()




#### if the site doesnt work in the plots change for the condition



#### no results because there is no DF genes

png("DF_12vs13.png", width = 8, height = 4, units = 'in', res = 300)
ggplot(geneCounts_12vs13, aes(x = "site", y = count, color = "site")) +
  scale_y_log10() +  geom_beeswarm(cex = 3)
dev.off()

png("DF2_12vs13.png", width = 8, height = 4, units = 'in', res = 300)
ggplot(geneCounts_12vs13, aes(x = metadata_12vs13$Site, y = count, color = metadata_12vs13$Site, group = metadata_12vs13$Site)) +
  scale_y_log10() + geom_point(size = 3) + geom_line()
dev.off()
######### gene clustering

# 20 genes with the highest variance across samples - Not necessarily statistical significant expressed

### regularized-logarithm transformation or rlog

# regularized-logarithm transformation
rlog_12vs13 <- rlog(dds_12vs13,blind = F)
# extract the rlog matris from the object
rlog_12vs13_mat <- assay(rlog_12vs13)
# compute pairwise correlation values
rlog_12vs13_cor <- cor(rlog_12vs13_mat)
# rlog_cor


topVarGenes_12vs13 <- head(order(rowVars(assay(rlog_12vs13)), decreasing = TRUE), 20)

mat_12vs13  <- assay(rlog_12vs13)[ topVarGenes_12vs13, ]
mat_12vs13  <- mat_12vs13 - rowMeans(mat_12vs13)
anno_12vs13 <- as.data.frame(colData(rlog_12vs13)[, c("Site","Preservation")])

png("top_20_genes_12vs13.png", width = 8, height = 4, units = 'in', res = 300)
pheatmap(mat_12vs13, annotation_col = anno_12vs13)
dev.off()
