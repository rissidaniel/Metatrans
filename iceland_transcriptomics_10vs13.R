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


###### ##################################### sample comparison_10vs13



head(dat1)
# "Cultive", "10", "10", "10", "11", "11", "11", "12", "12", "12", "13", "13", "13", "14", "14", "14", "Neg", "Neg", "Neg


### subsetting count data and metadata for the sample sites 
dat1_10vs13 <- dat1 %>% select(starts_with(c("X10","X13")))

metadata_10vs13 <- metadata[metadata$Site %in% c(10, 13), ]

## Deseq Matrix
dds_10vs13 <- DESeqDataSetFromMatrix(countData= dat1_10vs13,
                                     colData = metadata_10vs13,
                                     design = ~ 1)

# removing NA rows

keep_10vs13 <- rowSums(counts(dds_10vs13)) > 1
dds_10vs13 <- dds_10vs13[keep_10vs13,]
nrow(dds_10vs13)

dds_10vs13 <- estimateSizeFactors(dds_10vs13)



# Count variance is related to mean

## Computing mean and variance
norm.counts_10vs13 <- counts(dds_10vs13, normalized=TRUE)
mean.counts_10vs13 <- rowMeans(norm.counts_10vs13)
variance.counts_10vs13 <- apply(norm.counts_10vs13, 1, var)

## sum(mean.counts==0) # Number of completely undetected genes

norm.counts.stats_10vs13 <- data.frame(
  min=apply(norm.counts_10vs13, 2, min),
  mean=apply(norm.counts_10vs13, 2, mean),
  median=apply(norm.counts_10vs13, 2, median),
  max=apply(norm.counts_10vs13, 2, max),
  zeros=apply(norm.counts_10vs13==0, 2, sum),
  percent.zeros=100*apply(norm.counts_10vs13==0, 2, sum)/nrow(norm.counts_10vs13),
  perc05=apply(norm.counts_10vs13, 2, quantile, 0.05),
  perc10=apply(norm.counts_10vs13, 2, quantile, 0.10),
  perc90=apply(norm.counts_10vs13, 2, quantile, 0.90),
  perc95=apply(norm.counts_10vs13, 2, quantile, 0.95)
)

kable(norm.counts.stats_10vs13)


## Mean and variance relationship
mean.var.col_10vs13 <- densCols(x=log2(mean.counts_10vs13), y=log2(variance.counts_10vs13))

png("Mean_and variance_relationship_10vs13.png", width = 8, height = 4, units = 'in', res = 300)

plot(x=log2(mean.counts_10vs13), y=log2(variance.counts_10vs13), pch=16, cex=0.5, 
     col=mean.var.col_10vs13, main="Mean-variance relationship",
     xlab="Mean log2(normalized counts) per gene",
     ylab="Variance of log2(normalized counts)",
     panel.first = grid())
abline(a=0, b=1, col="brown")

dev.off()
## Modelling read counts through a negative binomial

## Performing estimation of dispersion parameter
dds.disp_10vs13 <- estimateDispersions(dds_10vs13)

## A diagnostic plot which
## shows the mean of normalized counts (x axis)
## and dispersion estimate for each genes
png("Dispersion_10vs13.png", width = 8, height = 4, units = 'in', res = 300)
plotDispEsts(dds.disp_10vs13)
dev.off()

##### Performing differential expression call

alpha <- 0.0001
wald.test_10vs13 <- nbinomWaldTest(dds.disp_10vs13)
res.DESeq2_10vs13 <- results(wald.test_10vs13, alpha=alpha, pAdjustMethod="BH")

res.DESeq2_10vs13 <- res.DESeq2_10vs13[order(res.DESeq2_10vs13$padj),]

## Draw an histogram of the p-values
png("histogram_10vs13.png", width = 8, height = 4, units = 'in', res = 300)
hist(res.DESeq2_10vs13$padj, breaks=20, col="grey", main="DESeq2 p-value distribution", xlab="DESeq2 P-value", ylab="Number of genes")
dev.off()
###########Volcano plot

alpha <- 0.01 # Threshold on the adjusted p-value
cols_10vs13 <- densCols(res.DESeq2_10vs13$log2FoldChange, -log10(res.DESeq2_10vs13$pvalue))

png("volcano_10vs13.png", width = 8, height = 4, units = 'in', res = 300)
plot(res.DESeq2_10vs13$log2FoldChange, -log10(res.DESeq2_10vs13$padj), col=cols_10vs13, panel.first=grid(),
     main="Volcano plot", xlab="Effect size: log2(fold-change)", ylab="-log10(adjusted p-value)",
     pch=20, cex=0.6)
abline(v=0)
abline(v=c(-1,1), col="brown")
abline(h=-log10(alpha), col="brown")
gn.selected_10vs13 <- abs(res.DESeq2_10vs13$log2FoldChange) > 2 & res.DESeq2_10vs13$padj < alpha 
text(res.DESeq2_10vs13$log2FoldChange[gn.selected_10vs13],
     -log10(res.DESeq2_10vs13$padj)[gn.selected_10vs13],
     lab=rownames(res.DESeq2_10vs13)[gn.selected_10vs13 ], cex=0.4)
dev.off()

### Check the expression levels of the most differentially expressed gene

gn.most.sign_10vs13 <- rownames(res.DESeq2_10vs13)[1]
gn.most.diff.val_10vs13 <- counts(dds_10vs13, normalized=T)[gn.most.sign_10vs13,]

png("DFgenes_10vs13.png", width = 8, height = 4, units = 'in', res = 300)
barplot(gn.most.diff.val_10vs13, col=metadata_10vs13, main=gn.most.sign_10vs13, las=2, cex.names=0.5)
dev.off()

### MA plot

## Draw a MA plot.
## Genes with adjusted p-values below 1% are shown
png("MA_10vs13.png", width = 8, height = 8, units = 'in', res = 300)

plotMA(res.DESeq2_10vs13, colNonSig = "blue")
abline(h=c(-1:1), col="red")
dev.off()


# functional enrichment analysis using the list of induced genes. 
# This step will be performed using the gProfileR R library.


epsilon <- 1 # pseudo-count to avoid problems with log(0)

res.DESeq2.df_10vs13 <- na.omit(data.frame(res.DESeq2_10vs13))
induced.sign_10vs13 <- rownames(res.DESeq2.df_10vs13)[res.DESeq2.df_10vs13$log2FoldChange >= 2 &  res.DESeq2.df_10vs13$padj < alpha]
# head(induced.sign)
# names(term.induced)

term.induced_10vs13 <- gprofiler(query=induced.sign_10vs13, organism="scerevisiae")
term.induced_10vs13 <- term.induced_10vs13[order(term.induced_10vs13$p.value),]
# term.induced$p.value
kable(term.induced_10vs13[1:10,c("term.name",
                                 "term.size",
                                 "query.size",
                                 "overlap.size",
                                 "recall",
                                 "precision",
                                 "p.value", 
                                 "intersection")], 
      format.args=c(engeneer=TRUE, digits=3), caption="**Table: functional analysis wit gProfileR. ** ")

## now, using the list of repressed genes.

res.DESeq2.df_10vs13 <- na.omit(data.frame(res.DESeq2_10vs13))
repressed.sign_10vs13 <- rownames(res.DESeq2.df_10vs13)[res.DESeq2.df_10vs13$log2FoldChange <= -2 &  res.DESeq2.df_10vs13$padj < alpha]
head(repressed.sign_10vs13)


term.repressed_10vs13 <- gprofiler(query=repressed.sign_10vs13, organism="scerevisiae")
term.repressed_10vs13 <- term.repressed_10vs13[order(term.repressed_10vs13$p.value),]
kable(head(term.induced_10vs13[,c("p.value", "term.name","intersection")], 10))



#design(dds_10vs13) <- ~ Preservation + Site

dds_10vs13 <- DESeq(dds_10vs13)

res_10vs13 <-results(dds_10vs13)


summary(res_10vs13)


resBigFC_10vs13 <- results(dds_10vs13, lfcThreshold=1, altHypothesis="greaterAbs")
resSort_10vs13 <- res_10vs13[order(res_10vs13$pvalue),]

head(resSort_10vs13)


## visualizing results of the DF genes

# Expression Heatmap

normlzd_dds_10vs13 <- counts(dds_10vs13, normalized=T)

res_sig_10vs13 <- data.frame(normlzd_dds_10vs13[resSort_10vs13$pvalue, ])

# write.csv(res_sig_10vs13,file="res_sig_10vs13")

log_dds_10vs13<-rlog(dds_10vs13)

# Densities of log2(counts). Each curve corresponds to one sample.


##### Scatter plots

nb.pairs <-6

## Define a function to draw a scatter plot for a pair of variables (samples) with density colors
plotFun <- function(x,y){ 
  dns <- densCols(x,y); 
  points(x,y, col=dns, pch=".", panel.first=grid());  
  #  abline(a=0, b=1, col="brown")
}

## Plot the scatter plot for a few pairs of variables selected at random
set.seed(123) # forces the random number generator to produce fixed results. Should generally not be used, except for the sake of demonstration with a particular selection. 
png("scatter_10vs13.png", width = 8, height = 8, units = 'in', res = 300)

pairs(log2(dat1_10vs13[,sample(ncol(dat1_10vs13), nb.pairs)] + epsilon), 
      panel=plotFun, lower.panel = NULL)
dev.off()

































# Top Differential Expressed Genes sorted by p-value

resSort_10vs13 <- res.DESeq2_10vs13[order(res.DESeq2_10vs13$pvalue),]
head(resSort_10vs13)
sum(res.DESeq2_10vs13$pvalue < 0.05, na.rm=TRUE)
sum(!is.na(res.DESeq2_10vs13$pvalue))
resSig_10vs13 <- subset(res.DESeq2_10vs13, padj < 0.5)

# No significant genes significantly different between Freeze x RNA Later

head(resSig_10vs13[ order(resSig_10vs13$log2FoldChange), ])
head(resSig_10vs13[ order(resSig_10vs13$log2FoldChange, decreasing = TRUE), ])

# visualizing results of the DF genes with plots  ------------------------------- No DF GENES THEREFORE, NO PLOTS

topGene_10vs13 <- rownames(res.DESeq2_10vs13)[which.min(res.DESeq2_10vs13$padj)]

geneCounts_10vs13 <- plotCounts(dds_10vs13, gene = topGene_10vs13, intgroup = c("Site"),
                                returnData = TRUE)
png("DF_genecounts_10vs13.png", width = 8, height = 4, units = 'in', res = 300)
ggplot(geneCounts_10vs13, aes(x = "site", y = count, color = "site")) +
  scale_y_log10() +  geom_beeswarm(cex = 3)
dev.off()




#### if the site doesnt work in the plots change for the condition



#### no results because there is no DF genes

png("DF_RLater_vs_Freeze_gene_ID.png", width = 8, height = 4, units = 'in', res = 300)
ggplot(geneCounts_10vs13, aes(x = "site", y = count, color = "site")) +
  scale_y_log10() +  geom_beeswarm(cex = 3)
dev.off()

png("DF2_RLater_vs_Freeze_gene_ID.png", width = 8, height = 4, units = 'in', res = 300)
ggplot(geneCounts_10vs13, aes(x = metadata_10vs13$Site, y = count, color = metadata_10vs13$Site, group = metadata_10vs13$Site)) +
  scale_y_log10() + geom_point(size = 3) + geom_line()
dev.off()
######### gene clustering

# 20 genes with the highest variance across samples - Not necessarily statistical significant expressed

### regularized-logarithm transformation or rlog

# regularized-logarithm transformation
rlog_10vs13 <- rlog(dds_10vs13,blind = F)
# extract the rlog matris from the object
rlog_10vs13_mat <- assay(rlog_10vs13)
# compute pairwise correlation values
rlog_10vs13_cor <- cor(rlog_10vs13_mat)
# rlog_cor


topVarGenes_10vs13 <- head(order(rowVars(assay(rlog_10vs13)), decreasing = TRUE), 20)

mat_10vs13  <- assay(rlog_10vs13)[ topVarGenes_10vs13, ]
mat_10vs13  <- mat_10vs13 - rowMeans(mat_10vs13)
anno_10vs13 <- as.data.frame(colData(rlog_10vs13)[, c("Site","Preservation")])

png("top_20_genes_10vs13.png", width = 8, height = 4, units = 'in', res = 300)
pheatmap(mat_10vs13, annotation_col = anno_10vs13)
dev.off()
