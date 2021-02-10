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


###### ##################################### sample comparison_11vs12



head(dat1)
# "Cultive", "10", "10", "10", "11", "11", "11", "12", "12", "12", "13", "13", "13", "14", "14", "14", "Neg", "Neg", "Neg


### subsetting count data and metadata for the sample sites 
dat1_11vs12 <- dat1 %>% select(starts_with(c("X11","X12")))

metadata_11vs12 <- metadata[metadata$Site %in% c(11, 12), ]

## Deseq Matrix
dds_11vs12 <- DESeqDataSetFromMatrix(countData= dat1_11vs12,
                                     colData = metadata_11vs12,
                                     design = ~ 1)

# removing NA rows

keep_11vs12 <- rowSums(counts(dds_11vs12)) > 1
dds_11vs12 <- dds_11vs12[keep_11vs12,]
nrow(dds_11vs12)

dds_11vs12 <- estimateSizeFactors(dds_11vs12)



# Count variance is related to mean

## Computing mean and variance
norm.counts_11vs12 <- counts(dds_11vs12, normalized=TRUE)
mean.counts_11vs12 <- rowMeans(norm.counts_11vs12)
variance.counts_11vs12 <- apply(norm.counts_11vs12, 1, var)

## sum(mean.counts==0) # Number of completely undetected genes

norm.counts.stats_11vs12 <- data.frame(
  min=apply(norm.counts_11vs12, 2, min),
  mean=apply(norm.counts_11vs12, 2, mean),
  median=apply(norm.counts_11vs12, 2, median),
  max=apply(norm.counts_11vs12, 2, max),
  zeros=apply(norm.counts_11vs12==0, 2, sum),
  percent.zeros=100*apply(norm.counts_11vs12==0, 2, sum)/nrow(norm.counts_11vs12),
  perc05=apply(norm.counts_11vs12, 2, quantile, 0.05),
  perc10=apply(norm.counts_11vs12, 2, quantile, 0.10),
  perc90=apply(norm.counts_11vs12, 2, quantile, 0.90),
  perc95=apply(norm.counts_11vs12, 2, quantile, 0.95)
)

kable(norm.counts.stats_11vs12)


## Mean and variance relationship
mean.var.col_11vs12 <- densCols(x=log2(mean.counts_11vs12), y=log2(variance.counts_11vs12))

png("Mean_and variance_relationship_11vs12.png", width = 8, height = 4, units = 'in', res = 300)

plot(x=log2(mean.counts_11vs12), y=log2(variance.counts_11vs12), pch=16, cex=0.5, 
     col=mean.var.col_11vs12, main="Mean-variance relationship",
     xlab="Mean log2(normalized counts) per gene",
     ylab="Variance of log2(normalized counts)",
     panel.first = grid())
abline(a=0, b=1, col="brown")

dev.off()
## Modelling read counts through a negative binomial

## Performing estimation of dispersion parameter
dds.disp_11vs12 <- estimateDispersions(dds_11vs12)

## A diagnostic plot which
## shows the mean of normalized counts (x axis)
## and dispersion estimate for each genes
png("Dispersion_11vs12.png", width = 8, height = 4, units = 'in', res = 300)
plotDispEsts(dds.disp_11vs12)
dev.off()

##### Performing differential expression call

alpha <- 0.0001
wald.test_11vs12 <- nbinomWaldTest(dds.disp_11vs12)
res.DESeq2_11vs12 <- results(wald.test_11vs12, alpha=alpha, pAdjustMethod="BH")

res.DESeq2_11vs12 <- res.DESeq2_11vs12[order(res.DESeq2_11vs12$padj),]

## Draw an histogram of the p-values
png("histogram_11vs12.png", width = 8, height = 4, units = 'in', res = 300)
hist(res.DESeq2_11vs12$padj, breaks=20, col="grey", main="DESeq2 p-value distribution", xlab="DESeq2 P-value", ylab="Number of genes")
dev.off()
###########Volcano plot

alpha <- 0.01 # Threshold on the adjusted p-value
cols_11vs12 <- densCols(res.DESeq2_11vs12$log2FoldChange, -log10(res.DESeq2_11vs12$pvalue))

png("volcano_11vs12.png", width = 8, height = 4, units = 'in', res = 300)
plot(res.DESeq2_11vs12$log2FoldChange, -log10(res.DESeq2_11vs12$padj), col=cols_11vs12, panel.first=grid(),
     main="Volcano plot", xlab="Effect size: log2(fold-change)", ylab="-log10(adjusted p-value)",
     pch=20, cex=0.6)
abline(v=0)
abline(v=c(-1,1), col="brown")
abline(h=-log10(alpha), col="brown")
gn.selected_11vs12 <- abs(res.DESeq2_11vs12$log2FoldChange) > 2 & res.DESeq2_11vs12$padj < alpha 
text(res.DESeq2_11vs12$log2FoldChange[gn.selected_11vs12],
     -log10(res.DESeq2_11vs12$padj)[gn.selected_11vs12],
     lab=rownames(res.DESeq2_11vs12)[gn.selected_11vs12 ], cex=0.4)
dev.off()

### Check the expression levels of the most differentially expressed gene

gn.most.sign_11vs12 <- rownames(res.DESeq2_11vs12)[1]
gn.most.diff.val_11vs12 <- counts(dds_11vs12, normalized=T)[gn.most.sign_11vs12,]

png("DFgenes_11vs12.png", width = 8, height = 4, units = 'in', res = 300)
barplot(gn.most.diff.val_11vs12, col=metadata_11vs12, main=gn.most.sign_11vs12, las=2, cex.names=0.5)
dev.off()

### MA plot

## Draw a MA plot.
## Genes with adjusted p-values below 1% are shown
png("MA_11vs12.png", width = 8, height = 8, units = 'in', res = 300)

plotMA(res.DESeq2_11vs12, colNonSig = "blue")
abline(h=c(-1:1), col="red")
dev.off()


# functional enrichment analysis using the list of induced genes. 
# This step will be performed using the gProfileR R library.


epsilon <- 1 # pseudo-count to avoid problems with log(0)

res.DESeq2.df_11vs12 <- na.omit(data.frame(res.DESeq2_11vs12))
induced.sign_11vs12 <- rownames(res.DESeq2.df_11vs12)[res.DESeq2.df_11vs12$log2FoldChange >= 2 &  res.DESeq2.df_11vs12$padj < alpha]
# head(induced.sign)
# names(term.induced)

term.induced_11vs12 <- gprofiler(query=induced.sign_11vs12, organism="scerevisiae")
term.induced_11vs12 <- term.induced_11vs12[order(term.induced_11vs12$p.value),]
# term.induced$p.value
kable(term.induced_11vs12[1:10,c("term.name",
                                 "term.size",
                                 "query.size",
                                 "overlap.size",
                                 "recall",
                                 "precision",
                                 "p.value", 
                                 "intersection")], 
      format.args=c(engeneer=TRUE, digits=3), caption="**Table: functional analysis wit gProfileR. ** ")

## now, using the list of repressed genes.

res.DESeq2.df_11vs12 <- na.omit(data.frame(res.DESeq2_11vs12))
repressed.sign_11vs12 <- rownames(res.DESeq2.df_11vs12)[res.DESeq2.df_11vs12$log2FoldChange <= -2 &  res.DESeq2.df_11vs12$padj < alpha]
head(repressed.sign_11vs12)


term.repressed_11vs12 <- gprofiler(query=repressed.sign_11vs12, organism="scerevisiae")
term.repressed_11vs12 <- term.repressed_11vs12[order(term.repressed_11vs12$p.value),]
kable(head(term.induced_11vs12[,c("p.value", "term.name","intersection")], 10))



#design(dds_11vs12) <- ~ Preservation + Site

dds_11vs12 <- DESeq(dds_11vs12)

res_11vs12 <-results(dds_11vs12)


summary(res_11vs12)


resBigFC_11vs12 <- results(dds_11vs12, lfcThreshold=1, altHypothesis="greaterAbs")
resSort_11vs12 <- res_11vs12[order(res_11vs12$pvalue),]

head(resSort_11vs12)


## visualizing results of the DF genes

# Expression Heatmap

normlzd_dds_11vs12 <- counts(dds_11vs12, normalized=T)

res_sig_11vs12 <- data.frame(normlzd_dds_11vs12[resSort_11vs12$pvalue, ])

# write.csv(res_sig_11vs12,file="res_sig_11vs12")

log_dds_11vs12<-rlog(dds_11vs12)

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
png("scatter_11vs12.png", width = 8, height = 8, units = 'in', res = 300)

pairs(log2(dat1_11vs12[,sample(ncol(dat1_11vs12), nb.pairs)] + epsilon), 
      panel=plotFun, lower.panel = NULL)
dev.off()

































# Top Differential Expressed Genes sorted by p-value

resSort_11vs12 <- res.DESeq2_11vs12[order(res.DESeq2_11vs12$pvalue),]
head(resSort_11vs12)
sum(res.DESeq2_11vs12$pvalue < 0.05, na.rm=TRUE)
sum(!is.na(res.DESeq2_11vs12$pvalue))
resSig_11vs12 <- subset(res.DESeq2_11vs12, padj < 0.5)

# No significant genes significantly different between Freeze x RNA Later

head(resSig_11vs12[ order(resSig_11vs12$log2FoldChange), ])
head(resSig_11vs12[ order(resSig_11vs12$log2FoldChange, decreasing = TRUE), ])

# visualizing results of the DF genes with plots  ------------------------------- No DF GENES THEREFORE, NO PLOTS

topGene_11vs12 <- rownames(res.DESeq2_11vs12)[which.min(res.DESeq2_11vs12$padj)]

geneCounts_11vs12 <- plotCounts(dds_11vs12, gene = topGene_11vs12, intgroup = c("Site"),
                                returnData = TRUE)
png("DF_genecounts_11vs12.png", width = 8, height = 4, units = 'in', res = 300)
ggplot(geneCounts_11vs12, aes(x = "site", y = count, color = "site")) +
  scale_y_log10() +  geom_beeswarm(cex = 3)
dev.off()




#### if the site doesnt work in the plots change for the condition



#### no results because there is no DF genes

png("DF_11vs12.png", width = 8, height = 4, units = 'in', res = 300)
ggplot(geneCounts_11vs12, aes(x = "site", y = count, color = "site")) +
  scale_y_log10() +  geom_beeswarm(cex = 3)
dev.off()

png("DF2_11vs12.png", width = 8, height = 4, units = 'in', res = 300)
ggplot(geneCounts_11vs12, aes(x = metadata_11vs12$Site, y = count, color = metadata_11vs12$Site, group = metadata_11vs12$Site)) +
  scale_y_log10() + geom_point(size = 3) + geom_line()
dev.off()
######### gene clustering

# 20 genes with the highest variance across samples - Not necessarily statistical significant expressed

### regularized-logarithm transformation or rlog

# regularized-logarithm transformation
rlog_11vs12 <- rlog(dds_11vs12,blind = F)
# extract the rlog matris from the object
rlog_11vs12_mat <- assay(rlog_11vs12)
# compute pairwise correlation values
rlog_11vs12_cor <- cor(rlog_11vs12_mat)
# rlog_cor


topVarGenes_11vs12 <- head(order(rowVars(assay(rlog_11vs12)), decreasing = TRUE), 20)

mat_11vs12  <- assay(rlog_11vs12)[ topVarGenes_11vs12, ]
mat_11vs12  <- mat_11vs12 - rowMeans(mat_11vs12)
anno_11vs12 <- as.data.frame(colData(rlog_11vs12)[, c("Site","Preservation")])

png("top_20_genes_11vs12.png", width = 8, height = 4, units = 'in', res = 300)
pheatmap(mat_11vs12, annotation_col = anno_11vs12)
dev.off()
