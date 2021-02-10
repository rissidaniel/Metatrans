library(DESeq2)
library(RColorBrewer)
library(pheatmap)
library(tidyverse)
library(gplots)
library(plyr)
library(affy)
library(glmpca)
library(knitr)
library(gProfileR)
library(knitr)

# setwd("preservation/")
data <- read.csv("../iceland_genes_new.tsv", sep = "\t", header = T)
metadata <-  read.csv("../metadata2.txt", sep = "\t")

dat1 <- ddply(data,"Gene.name",numcolwise(sum))


rownames(dat1) <- dat1[,1]
dat1 <- dat1[,-1]
# "10", "10", "10", "11", "11", "11", "12", "12", "12", "13", "13", "13", "14", "14", "14", "Neg", "Neg", "Neg

# condition =  factor(c("Cultive", "10", "10", "10", "11", "11", "11", "12", "12", "12", "13", "13", "13", "14", "14", "14", "Neg", "Neg", "Neg"))
# coldata <- data.frame(row.names=colnames(dat1), condition)


rownames(metadata)
colnames(dat1)

all(rownames(metadata) == colnames(dat1))


###### ##################################### sample comparison_RLater_vs_Freeze



head(dat1)
# "Cultive", "10", "10", "10", "11", "11", "11", "12", "12", "12", "13", "13", "13", "14", "14", "14", "Neg", "Neg", "Neg


### subsetting count data and metadata for the sample sites 
# dat1_Freeze <- dat1 %>% select(starts_with(c("F","L")))
dat1_RLater_vs_Freeze <- dat1[, -c(1,4,7,10,13,16,17,18,19)]

metadata_RLater_vs_Freeze <- metadata[metadata$Preservation %in% c("Freeze", "Rlater"), ]

## Deseq Matrix
dds_RLater_vs_Freeze <- DESeqDataSetFromMatrix(countData= dat1_RLater_vs_Freeze,
                                     colData = metadata_RLater_vs_Freeze,
                                     design = ~ 1)

# removing NA rows

keep_RLater_vs_Freeze <- rowSums(counts(dds_RLater_vs_Freeze)) > 1
dds_RLater_vs_Freeze <- dds_RLater_vs_Freeze[keep_RLater_vs_Freeze,]
nrow(dds_RLater_vs_Freeze)

dds_RLater_vs_Freeze <- estimateSizeFactors(dds_RLater_vs_Freeze)



# Count variance is related to mean

## Computing mean and variance
norm.counts_RLater_vs_Freeze <- counts(dds_RLater_vs_Freeze, normalized=TRUE)
mean.counts_RLater_vs_Freeze <- rowMeans(norm.counts_RLater_vs_Freeze)
variance.counts_RLater_vs_Freeze <- apply(norm.counts_RLater_vs_Freeze, 1, var)

## sum(mean.counts==0) # Number of completely undetected genes

norm.counts.stats_RLater_vs_Freeze <- data.frame(
  min=apply(norm.counts_RLater_vs_Freeze, 2, min),
  mean=apply(norm.counts_RLater_vs_Freeze, 2, mean),
  median=apply(norm.counts_RLater_vs_Freeze, 2, median),
  max=apply(norm.counts_RLater_vs_Freeze, 2, max),
  zeros=apply(norm.counts_RLater_vs_Freeze==0, 2, sum),
  percent.zeros=100*apply(norm.counts_RLater_vs_Freeze==0, 2, sum)/nrow(norm.counts_RLater_vs_Freeze),
  perc05=apply(norm.counts_RLater_vs_Freeze, 2, quantile, 0.05),
  perc10=apply(norm.counts_RLater_vs_Freeze, 2, quantile, 0.10),
  perc90=apply(norm.counts_RLater_vs_Freeze, 2, quantile, 0.90),
  perc95=apply(norm.counts_RLater_vs_Freeze, 2, quantile, 0.95)
)

kable(norm.counts.stats_RLater_vs_Freeze)


## Mean and variance relationship
mean.var.col_RLater_vs_Freeze <- densCols(x=log2(mean.counts_RLater_vs_Freeze), y=log2(variance.counts_RLater_vs_Freeze))

png("Mean_and variance_relationship_RLater_vs_Freeze.png", width = 8, height = 4, units = 'in', res = 300)

plot(x=log2(mean.counts_RLater_vs_Freeze), y=log2(variance.counts_RLater_vs_Freeze), pch=16, cex=0.5, 
     col=mean.var.col_RLater_vs_Freeze, main="Mean-variance relationship",
     xlab="Mean log2(normalized counts) per gene",
     ylab="Variance of log2(normalized counts)",
     panel.first = grid())
abline(a=0, b=1, col="brown")

dev.off()
## Modelling read counts through a negative binomial

## Performing estimation of dispersion parameter
dds.disp_RLater_vs_Freeze <- estimateDispersions(dds_RLater_vs_Freeze)

## A diagnostic plot which
## shows the mean of normalized counts (x axis)
## and dispersion estimate for each genes
png("Dispersion_RLater_vs_Freeze.png", width = 8, height = 4, units = 'in', res = 300)
plotDispEsts(dds.disp_RLater_vs_Freeze)
dev.off()

##### Performing differential expression call

alpha <- 0.0001
wald.test_RLater_vs_Freeze <- nbinomWaldTest(dds.disp_RLater_vs_Freeze)
res.DESeq2_RLater_vs_Freeze <- results(wald.test_RLater_vs_Freeze, alpha=alpha, pAdjustMethod="BH")

res.DESeq2_RLater_vs_Freeze <- res.DESeq2_RLater_vs_Freeze[order(res.DESeq2_RLater_vs_Freeze$padj),]

## Draw an histogram of the p-values
png("histogram_RLater_vs_Freeze.png", width = 8, height = 4, units = 'in', res = 300)
hist(res.DESeq2_RLater_vs_Freeze$padj, breaks=20, col="grey", main="DESeq2 p-value distribution", xlab="DESeq2 P-value", ylab="Number of genes")
dev.off()
###########Volcano plot

alpha <- 0.01 # Threshold on the adjusted p-value
cols_RLater_vs_Freeze <- densCols(res.DESeq2_RLater_vs_Freeze$log2FoldChange, -log10(res.DESeq2_RLater_vs_Freeze$pvalue))

png("volcano_RLater_vs_Freeze.png", width = 8, height = 4, units = 'in', res = 300)
plot(res.DESeq2_RLater_vs_Freeze$log2FoldChange, -log10(res.DESeq2_RLater_vs_Freeze$padj), col=cols_RLater_vs_Freeze, panel.first=grid(),
     main="Volcano plot", xlab="Effect size: log2(fold-change)", ylab="-log10(adjusted p-value)",
     pch=20, cex=0.6)
abline(v=0)
abline(v=c(-1,1), col="brown")
abline(h=-log10(alpha), col="brown")
gn.selected_RLater_vs_Freeze <- abs(res.DESeq2_RLater_vs_Freeze$log2FoldChange) > 2 & res.DESeq2_RLater_vs_Freeze$padj < alpha 
text(res.DESeq2_RLater_vs_Freeze$log2FoldChange[gn.selected_RLater_vs_Freeze],
     -log10(res.DESeq2_RLater_vs_Freeze$padj)[gn.selected_RLater_vs_Freeze],
     lab=rownames(res.DESeq2_RLater_vs_Freeze)[gn.selected_RLater_vs_Freeze ], cex=0.4)
dev.off()

### Check the expression levels of the most differentially expressed gene - no significant different expressed genes

gn.most.sign_RLater_vs_Freeze <- rownames(res.DESeq2_RLater_vs_Freeze)[1]
gn.most.diff.val_RLater_vs_Freeze <- counts(dds_RLater_vs_Freeze, normalized=T)[gn.most.sign_RLater_vs_Freeze,]

png("DFgenes_RLater_vs_Freeze.png", width = 8, height = 4, units = 'in', res = 300)
barplot(gn.most.diff.val_RLater_vs_Freeze, col=metadata_RLater_vs_Freeze, main=gn.most.sign_RLater_vs_Freeze, las=2, cex.names=0.5)
dev.off()

### MA plot

## Draw a MA plot.
## Genes with adjusted p-values below 1% are shown
png("MA_RLater_vs_Freeze.png", width = 8, height = 8, units = 'in', res = 300)

plotMA(res.DESeq2_RLater_vs_Freeze, colNonSig = "blue")
abline(h=c(-1:1), col="red")
dev.off()


# functional enrichment analysis using the list of induced genes. 
# This step will be performed using the gProfileR R library.


epsilon <- 1 # pseudo-count to avoid problems with log(0)

res.DESeq2.df_RLater_vs_Freeze <- na.omit(data.frame(res.DESeq2_RLater_vs_Freeze))
induced.sign_RLater_vs_Freeze <- rownames(res.DESeq2.df_RLater_vs_Freeze)[res.DESeq2.df_RLater_vs_Freeze$log2FoldChange >= 2 &  res.DESeq2.df_RLater_vs_Freeze$padj < alpha]
# head(induced.sign)
# names(term.induced)

term.induced_RLater_vs_Freeze <- gprofiler(query=induced.sign_RLater_vs_Freeze, organism="scerevisiae")
term.induced_RLater_vs_Freeze <- term.induced_RLater_vs_Freeze[order(term.induced_RLater_vs_Freeze$p.value),]
# term.induced$p.value
kable(term.induced_RLater_vs_Freeze[1:10,c("term.name",
                                 "term.size",
                                 "query.size",
                                 "overlap.size",
                                 "recall",
                                 "precision",
                                 "p.value", 
                                 "intersection")], 
      format.args=c(engeneer=TRUE, digits=3), caption="**Table: functional analysis wit gProfileR. ** ")

## now, using the list of repressed genes.

res.DESeq2.df_RLater_vs_Freeze <- na.omit(data.frame(res.DESeq2_RLater_vs_Freeze))
repressed.sign_RLater_vs_Freeze <- rownames(res.DESeq2.df_RLater_vs_Freeze)[res.DESeq2.df_RLater_vs_Freeze$log2FoldChange <= -2 &  res.DESeq2.df_RLater_vs_Freeze$padj < alpha]
head(repressed.sign_RLater_vs_Freeze)


term.repressed_RLater_vs_Freeze <- gprofiler(query=repressed.sign_RLater_vs_Freeze, organism="scerevisiae")
term.repressed_RLater_vs_Freeze <- term.repressed_RLater_vs_Freeze[order(term.repressed_RLater_vs_Freeze$p.value),]
kable(head(term.induced_RLater_vs_Freeze[,c("p.value", "term.name","intersection")], 10))



#design(dds_RLater_vs_Freeze) <- ~ Preservation + Site

dds_RLater_vs_Freeze <- DESeq(dds_RLater_vs_Freeze)

res_RLater_vs_Freeze <-results(dds_RLater_vs_Freeze)


summary(res_RLater_vs_Freeze)


resBigFC_RLater_vs_Freeze <- results(dds_RLater_vs_Freeze, lfcThreshold=1, altHypothesis="greaterAbs")
resSort_RLater_vs_Freeze <- res_RLater_vs_Freeze[order(res_RLater_vs_Freeze$pvalue),]

head(resSort_RLater_vs_Freeze)


## visualizing results of the DF genes

# Expression Heatmap

normlzd_dds_RLater_vs_Freeze <- counts(dds_RLater_vs_Freeze, normalized=T)

res_sig_RLater_vs_Freeze <- data.frame(normlzd_dds_RLater_vs_Freeze[resSort_RLater_vs_Freeze$pvalue, ])

# write.csv(res_sig_RLater_vs_Freeze,file="res_sig_RLater_vs_Freeze")

log_dds_RLater_vs_Freeze<-rlog(dds_RLater_vs_Freeze)

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
png("scatter_RLater_vs_Freeze.png", width = 8, height = 8, units = 'in', res = 300)

pairs(log2(dat1_RLater_vs_Freeze[,sample(ncol(dat1_RLater_vs_Freeze), nb.pairs)] + epsilon), 
      panel=plotFun, lower.panel = NULL)
dev.off()

































# Top Differential Expressed Genes sorted by p-value

resSort_RLater_vs_Freeze <- res.DESeq2_RLater_vs_Freeze[order(res.DESeq2_RLater_vs_Freeze$pvalue),]
head(resSort_RLater_vs_Freeze)
sum(res.DESeq2_RLater_vs_Freeze$pvalue < 0.05, na.rm=TRUE)
sum(!is.na(res.DESeq2_RLater_vs_Freeze$pvalue))
resSig_RLater_vs_Freeze <- subset(res.DESeq2_RLater_vs_Freeze, padj < 0.5)

# No significant genes significantly different between Freeze x RNA Later

head(resSig_RLater_vs_Freeze[ order(resSig_RLater_vs_Freeze$log2FoldChange), ])
head(resSig_RLater_vs_Freeze[ order(resSig_RLater_vs_Freeze$log2FoldChange, decreasing = TRUE), ])

# visualizing results of the DF genes with plots  ------------------------------- No DF GENES THEREFORE, NO PLOTS

topGene_RLater_vs_Freeze <- rownames(res.DESeq2_RLater_vs_Freeze)[which.min(res.DESeq2_RLater_vs_Freeze$padj)]

geneCounts_RLater_vs_Freeze <- plotCounts(dds_RLater_vs_Freeze, gene = topGene_RLater_vs_Freeze, intgroup = c("Site"),
                                returnData = TRUE)
png("DF_genecounts_RLater_vs_Freeze.png", width = 8, height = 4, units = 'in', res = 300)
ggplot(geneCounts_RLater_vs_Freeze, aes(x = "site", y = count, color = "site")) +
  scale_y_log10() +  geom_beeswarm(cex = 3)
dev.off()




#### if the site doesnt work in the plots change for the condition



#### no results because there is no DF genes

png("DF_RLater_vs_Freeze_gene_ID.png", width = 8, height = 4, units = 'in', res = 300)
ggplot(geneCounts_RLater_vs_Freeze, aes(x = "site", y = count, color = "site")) +
  scale_y_log10() +  geom_beeswarm(cex = 3)
dev.off()

png("DF2_RLater_vs_Freeze_gene_ID.png", width = 8, height = 4, units = 'in', res = 300)
ggplot(geneCounts_RLater_vs_Freeze, aes(x = metadata_RLater_vs_Freeze$Site, y = count, color = metadata_RLater_vs_Freeze$Site, group = metadata_RLater_vs_Freeze$Site)) +
  scale_y_log10() + geom_point(size = 3) + geom_line()
dev.off()
######### gene clustering

# 20 genes with the highest variance across samples - Not necessarily statistical significant expressed

### regularized-logarithm transformation or rlog

# regularized-logarithm transformation
rlog_RLater_vs_Freeze <- rlog(dds_RLater_vs_Freeze,blind = F)
# extract the rlog matris from the object
rlog_RLater_vs_Freeze_mat <- assay(rlog_RLater_vs_Freeze)
# compute pairwise correlation values
rlog_RLater_vs_Freeze_cor <- cor(rlog_RLater_vs_Freeze_mat)
# rlog_cor


topVarGenes_RLater_vs_Freeze <- head(order(rowVars(assay(rlog_RLater_vs_Freeze)), decreasing = TRUE), 20)

mat_RLater_vs_Freeze  <- assay(rlog_RLater_vs_Freeze)[ topVarGenes_RLater_vs_Freeze, ]
mat_RLater_vs_Freeze  <- mat_RLater_vs_Freeze - rowMeans(mat_RLater_vs_Freeze)
anno_RLater_vs_Freeze <- as.data.frame(colData(rlog_RLater_vs_Freeze)[, c("Site","Preservation")])

png("top_20_genes_RLater_vs_Freeze.png", width = 8, height = 4, units = 'in', res = 300)
pheatmap(mat_RLater_vs_Freeze, annotation_col = anno_RLater_vs_Freeze)
dev.off()
