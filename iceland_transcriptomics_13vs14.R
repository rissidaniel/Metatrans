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


###### ##################################### sample comparison_13vs14



head(dat1)
# "Cultive", "10", "10", "10", "11", "11", "11", "12", "12", "12", "13", "13", "13", "14", "14", "14", "Neg", "Neg", "Neg


### subsetting count data and metadata for the sample sites 
dat1_13vs14 <- dat1 %>% select(starts_with(c("X13","X14")))

metadata_13vs14 <- metadata[metadata$Site %in% c(13, 14), ]

## Deseq Matrix
dds_13vs14 <- DESeqDataSetFromMatrix(countData= dat1_13vs14,
                                     colData = metadata_13vs14,
                                     design = ~ 1)

# removing NA rows

keep_13vs14 <- rowSums(counts(dds_13vs14)) > 1
dds_13vs14 <- dds_13vs14[keep_13vs14,]
nrow(dds_13vs14)

dds_13vs14 <- estimateSizeFactors(dds_13vs14)



# Count variance is related to mean

## Computing mean and variance
norm.counts_13vs14 <- counts(dds_13vs14, normalized=TRUE)
mean.counts_13vs14 <- rowMeans(norm.counts_13vs14)
variance.counts_13vs14 <- apply(norm.counts_13vs14, 1, var)

## sum(mean.counts==0) # Number of completely undetected genes

norm.counts.stats_13vs14 <- data.frame(
  min=apply(norm.counts_13vs14, 2, min),
  mean=apply(norm.counts_13vs14, 2, mean),
  median=apply(norm.counts_13vs14, 2, median),
  max=apply(norm.counts_13vs14, 2, max),
  zeros=apply(norm.counts_13vs14==0, 2, sum),
  percent.zeros=100*apply(norm.counts_13vs14==0, 2, sum)/nrow(norm.counts_13vs14),
  perc05=apply(norm.counts_13vs14, 2, quantile, 0.05),
  perc10=apply(norm.counts_13vs14, 2, quantile, 0.10),
  perc90=apply(norm.counts_13vs14, 2, quantile, 0.90),
  perc95=apply(norm.counts_13vs14, 2, quantile, 0.95)
)

kable(norm.counts.stats_13vs14)


## Mean and variance relationship
mean.var.col_13vs14 <- densCols(x=log2(mean.counts_13vs14), y=log2(variance.counts_13vs14))

png("Mean_and variance_relationship_13vs14.png", width = 8, height = 4, units = 'in', res = 300)

plot(x=log2(mean.counts_13vs14), y=log2(variance.counts_13vs14), pch=16, cex=0.5, 
     col=mean.var.col_13vs14, main="Mean-variance relationship",
     xlab="Mean log2(normalized counts) per gene",
     ylab="Variance of log2(normalized counts)",
     panel.first = grid())
abline(a=0, b=1, col="brown")

dev.off()
## Modelling read counts through a negative binomial

## Performing estimation of dispersion parameter
dds.disp_13vs14 <- estimateDispersions(dds_13vs14)

## A diagnostic plot which
## shows the mean of normalized counts (x axis)
## and dispersion estimate for each genes
png("Dispersion_13vs14.png", width = 8, height = 4, units = 'in', res = 300)
plotDispEsts(dds.disp_13vs14)
dev.off()

##### Performing differential expression call

alpha <- 0.0001
wald.test_13vs14 <- nbinomWaldTest(dds.disp_13vs14)
res.DESeq2_13vs14 <- results(wald.test_13vs14, alpha=alpha, pAdjustMethod="BH")

res.DESeq2_13vs14 <- res.DESeq2_13vs14[order(res.DESeq2_13vs14$padj),]

## Draw an histogram of the p-values
png("histogram_13vs14.png", width = 8, height = 4, units = 'in', res = 300)
hist(res.DESeq2_13vs14$padj, breaks=20, col="grey", main="DESeq2 p-value distribution", xlab="DESeq2 P-value", ylab="Number of genes")
dev.off()
###########Volcano plot

alpha <- 0.01 # Threshold on the adjusted p-value
cols_13vs14 <- densCols(res.DESeq2_13vs14$log2FoldChange, -log10(res.DESeq2_13vs14$pvalue))

png("volcano_13vs14.png", width = 8, height = 4, units = 'in', res = 300)
plot(res.DESeq2_13vs14$log2FoldChange, -log10(res.DESeq2_13vs14$padj), col=cols_13vs14, panel.first=grid(),
     main="Volcano plot", xlab="Effect size: log2(fold-change)", ylab="-log10(adjusted p-value)",
     pch=20, cex=0.6)
abline(v=0)
abline(v=c(-1,1), col="brown")
abline(h=-log10(alpha), col="brown")
gn.selected_13vs14 <- abs(res.DESeq2_13vs14$log2FoldChange) > 2 & res.DESeq2_13vs14$padj < alpha 
text(res.DESeq2_13vs14$log2FoldChange[gn.selected_13vs14],
     -log10(res.DESeq2_13vs14$padj)[gn.selected_13vs14],
     lab=rownames(res.DESeq2_13vs14)[gn.selected_13vs14 ], cex=0.4)
dev.off()

### Check the expression levels of the most differentially expressed gene

gn.most.sign_13vs14 <- rownames(res.DESeq2_13vs14)[1]
gn.most.diff.val_13vs14 <- counts(dds_13vs14, normalized=T)[gn.most.sign_13vs14,]

png("DFgenes_13vs14.png", width = 8, height = 4, units = 'in', res = 300)
barplot(gn.most.diff.val_13vs14, col=metadata_13vs14, main=gn.most.sign_13vs14, las=2, cex.names=0.5)
dev.off()

### MA plot

## Draw a MA plot.
## Genes with adjusted p-values below 1% are shown
png("MA_13vs14.png", width = 8, height = 8, units = 'in', res = 300)

plotMA(res.DESeq2_13vs14, colNonSig = "blue")
abline(h=c(-1:1), col="red")
dev.off()


# functional enrichment analysis using the list of induced genes. 
# This step will be performed using the gProfileR R library.


epsilon <- 1 # pseudo-count to avoid problems with log(0)

res.DESeq2.df_13vs14 <- na.omit(data.frame(res.DESeq2_13vs14))
induced.sign_13vs14 <- rownames(res.DESeq2.df_13vs14)[res.DESeq2.df_13vs14$log2FoldChange >= 2 &  res.DESeq2.df_13vs14$padj < alpha]
# head(induced.sign)
# names(term.induced)

term.induced_13vs14 <- gprofiler(query=induced.sign_13vs14, organism="scerevisiae")
term.induced_13vs14 <- term.induced_13vs14[order(term.induced_13vs14$p.value),]
# term.induced$p.value
kable(term.induced_13vs14[1:10,c("term.name",
                                 "term.size",
                                 "query.size",
                                 "overlap.size",
                                 "recall",
                                 "precision",
                                 "p.value", 
                                 "intersection")], 
      format.args=c(engeneer=TRUE, digits=3), caption="**Table: functional analysis wit gProfileR. ** ")

## now, using the list of repressed genes.

res.DESeq2.df_13vs14 <- na.omit(data.frame(res.DESeq2_13vs14))
repressed.sign_13vs14 <- rownames(res.DESeq2.df_13vs14)[res.DESeq2.df_13vs14$log2FoldChange <= -2 &  res.DESeq2.df_13vs14$padj < alpha]
head(repressed.sign_13vs14)


term.repressed_13vs14 <- gprofiler(query=repressed.sign_13vs14, organism="scerevisiae")
term.repressed_13vs14 <- term.repressed_13vs14[order(term.repressed_13vs14$p.value),]
kable(head(term.induced_13vs14[,c("p.value", "term.name","intersection")], 10))



#design(dds_13vs14) <- ~ Preservation + Site

dds_13vs14 <- DESeq(dds_13vs14)

res_13vs14 <-results(dds_13vs14)


summary(res_13vs14)


resBigFC_13vs14 <- results(dds_13vs14, lfcThreshold=1, altHypothesis="greaterAbs")
resSort_13vs14 <- res_13vs14[order(res_13vs14$pvalue),]

head(resSort_13vs14)


## visualizing results of the DF genes

# Expression Heatmap

normlzd_dds_13vs14 <- counts(dds_13vs14, normalized=T)

res_sig_13vs14 <- data.frame(normlzd_dds_13vs14[resSort_13vs14$pvalue, ])

# write.csv(res_sig_13vs14,file="res_sig_13vs14")

log_dds_13vs14<-rlog(dds_13vs14)

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
png("scatter_13vs14.png", width = 8, height = 8, units = 'in', res = 300)

pairs(log2(dat1_13vs14[,sample(ncol(dat1_13vs14), nb.pairs)] + epsilon), 
      panel=plotFun, lower.panel = NULL)
dev.off()

































# Top Differential Expressed Genes sorted by p-value

resSort_13vs14 <- res.DESeq2_13vs14[order(res.DESeq2_13vs14$pvalue),]
head(resSort_13vs14)
sum(res.DESeq2_13vs14$pvalue < 0.05, na.rm=TRUE)
sum(!is.na(res.DESeq2_13vs14$pvalue))
resSig_13vs14 <- subset(res.DESeq2_13vs14, padj < 0.5)

# No significant genes significantly different between Freeze x RNA Later

head(resSig_13vs14[ order(resSig_13vs14$log2FoldChange), ])
head(resSig_13vs14[ order(resSig_13vs14$log2FoldChange, decreasing = TRUE), ])

# visualizing results of the DF genes with plots  ------------------------------- No DF GENES THEREFORE, NO PLOTS

topGene_13vs14 <- rownames(res.DESeq2_13vs14)[which.min(res.DESeq2_13vs14$padj)]

geneCounts_13vs14 <- plotCounts(dds_13vs14, gene = topGene_13vs14, intgroup = c("Site"),
                                returnData = TRUE)
png("DF_genecounts_13vs14.png", width = 8, height = 4, units = 'in', res = 300)
ggplot(geneCounts_13vs14, aes(x = "site", y = count, color = "site")) +
  scale_y_log10() +  geom_beeswarm(cex = 3)
dev.off()




#### if the site doesnt work in the plots change for the condition



#### no results because there is no DF genes

png("DF_RLater_vs_Freeze_gene_ID.png", width = 8, height = 4, units = 'in', res = 300)
ggplot(geneCounts_13vs14, aes(x = "site", y = count, color = "site")) +
  scale_y_log10() +  geom_beeswarm(cex = 3)
dev.off()

png("DF2_RLater_vs_Freeze_gene_ID.png", width = 8, height = 4, units = 'in', res = 300)
ggplot(geneCounts_13vs14, aes(x = metadata_13vs14$Site, y = count, color = metadata_13vs14$Site, group = metadata_13vs14$Site)) +
  scale_y_log10() + geom_point(size = 3) + geom_line()
dev.off()
######### gene clustering

# 20 genes with the highest variance across samples - Not necessarily statistical significant expressed

### regularized-logarithm transformation or rlog

# regularized-logarithm transformation
rlog_13vs14 <- rlog(dds_13vs14,blind = F)
# extract the rlog matris from the object
rlog_13vs14_mat <- assay(rlog_13vs14)
# compute pairwise correlation values
rlog_13vs14_cor <- cor(rlog_13vs14_mat)
# rlog_cor


topVarGenes_13vs14 <- head(order(rowVars(assay(rlog_13vs14)), decreasing = TRUE), 20)

mat_13vs14  <- assay(rlog_13vs14)[ topVarGenes_13vs14, ]
mat_13vs14  <- mat_13vs14 - rowMeans(mat_13vs14)
anno_13vs14 <- as.data.frame(colData(rlog_13vs14)[, c("Site","Preservation")])

png("top_20_genes_13vs14.png", width = 8, height = 4, units = 'in', res = 300)
pheatmap(mat_13vs14, annotation_col = anno_13vs14)
dev.off()
