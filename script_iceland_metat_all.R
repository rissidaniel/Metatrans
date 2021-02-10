library(DESeq2)
library(RColorBrewer)
library(pheatmap)
library(tidyverse)
library(gplots)
library(plyr)
library(affy)
library(glmpca)
library(knitr)

# setwd("../Desktop/iceland_metat/")
data <- read.csv("iceland_genes_new.tsv", sep = "\t", header = T)
metadata <-  read.csv("metadata2.txt", sep = "\t", row.names = 1)






dat1 <- ddply(data,"Gene.name",numcolwise(sum))


rownames(dat1) <- dat1[,1]
dat1 <- dat1[,-1]
# "10", "10", "10", "11", "11", "11", "12", "12", "12", "13", "13", "13", "14", "14", "14", "Neg", "Neg", "Neg

# condition =  factor(c("Cultive", "10", "10", "10", "11", "11", "11", "12", "12", "12", "13", "13", "13", "14", "14", "14", "Neg", "Neg", "Neg"))
# coldata <- data.frame(row.names=colnames(dat1), condition)


rownames(metadata)
colnames(dat1)

all(rownames(metadata) == colnames(dat1))


dds <- DESeqDataSetFromMatrix(countData= dat1,
                              colData = metadata,
                              design = ~ 1)

dds <- estimateSizeFactors(dds)


sizeFactors(dds)


plot(sizeFactors(dds), colSums(counts(dds)))
abline(lm(colSums(counts(dds)) ~ sizeFactors(dds) + 0))


normlzd_dds <- counts(dds, normalized=T)

head(normlzd_dds)


plot(hclust(dist(t(normlzd_dds))), labels=colData(dds)$site)

plot(hclust(dist(t(normlzd_dds))), labels=colData(dds)$Preservation)


plot(log(normlzd_dds[,1])+1, log(normlzd_dds[,2])+1, cex =.1)


# Varaiance Stabilizing transformation
vsd <- rlog(dds, blind = T)

# extract the vst matris from the object
vsd_mat <- assay(vsd)

# compute pairwise correlation values
vsd_cor <- cor(vsd_mat)

vsd_cor


pheatmap(vsd_cor)


dds <- DESeqDataSetFromMatrix(countData= dat1,
                                     colData = metadata,
                                     design = ~ 1)
keep <- rowSums(counts(dds)) > 1
dds <- dds[keep,]
nrow(dds)

### adding color to sites

col.strain <- c("10"="green","11"="orange", "12"="red", "13"="blue", "14"="yellow", "Neg"="white", "Cultive"="black" ) # Choose one color per strain
metadata$color <- col.strain[as.vector(metadata$Site)]

## Boxplots of gene count distributions per sample

boxplot(log2(dat1 + epsilon), col=metadata$color, pch=".", 
        horizontal=TRUE, cex.axis=0.5,
        las=1, ylab="Samples", xlab="Box plots of non-normalized log2(counts) per sample")

# Density plots

# Densities of log2(counts). Each curve corresponds to one sample.

epsilon <- 1

plotDensity(log2(dat1 + epsilon), lty=1, col=metadata$color, lwd=2)
grid()
legend("topright", legend=names(col.strain), col=col.strain, lwd=2)

# PCA plot


gpca <- glmpca(counts(dds), L=2)
gpca.dat <- gpca$factors
gpca.dat$Preservation <- dds$Preservation
gpca.dat$Site <- dds$Site
ggplot(gpca.dat, aes(x = dim1, y = dim2, color = Site, shape = Preservation)) +
  geom_point(size =3) + coord_fixed() + ggtitle("glmpca - Generalized PCA")

## histogram 
# http://pedagogix-tagc.univ-mrs.fr/courses/ASG1/practicals/rnaseq_diff_Snf2/rnaseq_diff_Snf2.html
hist(as.matrix(dat1), col="blue", border="white", breaks=100)

hist(as.matrix(dat1), col="blue", border="white",
     breaks=20000, xlim=c(0,2000), main="Counts per gene",
     xlab="Counts (truncated axis)", ylab="Number of genes", 
     las=1, cex.axis=0.7)

epsilon <- 1 # pseudo-count to avoid problems with log(0)
hist(as.matrix(log2(dat1 + epsilon)), breaks=100, col="blue", border="white",
     main="Log2-transformed counts per gene", xlab="log2(counts+1)", ylab="Number of genes", 
     las=1, cex.axis=0.7)

par(mfrow=c(1,1))





## normalization

################################################ FUNCTION
### Let's implement such a function
### cds is a countDataset
estimSf <- function (cds){
  # Get the count matrix
  cts <- counts(cds)
  
  # Compute the geometric mean
  geomMean <- function(x) prod(x)^(1/length(x))
  
  # Compute the geometric mean over the line
  gm.mean  <-  apply(cts, 1, geomMean)
  
  # Zero values are set to NA (avoid subsequentcdsdivision by 0)
  gm.mean[gm.mean == 0] <- NA
  
  # Divide each line by its corresponding geometric mean
  # sweep(x, MARGIN, STATS, FUN = "-", check.margin = TRUE, ...)
  # MARGIN: 1 or 2 (line or columns)
  # STATS: a vector of length nrow(x) or ncol(x), depending on MARGIN
  # FUN: the function to be applied
  cts <- sweep(cts, 1, gm.mean, FUN="/")
  
  # Compute the median over the columns
  med <- apply(cts, 2, median, na.rm=TRUE)
  
  # Return the scaling factor
  return(med)
}


## Normalizing using the method for an object of class"CountDataSet" 
dds.norm <-  estimateSizeFactors(dds)
sizeFactors(dds.norm)

## Now get the scaling factor with our homemade function.cds.norm
head(estimSf(dds)) 

all(round(estimSf(dds),6) == round(sizeFactors(dds.norm), 6))

## Checking the normalization
par(mfrow=c(2,2),cex.lab=0.7)
boxplot(log2(counts(dds.norm)+epsilon),  col=metadata$color, cex.axis=0.7, 
        las=1, xlab="log2(counts)", horizontal=TRUE, main="Raw counts")
boxplot(log2(counts(dds.norm, normalized=TRUE)+epsilon),  col=metadata$color, cex.axis=0.7, 
        las=1, xlab="log2(normalized counts)", horizontal=TRUE, main="Normalized counts") 
plotDensity(log2(counts(dds.norm)+epsilon),  col=metadata$color, 
            xlab="log2(counts)", cex.lab=0.7, panel.first=grid()) 

plotDensity(log2(counts(dds.norm, normalized=TRUE)+epsilon), col=metadata$color, 
            xlab="log2(normalized counts)", cex.lab=0.7, panel.first=grid()) 


####### Count variance is related to mean

## Computing mean and variance
norm.counts <- counts(dds.norm, normalized=TRUE)
mean.counts <- rowMeans(norm.counts)
variance.counts <- apply(norm.counts, 1, var)

## sum(mean.counts==0) # Number of completely undetected genes

norm.counts.stats <- data.frame(
  min=apply(norm.counts, 2, min),
  mean=apply(norm.counts, 2, mean),
  median=apply(norm.counts, 2, median),
  max=apply(norm.counts, 2, max),
  zeros=apply(norm.counts==0, 2, sum),
  percent.zeros=100*apply(norm.counts==0, 2, sum)/nrow(norm.counts),
  perc05=apply(norm.counts, 2, quantile, 0.05),
  perc10=apply(norm.counts, 2, quantile, 0.10),
  perc90=apply(norm.counts, 2, quantile, 0.90),
  perc95=apply(norm.counts, 2, quantile, 0.95)
)

kable(norm.counts.stats)


## Mean and variance relationship
mean.var.col <- densCols(x=log2(mean.counts), y=log2(variance.counts))
plot(x=log2(mean.counts), y=log2(variance.counts), pch=16, cex=0.5, 
     col=mean.var.col, main="Mean-variance relationship",
     xlab="Mean log2(normalized counts) per gene",
     ylab="Variance of log2(normalized counts)",
     panel.first = grid())
abline(a=0, b=1, col="brown")

#### Figure: variance/mean plot. 
#The brown line highlights x=y, which corresponds to the expected relationship between 
#mean and variance for a Poisson distribution.


###### ##################################### sample comparison_10vs11

head(dat1)
# "Cultive", "10", "10", "10", "11", "11", "11", "12", "12", "12", "13", "13", "13", "14", "14", "14", "Neg", "Neg", "Neg
dat1_10vs11 <- dat1 %>% select(2:7)
metadata_10vs11 <- slice(metadata, (2:7))


dds_10vs11 <- DESeqDataSetFromMatrix(countData= dat1_10vs11,
                              colData = metadata_10vs11,
                              design = ~ 1)

keep <- rowSums(counts(dds_10vs11)) > 1
dds_10vs11 <- dds_10vs11[keep,]
nrow(dds_10vs11)

dds_10vs11 <- estimateSizeFactors(dds_10vs11)

#design(dds_10vs11) <- ~ Preservation + Site

dds_10vs11 <- DESeq(dds_10vs11)

res_10vs11 <-results(dds_10vs11)


summary(res_10vs11)


plotMA(res_10vs11, ylim=c(-5,5) )


resBigFC_10vs11 <- results(dds_10vs11, lfcThreshold=1, altHypothesis="greaterAbs")
plotMA(resBigFC_10vs11, ylim=c(-5,5))
abline(h=c(-1,1),lwd=5)


plotDispEsts(dds_10vs11)



resSort_10vs11 <- res_10vs11[order(res_10vs11$pvalue),]

head(resSort_10vs11)


## visualizing results of the DF genes

# Expression Heatmap

normlzd_dds_10vs11 <- counts(dds_10vs11, normalized=T)

res_sig_10vs11 <- data.frame(normlzd_dds_10vs11[resSort_10vs11$pvalue, ])

# write.csv(res_sig_10vs11,file="res_sig_10vs11")


plotDispEsts(dds)


log_dds_10vs11<-rlog(dds_10vs11)
# plotPCAWithSampleNames(log_dds_10vs11, intgroup="Site")



# Densities of log2(counts). Each curve corresponds to one sample.

plotDensity(log2(dat1_10vs11 + epsilon), lty=1, col=metadata_10vs11$color, lwd=2)
grid()
legend("topright", legend=names(col.strain), col=col.strain, lwd=2)


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
pairs(log2(dat1_10vs11[,sample(ncol(dat1_10vs11), nb.pairs)] + epsilon), 
      panel=plotFun, lower.panel = NULL)





