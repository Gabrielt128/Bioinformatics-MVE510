setwd('D:/CTH/bioinformatics/exercise4')

library(ggplot2)
# install.packages("reshape2")
library(reshape2)
library(gplots)
library(RColorBrewer)
library(dplyr)

counts = read.table('16s_counts.txt', comment.char = '', quote = '', sep = '\t', header = TRUE)
# head(counts)
annotation = read.table('16s_annotation.txt', comment.char = '', quote = '', sep = '\t', header = TRUE)
# head(annotation)

# 4.1 -------------------------------------------------------------------
# sum(counts[2,])
# sum(colSums(counts))
sum(counts)
colSums(counts)
dim(annotation)
head(annotation)
# the answer is in Jan21_6
head(counts)
counts.filter = counts[rowSums(counts) >= 5,]
# dim(counts[rowSums(counts) >= 5,])
# dim(counts)

# 4,2 --------------------------------------------------------------------
# trying to conduct clustering first
dist.mat = dist(t(counts.filter))
res.cluster = hclust(dist.mat)
plot(res.cluster)

# then the PCA
res.pca = prcomp(t(counts.filter), center = FALSE)
summary(res.pca)
res.pca$x
plot(res.pca)
# we only keep the first 3 PC
top3pc = as.data.frame(res.pca$x[,1:3]) 
top3pc$group = c('HC', 'HC', 'HC', 'LC', 'LC', 'LC')
top3pc$groupname = colnames(counts)
p.pca = (
  # ggplot(top3pc, aes(x = PC2, y = PC3, color = group, label = groupname))
  ggplot(top3pc, aes(x = PC1, y = PC2, color = group, label = groupname))
         + geom_point()
         # + geom_text(data, color = 'black')
         + geom_text(aes(label = groupname), nudge_y = 200)
         + theme_bw()
         # + theme_minimal()
         )
# dim(top2pc)
p.pca
counts.filter.log = log(counts.filter+1)
# redo the cluster
dist.mat.log = dist(t(counts.filter.log))
res.cluster2 = hclust(dist.mat.log)
plot(res.cluster2)
# yes it is better now, we can calculate log(8000) to demonstrate it

# 4.3 ----------------------------------
# the function take 
# 1)a character vector of OTU name
# 2)a vector of counts for the OTUs
# 3)number indicating the resulting seq depth
# lets set the sample size 20000
sampling = function(OTUs, counts, size) {
  reads = rep(OTUs, times = counts)
  reads.sample = sample(reads, size = size, replace = FALSE)
  # kk = as.data.frame(table(readssample))
  return(as.data.frame(table(reads.sample)))
}
 # length(rep(rownames(counts.filter),counts.filter$HC1))
 # colSums(counts.filter)
test = sampling(rownames(counts.filter),counts.filter$HC1, 20000)
dim(test)
# try to use it to all six samples
counts.filter.rarefy = apply(counts.filter, 2, sampling, OTUs = rownames(counts.filter), size = 20000)
class(counts.filter.rarefy$HC1)
sapply(counts.filter.rarefy, dim)
rowCounts(counts.filter != 0)

# 4.4 richness and evenness

sapply(counts.filter.rarefy, nrow)
# head(counts.filter.rarefy$HC1)

diversity = function(count){
  rich = nrow(count)
  whole = sum(count[,2])
  # whole = colSums(count)[2]
  # we use shannon here
  H = - sum(count[,2]/whole*log(count[,2]/whole))
  return(c(rich, H))
}

sapply(counts.filter.rarefy, diversity)
colSums(counts.filter != 0)
# for double check richness
sapply(counts.filter.rarefy, nrow)
# for double check the validity of the data
colSums(counts.filter != 0)
colSums(counts.filter)
colSums(counts != 0)

# 4.5 
library(DESeq2)
design.matrix = data.frame(exposure = factor(c(1,1,1,0,0,0)))
counts.ds = DESeqDataSetFromMatrix(countData = counts.filter, design.matrix, design = ~exposure)
res.ds = DESeq(counts.ds)
result = results(res.ds, independentFiltering = FALSE, cooksCutoff = FALSE)
result.df = as.data.frame(result)

annotation.filter = annotation[rowSums(counts) >= 5,]
rownames(annotation.filter) = rownames(result.df)
# comb16s = rbind(annotation.filter, result.df)
# need to work out how to combine them together

sum(row.names(annotation.filter) == row.names(result.df))
sum(rownames(annotation.filter) == annotation.filter$OTU.ID)
head(annotation.filter)
head(result.df)

head(merge(annotation.filter, result.df, by = 0))
OTU.full = merge(annotation.filter, result.df, by = 0)
# dim(OTU.full)
# the meaning of adjustp
head(OTU.full$pvalue * 4020/rank(OTU.full$pvalue))
head(OTU.full$padj)
# how do you interpret padjust?
# not exactly BH because seems like 省略了一些数据
# because several versions of BH 
# https://www.r-bloggers.com/2023/07/the-benjamini-hochberg-procedure-fdr-and-p-value-adjusted-explained/
otu.da.sum = sum(result$padj <= 0.05)


# 4.6 top10
OTU.full[OTU.full$padj<=0.05,][1:10,]
OTU.full[OTU.full$padj<=0.2,][(order(OTU.full[OTU.full$padj<=0.2,]$padj)[1:57]),][,c('Order','Family','log2FoldChange','padj')]
df[, c("Column1", "Column2", "Column3")]
# What is log2foldchange?
# why otu4343 less than otu4325?
OTU.full[OTU.full$padj<=0.2,][(order(OTU.full[OTU.full$padj<=0.2,]$padj)[1:57]),][
  ,c('OTU.ID','Order','Family','log2FoldChange','padj')][c(1,3,19,30,44),]

# 4.7 
# shotgun dataset is a lot bigger than the amplicon amplicon只是shotgun众多基因里面
# 一个基因的各种子集 so the data volume of the two files
# do make sense
counts.shotgun = read.table('gene_counts.txt')
annotation.shotgun = read.table('gene_annotation.txt', comment.char = '', quote = '', sep = '\t', header = TRUE)
dim(counts.shotgun)
dim(annotation.shotgun)
colSums(counts.shotgun)
# the threshold for filtering is not set and 5 is too small here
# one way to conduct is the rate the ratio
sum(counts)
sum(counts.shotgun)
threshold = sum(counts.shotgun)/sum(counts) * 5
# threshold = 5
counts.shotgun.filter = counts.shotgun[rowSums(counts.shotgun) >= threshold,]

# counts.filter = counts[rowSums(counts) >= 5,]

# lets try cluster
# dist.mat = dist(t(counts.filter))
# res.cluster = hclust(dist.mat)
# plot(res.cluster)
dist.mat.shotgun = dist(t(counts.shotgun.filter))
res.cluster.shotgun = hclust(dist.mat.shotgun)
plot(res.cluster.shotgun)
# try log transfer and see if its getting better
counts.shotgun.filter.log = log(counts.shotgun.filter + 1)
dist.mat.shotgun.log = dist(t(counts.shotgun.filter.log))
res.cluster.shotgun.log = hclust(dist.mat.shotgun.log)
plot(res.cluster.shotgun.log)


# 4.8
colSums(counts.shotgun.filter)
counts.shotgun.filter.rarefy = apply(
  counts.shotgun.filter, 2, sampling, OTUs = rownames(counts.shotgun.filter), size = 530000)
sapply(counts.shotgun.filter.rarefy, diversity)
dim(counts.shotgun.filter.rarefy$hc3)

# 4.9
design.matrix.shotgun = data.frame(exposure = factor(c(1,1,1,0,0,0)))
count.ds.shotgun = DESeqDataSetFromMatrix(
  countData = counts.shotgun.filter, design.matrix.shotgun, design = ~exposure)
res.ds.shotgun = DESeq(count.ds.shotgun)
result.shotgun = results(res.ds.shotgun, independentFiltering = FALSE, cooksCutoff = FALSE)
result.df.shotgun = as.data.frame(result.shotgun)
annotation.shotgun.filter = annotation.shotgun[rowSums(counts.shotgun) >= threshold,]
rownames(annotation.shotgun.filter) = rownames(result.df.shotgun)
gene.full = merge(annotation.shotgun.filter, result.df.shotgun, by = 0)

sum(gene.full$padj <= 0.05)
gene.full[order(gene.full$padj),][1,]
gene.full[grep("pdxA", gene.full$Description), ]
annotation.shotgun[grep("pdxA", annotation.shotgun$Description),]
