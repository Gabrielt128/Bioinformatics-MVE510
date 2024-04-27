setwd("D:/CTH/bioinformatics/exercise3/")
library(ggplot2)
# install.packages("reshape2")
library(reshape2)
library(gplots)
library(RColorBrewer)

x = read.table('counts_matrix.txt')

# number of genes, sample
num_genes = length(x[,1]) 
num_samples = length(x[1,])

metadata = read.table('metadata.txt', sep = '\t', header = TRUE)
head(metadata)
summary(metadata)

# how many male/female how many have the disease?
table(metadata$Sex)
table(metadata$diagnosis)

rownames(x)
colnames(x)
metadata$patient.id
colnames(x) == metadata$patient.id

# two specific gene
x['ENSG00000002586',]
x['ENSG00000004809',]

# filter function pt
count.nonzero = function (gene) {
  return(sum(gene > 0))
}
selected_rows = (apply(x, 1, count.nonzero) >= 20)
x.filtered = x[selected_rows,]
sum(selected_rows)

# 3.2 CPM
CPM = function (sample) {
  sum_sample = sum(sample)
  cpm = (sample + 1)/sum_sample * 10^6
  return(cpm)
}
# x.new = x.filtered + 1
x.new = apply(x.filtered, 2, CPM)
x.new.log = data.frame(log(x.new))
x.filtered.log = log(x.filtered+1)

data.new = melt(x.new.log)
data.old = melt(x.filtered.log)

# 3.3 plot
p1 = (ggplot(data.new, aes(x = variable, y = value))
      + geom_boxplot()
      + theme_bw()
      + labs(title = 'gene expression', y = 'LogCPM')
      + theme(axis.text.x = element_blank(),  
              axis.ticks.x = element_blank(),  
              axis.title.x = element_blank(),
              plot.margin = margin(1, 1, 1, 1, "cm")
      )
      + ylim(-10, 20)
)
# p1

p2 = (ggplot(data.old, aes(x = variable, y = value))
      + geom_boxplot()
      + theme_bw()
      + labs(title = 'gene expression', y = 'LogCount')
      + theme(axis.text.x = element_blank(),  
              axis.ticks.x = element_blank(),  
              axis.title.x = element_blank(),
              plot.margin = margin(1, 1, 1, 1, "cm")
      )
      + ylim(-10, 20)
)
# p2

# p2.5 combine p1 and p2 together
data.newnold = data.frame(data.new, data.old$value)
head(data.newnold)
colnames(data.newnold) = c('sample', 'logCPM', 'logCount')
data.newnold = melt(data.newnold)
# colnames(data.newnold)
p2.5 = (ggplot(data.newnold, aes(x = sample, y = value, color = variable))
        + geom_boxplot()
        + theme_bw()
        + labs(title = 'gene expression per sample', y = 'value of the expression')
        + theme(axis.text.x = element_blank(),  
                axis.ticks.x = element_blank(),  
                axis.title.x = element_blank(),
                plot.margin = margin(0.5, 1, 1, 1, "cm"),
                # legend.position = 'top'
                # legend.title = 'expression choice'
        )
        + ylim(-10, 20)
)
p2.5

p3 = (ggplot(x.new.log, aes(x = SRR1782694, y =SRR1782687))
      + geom_point()
      + theme_bw())
# p3

# 3.4 linear model for gene1
counts.gene1 = x.filtered.log[1,]
CPM.gene1 = x.new.log[1,]

# 两个dataframe都带有rownames 不为NULL，因此不能直接作为元素填充进dataframe
# function进行合并，因为dataframe function 会将列表中的元素作为列添加到数据框中
# ，如果第一个数据框没有设置行名（rownames），直接使用 data.frame(selected_row,
# selected_column) 也是可以的，因为在这种情况下没有需要对齐的行名。

# # scheme 1 delet the row names like in the slides (not working at all)
# row.names(counts.gene1) = NULL
# row.names(CPM.gene1) = NULL
# df.gene1 = data.frame(counts.gene1, CPM.gene1, metadata$diagnosis)
# df.gene1

# scheme 2 我们把dataframe元素中的行向量转换为列向量 然后再代入
# data.frame(t(counts.gene1), t(CPM.gene1))
# data.frame(metadata$Sex, metadata$diagnosis)

df.gene1 = data.frame(t(counts.gene1), t(CPM.gene1), metadata$diagnosis)
# df.gene1
dim(df.gene1)
colnames(df.gene1) = c('logcount', 'logCPM', 'diagnosis')

# this std form df is just for the ggplot
df.gene1.std = melt(df.gene1)
# colnames(df.gene1.std) = c('diagnosis', 'variable', 'value')
colnames(df.gene1.std)

p4 = (ggplot(df.gene1.std, aes( x = variable, y = value, color = diagnosis))
      + geom_boxplot()
      + theme_bw()
      + labs(x = '', y = 'value', title = 'normalized and log-transformed 
counts grouped by diagnosis (gene TSPAN6)')
      )
p4

# now the lm time
df.gene1
df.gene1.full = data.frame(df.gene1, metadata$age.at.diagnosis, metadata$Sex)
df.gene1.full$age = df.gene1.full$metadata.age.at.diagnosis
df.gene1.full$metadata.age.at.diagnosis = NULL
df.gene1.full$sex = df.gene1.full$metadata.Sex
df.gene1.full$metadata.Sex = NULL

colnames(df.gene1.full)

df.gene1.full$diagnosis = as.factor(df.gene1.full$diagnosis)
fit1 = lm(data = df.gene1.full, logCPM ~ diagnosis)
fit1.summary = summary(fit1)

df.gene1.full$sex = as.factor(df.gene1.full$sex)
fit2 = lm(data = df.gene1.full, logCPM ~ diagnosis + sex + age)
fit2.summary = summary(fit2)

# this is how to retrieve the data
# fit2.summary
# fit1.summary$coefficients[2,]
# fit2.summary$coefficients[2,4]

# 3.5 Identification
# for the key variable
ngenes = nrow(x.new.log)
pvalue.fit1 = numeric(ngenes)
pvalue.fit2 = numeric(ngenes)

coefficient1 = numeric(ngenes)
coefficient2 = numeric(ngenes)
# for less key variable
sex.pvalue.fit2 = numeric(ngenes)
sex.coefficient2 = numeric(ngenes)
age.pvalue.fit2 = numeric(ngenes)
age.coefficient2 = numeric(ngenes)
# colnames(metadata)

for ( i in 1:ngenes){
  current.gene= x.new.log[i,]
  current.df = data.frame(t(current.gene), metadata$diagnosis, metadata$Sex, metadata$age.at.diagnosis)

  colnames(current.df) = c('expression', 'diagnosis', 'sex', 'age')
  current.df$diagnosis = as.factor(current.df$diagnosis)
  current.df$sex = as.factor(current.df$sex)

  current.fit1 = lm(data = current.df, expression ~ diagnosis)
  current.fit2 = lm(data = current.df, expression ~ diagnosis + age + sex)

  coefficient1[i] = summary(current.fit1)$coefficients[2,1]
  pvalue.fit1[i] = summary(current.fit1)$coefficients[2,4]
  coefficient2[i] = summary(current.fit2)$coefficients[2,1]
  pvalue.fit2[i] = summary(current.fit2)$coefficients[2,4]
  
  # for sex and age
  sex.coefficient2[i] = summary(current.fit2)$coefficients[4,1]
  sex.pvalue.fit2[i] = summary(current.fit2)$coefficients[4,4]
  age.pvalue.fit2[i] = summary(current.fit2)$coefficients[3,4]
  age.coefficient2[i] = summary(current.fit2)$coefficients[3,1]
}
# head(pvalue.fit1)
# head(pvalue.fit2)
# head(coefficient1)
# head(coefficient2)
# fit2.summary

# ?p.adjust
df.fit1 = data.frame(coefficient1, pvalue.fit1)
rownames(df.fit1) = rownames(x.filtered)
df.fit1$padjust = p.adjust(df.fit1$pvalue.fit1, 'BH')
# df.fit1$padjust = p.adjust(df.fit1$pvalue.fit1, 'bonferroni')
head(df.fit1)

df.fit2 = data.frame(coefficient2, pvalue.fit2)
rownames(df.fit2) = rownames(x.filtered)
df.fit2$padjust = p.adjust(df.fit2$pvalue.fit2, 'BH')
# head(df.fit2)

# identification step fit1 model
p.threshold = 0.05
# head(df.fit1)
# sum(df.fit1$padjust < 0.05)
de.fit1 = df.fit1[df.fit1$padjust < p.threshold,]
de.fit2 = df.fit2[df.fit2$padjust < p.threshold,]
# sum(de.fit2$padjust >= 0.05)
# check一下meaning of the coefficient considering it's categorical data anyway
# again the baseline is CD, and the question use the Not IBD as baseline, so 
# reverse needed here
# levels(current.df$diagnosis)
# summary(current.fit1)
# fit1.summary
# p4
upde.fit1 = de.fit1[de.fit1$coefficient1 < 0,]
downde.fit1 = de.fit1[de.fit1$coefficient1 >= 0,]
# head(upde.fit1)
# head(downde.fit1)

upde.fit2 = de.fit2[de.fit2$coefficient2 < 0,]
downde.fit2 = de.fit2[de.fit2$coefficient2 >= 0,]
# head(upde.fit2)
# head(downde.fit2)

# the most significant ones
head(de.fit1)
top5.fit1 = de.fit1[order(de.fit1$padjust)[1:5],]
top5.fit2 = de.fit2[order(de.fit2$padjust)[1:5],]

# for sex pt
# sex.pvalue.fit2
# sex.coefficient2
df.sex2 = data.frame(sex.coefficient2, sex.pvalue.fit2)
dim(df.sex2)
rownames(df.sex2) = rownames(x.filtered)
head(df.sex2)
df.sex2$padjust = p.adjust(df.sex2$sex.pvalue.fit2, 'BH')
de.sex.fit2 = df.sex2[df.sex2$padjust < p.threshold,]
de.sex.fit2[order(de.sex.fit2$padjust)[1:5],]
top1gene.sex = de.sex.fit2[order(de.sex.fit2$padjust)[1],]
top1gene.sex
dim(de.sex.fit2)

# age.pvalue.fit2 = numeric(ngenes)
# age.coefficient2 = numeric(ngenes)
df.age2 = data.frame(age.coefficient2, age.pvalue.fit2)
rownames(df.age2) = rownames(x.filtered)
df.age2$padjust = p.adjust(df.age2$age.pvalue.fit2, 'BH')
de.age.fit2 = df.age2[df.age2$padjust < p.threshold,]
de.age.fit2[order(de.age.fit2$padjust)[1:5],]

# 3.6 Clustering
# head(de.)
top100gene = rownames(de.fit2[order(de.fit2$padjust)[1:100],])
xsig = x.new.log[top100gene,]
xsig = as.matrix(xsig)

mypalette = brewer.pal(11, 'RdYlBu')
morecols = colorRampPalette(mypalette)
mycols = rev(morecols(255))

column.cols = c('red', 'blue')[factor(metadata$diagnosis)]
# levels(factor(metadata$diagnosis))
column.cols
pdf('top100sigGenesHeatmap.pdf', height = 18, width = 12)
heatmap.2(xsig, trace = 'none', col = mycols, main = 'The 100 most significant genes',
          ColSideColors = column.cols)
dev.off()

# 3.7PCA
head(x.new.log)
pca = prcomp(t(x.new.log))
summary(pca)
pca$x
top3pc = pca$x[,1:3]
df.top3pc = as.data.frame(top3pc)
class(top3pc)
df.top3pc$diagnosis = metadata$diagnosis
df.top3pc$sex = metadata$Sex
head(df.top3pc)
p5 = (ggplot(df.top3pc, aes( x = PC2, y = PC3, color = diagnosis ))
      + geom_point()
      + theme_bw()
)
p5

p6 = (ggplot(df.top3pc, aes( x = PC2, y = PC3, color = sex))
      + geom_point()
      + theme_bw())
p6

# effect size and padjust for the other top 5
log2(-top5.fit1$coefficient1)
top5.fit1$pvalue.fit1

log2(-top5.fit2$coefficient2)
top5.fit2$pvalue.fit2
