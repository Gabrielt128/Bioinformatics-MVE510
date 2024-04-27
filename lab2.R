library(ggplot2)
load('genome1.rdata')
load('genome2.rdata')
load('genome3.rdata')

# 2.6 converage
cov1 = apply(genome1[,2:5], 1, sum)
cov2 = apply(genome2[,2:5], 1, sum)
cov3 = apply(genome3[,2:5], 1, sum)

mean1 = mean(cov1)
max1 = max(cov1)
max2 = max(cov2)
max3 = max(cov3)

data_range = numeric(3)
data_range[1] = max(cov1) - min(cov1)
data_range[2] = max(cov2) - min(cov2)
data_range[3] = max(cov3) - min(cov3)

sd = numeric(3)
sd[1] = sd(cov1)
sd[2] = sd(cov2)
sd[3] = sd(cov3)

mylist = list(cov1, cov2, cov3)
plotlist = list()

indexOfMax = which(cov1 == max1)

for(i in 1:3){
  subset1 = mylist[[i]][1:1000]
  subset2 = mylist[[i]][4640653:4641652]
  subset = c(subset1, subset2)
  
  dfsubset = as.data.frame(subset)
  indexOfSubset = c(rep(1,1000), rep(2,1000))
  dfsubset$index = indexOfSubset
  dfsubset$index = as.character(dfsubset$index)
  dfsubset$number = c(1:(length(subset)/2), 1:(length(subset)/2))
  dfsubset$coverage = dfsubset$subset
  
  p = (ggplot(dfsubset, aes(x = number, y = coverage, fill = index))
        + geom_bar(stat = 'identity')
        + facet_grid(.~index)
        + scale_fill_discrete('red', 'blue')
        + theme_classic()
        + labs(x = 'index'))
  
  plotlist[[i]] = p
}

# 2.7 error rate pt.
# Calculate the length of genome1
genome1.length = nrow(genome1)
genome2.length = genome1.length
genome3.length = genome1.length
# Allocate an empty vector for matches
matches = vector(length = genome1.length)
matches2 = numeric(genome2.length)
matches3 = numeric(genome3.length)

# Syntax of for-loop in R
for (pos in 1:genome1.length){
  # Get the number of reads that have the same nucleotide
  # in position ‘pos’ as the reference
  matches[pos]=genome1[pos,reference[pos]]
} 

for (pos in 1:genome2.length){
  # Get the number of reads that have the same nucleotide
  # in position ‘pos’ as the reference
  matches2[pos]=genome2[pos,reference[pos]]
}

for (pos in 1:genome1.length){
  # Get the number of reads that have the same nucleotide
  # in position ‘pos’ as the reference
  matches3[pos]=genome3[pos,reference[pos]]
}


# the proportion of reads that do not match the reference in genome1.
erate = (cov1 - matches) / cov1

# How many positions have at least one read with a mismatch?
length(matches[matches < cov1])

# Visualize the proportion of mismatching reads over a region covering 1,000 positions.
subset_erate = erate[15680:16679]
dfsubset_erate = as.data.frame(subset_erate)
dfsubset_erate$index = 1:length(subset_erate)
dfsubset_erate$erate = dfsubset_erate$subset_erate
head(dfsubset_erate)

p2 = (ggplot(dfsubset_erate, aes(x = index, y = erate))
      + geom_point(stat = 'identity')
      + labs(y = 'mismatching rate'))

p2

# -----------------
# the btest function
p0 = 0.01
mismatch1 = cov1 - matches
mismatch2 = cov2 - matches2
mismatch3 = cov3 - matches3

head(mismatch1)

btest = function(mismatch, coverage){
  if(coverage == 0){
    return(0)
  }
  
  b = binom.test(mismatch, coverage, p0, alternative = "greater")
  pvalue = b$p.value
  return(pvalue) 
}

# till the end of genome1
pvalue1 = numeric(genome1.length)
for (pos in 1:genome1.length) {
  pvalue1[pos] = btest(mismatch1[pos], cov1[pos])
}
head(pvalue1)

# till the end of genome2
pvalue2 = numeric(genome2.length)
for (pos in 1:genome2.length) {
  pvalue2[pos] = btest(mismatch2[pos], cov2[pos])
}
head(pvalue2)

# till the end of genome3
pvalue3 = numeric(genome3.length)
for (pos in 1:genome3.length) {
  pvalue3[pos] = btest(mismatch3[pos], cov3[pos])
}
head(pvalue3)

# the mutations of 1, 2, 3
plessthan5 = numeric(3)
theone = numeric(3)

plessthan5[1] = length(pvalue1[pvalue1 <= 0.05])
theone[1] = order(pvalue1)[1]

plessthan5[2] = length(pvalue2[pvalue2 <= 0.05])
theone[2] = order(pvalue2)[1]

plessthan5[3] = length(pvalue3[pvalue3 <= 0.05])
theone[3] = order(pvalue3)[1]

# double check
# pvalue1[theone[1]]
# order(pvalue1)[1] == theone[1]
# genome1[theone[1],]
# reference[theone[1]]
# binom.test(17, 17, p = 0.01, alternative = 'greater')
# binom.test(17, 17, p = 0.01, alternative = 'greater')$p.value
# pvalue1[order(pvalue1)[1]]

# 2.10
# genom2
# numbers <- seq(2339420, 2336793, 1, by = -1)`
# 
# index <- which(numbers == 2339173)
# 
# start_index <- ((length(numbers) - index) %/% 3) * 3 + 1
# 
# unit_numbers <- numbers[start_index:(start_index + 2)]
# 
# print(unit_numbers)


