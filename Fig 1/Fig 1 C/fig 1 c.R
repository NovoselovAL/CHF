
###########################################################################################################################################################
#                                                                                                                                                         #
#                                                                Script for "figure 1 c" of paper                                                         # 
#  Microbiome changes in chronic heart failure with preserved ejection fraction patients correlate with fibrosis markers: description of Russian cohort   #
# R version 3.6.1 (2019-07-05)                                                                                                                            #
###########################################################################################################################################################


install.packages("RCurl")

library( RCurl )
library(ggplot2)

setwd("C:/R project/")
getwd()


#Load the abundance table 

data1 <- read.csv("C:/path to/abundance table_1c.csv", header=TRUE, sep=";", row.name = 1,)



dim(data1)
data1[1:4,1:6]

#Normalise to library size
data <- data1/colSums(data1)[col(data1)] 
data[1:4,1:6]

#replace 0 to 1

data[data == 0] <- 1

#Data -Log2 transformation

LogOTU_norm <- log2(data)
dim(LogOTU_norm)
LogOTU_norm[1:4,1:6]


#Subtract data

CHFind <- read.csv("C:/path to/metadata CHF.csv", header=TRUE, sep=";")
Helthind <- read.csv("C:/path to/metadata Healthy.csv", header=TRUE, sep=";")
CHFmeta <- read.csv("C:/path to/metadata Groups.csv", header=TRUE, sep=";")



CHFset <- as.matrix(CHFind)
CHFpath  <- LogOTU_norm[ ,CHFset]

Helthset <- as.matrix(Helthind)
Helthpath <- LogOTU_norm[ ,Helthset]




Helthpath.mean = apply(Helthpath, 1, mean)
CHFpath.mean = apply(CHFpath, 1, mean)

fold = CHFpath.mean - Helthpath.mean

pvalue = NULL 
tstat = NULL

for(i in 1 : nrow(LogOTU_norm)) { # For each gene : 
  x = Helthpath[i,] # WT of gene number i
  y = CHFpath[i,] # KO of gene number i
  
  # Compute t-test between the two conditions
  t = t.test(x, y)
  
  # Put the current p-value in the pvalues list
  pvalue[i] = t$p.value
  # Put the current t-statistic in the tstats list
  tstat[i] = t$statistic
}

head(pvalue)

# Histogram of p-values (-log10)
hist(-log10(pvalue), col = "gray")

plot(fold, -log10(pvalue), main = "Volcano")

fold_cutoff = 1
pvalue_cutoff = 0.05
abline(v = fold_cutoff, col = "blue", lwd = 3)
abline(v = -fold_cutoff, col = "red", lwd = 3)
abline(h = -log10(pvalue_cutoff), col = "green", lwd = 3)

# Fold-change filter for "biological" significance
filter_by_fold = abs(fold) >= fold_cutoff
dim(CHFmeta[filter_by_fold, ])

# P-value filter for "statistical" significance
filter_by_pvalue = pvalue <= pvalue_cutoff
dim(CHFmeta[filter_by_pvalue, ])

# Combined filter (both biological and statistical)
filter_combined = filter_by_fold & filter_by_pvalue

filtered = CHFmeta[filter_combined,]

dim(filtered)

head(filtered)

# Let's generate the volcano plot again,
# highlighting the significantly differentially expressed genera
plot(fold, -log10(pvalue), main = "Volcano #2")
points (fold[filter_combined], -log10(pvalue[filter_combined]),
        pch = 16, col = "red")
abline(v = fold_cutoff, col = "blue", lwd = 3)
abline(v = -fold_cutoff, col = "red", lwd = 3)
abline(h = -log10(pvalue_cutoff), col = "green", lwd = 3)
# Highlighting up-regulated in blue and down-regulated in red
plot(fold, -log10(pvalue), main = "Volcano #3")
points (fold[filter_combined & fold < 0],
        -log10(pvalue[filter_combined & fold < 0]),
        pch = 16, col = "red")
points (fold[filter_combined & fold > 0],
        -log10(pvalue[filter_combined & fold > 0]),
        pch = 16, col = "blue")
abline(v = fold_cutoff, col = "blue", lwd = 3)
abline(v = -fold_cutoff, col = "red", lwd = 3)
abline(h = -log10(pvalue_cutoff), col = "green", lwd = 3)

################

#Generate output table

autput.table <- data.frame(Helthpath.mean,CHFpath.mean,fold,pvalue)
dim(autput.table)
autput.table[1:8,1:4]

#Filter Up regulated genera

UpGeneAuT <- autput.table[(autput.table$fold > fold_cutoff) & (autput.table$pvalue < pvalue_cutoff),]
dim(UpGeneAuT)
UpGeneAuT[1:15,1:4]

#Filter Down regulated genera

DownGeneAuT <- autput.table[(autput.table$fold < -fold_cutoff) & (autput.table$pvalue < pvalue_cutoff),]
dim(DownGeneAuT)
DownGeneAuT[1:15,1:4]

# Make list for labeling in volcano plot

LableSet <- rbind(DownGeneAuT,UpGeneAuT)

dim(LableSet)
LableSet[1:10,1:4]

# Highlighting up-regulated in red and down-regulated in blue with lable spesific points
plot(-log10(pvalue)~fold, xlab = 'log2(fold change)', ylab = '-log10(pvalue)', main = 'Health VS HCF Volcano #4', data = autput.table,xlim=c(-3,3))
points (fold[filter_combined & fold < 0],
                   -log10(pvalue[filter_combined & fold < 0]),
                   pch = 16, col = "blue")
points (fold[filter_combined & fold > 0],
                   -log10(pvalue[filter_combined & fold > 0]),
                   pch = 16, col = "red")
abline(v = fold_cutoff, col = "red", lwd = 3)
abline(v = -fold_cutoff, col = "red", lwd = 3)
abline(h = -log10(pvalue_cutoff), col = "green", lwd = 3)
with(LableSet, text(-log10(pvalue)~fold, labels = row.names(LableSet), pos = 2))

############
