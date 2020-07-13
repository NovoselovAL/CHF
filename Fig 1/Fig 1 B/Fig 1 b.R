

###########################################################################################################################################################
#                                                                                                                                                         #
#                                                                Script for "figure 1 b" of paper                                                         # 
#  Microbiome changes in chronic heart failure with preserved ejection fraction patients correlate with fibrosis markers: description of Russian cohort   #
# R version 3.6.1 (2019-07-05)                                                                                                                            #
# Script for plotting Heatmap with OTU filtering using source script                                                                                      #
# https://www.molecularecologist.com/2013/08/making-heatmaps-with-r-for-microbiome-analysis/                                                              #
# cut of 0.1 (less than 10% cutted)                                                                                                                       #
###########################################################################################################################################################

# Instal librarys
install.packages("RCurl")
install.packages("dplyr")
install.packages("data.table")
install.packages("vegan")
install.packages("RColorBrewer")
install.packages("gplots")
# to install packages from Bioconductor:
source("http://bioconductor.org/biocLite.R")
biocLite("Heatplus")  # annHeatmap or annHeatmap2


# Load librarys
library( RCurl )
library(ggplot2)
library(dplyr)
library(data.table)
library(vegan)
library(gplots)  # for heatmap.2
library(Heatplus)
library(RColorBrewer)# load the RColorBrewer package for better colour options


# Set working directory
setwd("C:/R project/")
getwd()

all.data <- read.csv("C:/path to/abundance table_Fig_1b.csv", header=TRUE, sep=";", row.names=1)  # load the data 
dim(all.data)

all.data[1:3, 1:4]
row.names(all.data) <- all.data$sample # Probably wy do not need this, depend on input data

all.data <- all.data[, -1]

data.prop <- all.data/rowSums(all.data)

data.prop[1:3, 1:3]
# colorRampPalette is in the RColorBrewer package.  This creates a colour palette that shades from light yellow to red in RGB space with 100 unique colours
scaleyellowred <- colorRampPalette(c("lightyellow", "red"), space = "rgb")(100)

# The most basic of heatmaps:
  heatmap(as.matrix(data.prop), Rowv = NA, Colv = NA, col = scaleyellowred)
# determine the maximum relative abundance for each column
maxab <- apply(data.prop, 2, max)
head(maxab)

# remove the genera with less than 1% as their maximum relative abundance
n1 <- names(which(maxab < 0.1))
data.prop.1 <- data.prop[, -which(names(data.prop) %in% n1)]
# the margins command sets the width of the white space around the plot. The first element is the bottom margin and the second is the right margin
heatmap(as.matrix(data.prop.1), Rowv = NA, Colv = NA, col = scaleyellowred, margins = c(10, 2))
# calculate the Bray-Curtis dissimilarity matrix on the full dataset:
data.dist <- vegdist(data.prop, method = "bray")

# Do average linkage hierarchical clustering. Other options are 'complete' or 'single'. You'll need to choose the one that best fits the needs of your situation and your data.

row.clus <- hclust(data.dist, "aver")

# make the heatmap with Rowv = as.dendrogram(row.clus)

heatmap(as.matrix(data.prop.1), Rowv = as.dendrogram(row.clus), Colv = NA, col = scaleyellowred, margins = c(10, 3))

# you have to transpose the dataset to get the genera as rows
data.dist.g <- vegdist(t(data.prop.1), method = "bray")
col.clus <- hclust(data.dist.g, "aver")
# make the heatmap with Rowv = as.dendrogram(row.clus)
heatmap(as.matrix(data.prop.1), Rowv = as.dendrogram(row.clus), Colv = as.dendrogram(col.clus), col = scaleyellowred, margins = c(10, 3))

# Label CHF and Control with color 
SampleGroup <- read.csv("C:/path to/metadata Groups_Fig_1b.csv", header=TRUE, sep=";" )  # load the data (table from Maria)
dim(SampleGroup)

# just to remove weird column
SampleGroup <- SampleGroup[,-1]
SampleGroup <- SampleGroup[,-2]
# just to remove weird column

str(SampleGroup)
# Create vector from input data
var1 <- as.vector(SampleGroup)


heatmap.2(as.matrix(data.prop.1),Rowv = as.dendrogram(row.clus), Colv = as.dendrogram(col.clus), col = scaleyellowred, RowSideColors = var1,) # this puts in the annotation for the samples margins = c(10, 3))
heatmap.2(as.matrix(data.prop.1), Rowv = as.dendrogram(row.clus), Colv = as.dendrogram(col.clus), col = scaleyellowred, RowSideColors = var1, margins = c(11, 5), trace = "none", density.info = "none", xlab = "genera", ylab = "Samples", main = "Heatmap", lhei = c(2, 8)) # this makes the colour-key legend a little thinner

# the annHeatmap2 function needs to be wrapped in the plot function in order to display the results
# the colours have to be fiddled with a little more than with the other functions. With 50 breaks in the heatmap densities, there needs to be 51 colours in the scale palette. Error messages from the plotting should help to determine how many colours are needed for a given number of breaks
# dendrogram controls how dendrograms are made and can be calculatated witin the function. The internal calculation can be overridden and externally calculated dendrograms can be displayed.
# this puts the colour-scale legend on the plot. The number indicates the side on which to plot it (1 = bottom, 2 = left, 3 = top, 4 = right)
# gives more space for the Genus names

plot(annHeatmap2(as.matrix(data.prop.1), col = colorRampPalette(c("darkred", "darkgreen"), space = "rgb")(51), breaks = 50, dendrogram = list(Row = list(dendro = as.dendrogram(row.clus)), Col = list(dendro = as.dendrogram(col.clus))), legend = 3, labels = list(Col = list(nrow = 12))),xaxis.rot = 45) 
                 
                 
heatmap.2(as.matrix(data.prop.1), Rowv = as.dendrogram(row.clus), Colv = as.dendrogram(col.clus), col = scaleyellowred, srtCol=45, RowSideColors = var1,margins = c(11, 5), trace = "none", density.info = "none", xlab = "genera", ylab = "Samples", main = "Heatmap", lhei = c(2, 8))  # this makes the colour-key legend a little thinner
  
heatmap.2(as.matrix(t(data.prop.1)), Rowv = as.dendrogram(col.clus), Colv = as.dendrogram(row.clus), col = scaleyellowred, srtCol=45, ColSideColors = var1, margins = c(11, 5), trace = "none", density.info = "none", xlab = "Samples", ylab = "genera", main = "Heatmap", lhei = c(2, 8) )        


