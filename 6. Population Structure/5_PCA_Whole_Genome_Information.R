#################################################################
##                          PCA PLOTS                          ##
#################################################################

## Carla Coll Costa, C3. February 2021.
## This script plots PCA and saves the plots automatically.
## The input files are .eigenvec and .eigenval files calculated with plink and
## a file with the populations and individuals that are being studied.
## The script plots the PCA labelling and colouring the groups by different
## classifications specified in columns in the file with the populations studied.

# Upload the needed libraries
library("vcfR")
library("OutFLANK")
library("scales")
library("adegenet")
library("qqman")
library("PopGenome")
library("MASS")
library("ggplot2")
library("ggfortify")
library("tidyverse")
library("easyGgplot2")

# Set the working directory
setwd("<Working_Directory>")

# Get population classification
pops_all = read.table("<Information_All_Populations>", header=T, strip.white = T)
pops_notpacific = read.table("<Information_European_Populations>", heade=T, strip.white = T)

# Load the eigenvalues and eigenvectors
eigenvec_all = read_table2("<Eigenvectors_All_Populations.eigenvec>", col_names = FALSE)
eigenval_all = scan("<Eigenvalues_All_Populations.eigenval>")
eigenvec_notpacific = read_table2("<Eigenvectors_European_Populations.eigenvec>", col_names = FALSE)
eigenval_notpacific = scan("<Eigenvalues_European_Populations.eigenval>")

# Remove the repeated column with the names
eigenvec_all = eigenvec_all[,-1]
eigenvec_notpacific = eigenvec_notpacific[,-1]

# Set proper column names according to the individual
names(eigenvec_all)[1] = "Individual"
names(eigenvec_all)[2:ncol(eigenvec_all)] = paste0("PC", 1:(ncol(eigenvec_all)-1))
names(eigenvec_notpacific)[1] = "Individual"
names(eigenvec_notpacific)[2:ncol(eigenvec_notpacific)] = paste0("PC", 1:(ncol(eigenvec_notpacific)-1))

# Remake the dataframe adding names and type of populations. For this we sort the two dataframes
# by the individual names and then remove the column "individual" from one of the dataframes so
# it is not repeated.
pops_all = pops_all[order(pops_all$Individual),]
eigenvec_all = eigenvec_all[order(eigenvec_all$Individual),]
eigenvec_all = subset(eigenvec_all, select = -c(Individual))
pca_all = data.frame(pops_all, eigenvec_all)
pops_notpacific = pops_notpacific[order(pops_notpacific$Individual),]
eigenvec_notpacific = eigenvec_notpacific[order(eigenvec_notpacific$Individual),]
eigenvec_notpacific = subset(eigenvec_notpacific, select = -c(Individual))
pca_notpacific = data.frame(pops_notpacific, eigenvec_notpacific)

# Translate eigenvalues to percentage of variance explained
pvariance_all = data.frame(PC=1:20, pvariance=eigenval_all/sum(eigenval_all)*100)
pvariance_notpacific = data.frame(PC=1:20, pvariance=eigenval_notpacific/sum(eigenval_notpacific)*100)

# Percentatge of variance explained by each PC 
ggplot(pvariance_all, aes(PC, pvariance)) + 
  geom_bar(stat="identity") + 
  ylab("Percentatge of Variance Explained") + 
  ggtitle("PCA All Samples") +
  theme(plot.title = element_text(hjust=0.5))

ggplot(pvariance_notpacific, aes(PC, pvariance)) + 
  geom_bar(stat="identity") + 
  ylab("Percentatge of Variance Explained") + 
  ggtitle("PCA European Samples") +
  theme(plot.title = element_text(hjust=0.5))

# Cumulative sum of variance explained by PC
cumsum(pvariance_all$pvariance)
cumsum(pvariance_notpacific$pvariance)

# Vector of colours to use for plots
colour_vector = c("#7FC97F", "#BEAED4", "#FDC086", "#FFFF99", "#386CB0", "#F0027F", 
                  "#BF5B17", "#666666", "#1B9E77", "#D95F02", "#7570B3", "#E7298A", 
                  "#66A61E", "#E6AB02", "#A6761D", "#666666", "#A6CEE3", "#1F78B4", 
                  "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", 
                  "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928", "#FBB4AE", "#B3CDE3", 
                  "#CCEBC5", "#DECBE4", "#FED9A6", "#FFFFCC", "#E5D8BD", "#FDDAEC", 
                  "#F2F2F2", "#B3E2CD", "#FDCDAC", "#CBD5E8", "#F4CAE4", "#E6F5C9", 
                  "#FFF2AE", "#F1E2CC", "#CCCCCC", "#E41A1C", "#377EB8", "#4DAF4A", 
                  "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999", 
                  "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", 
                  "#E5C494", "#B3B3B3", "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", 
                  "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", "#D9D9D9", "#BC80BD", 
                  "#CCEBC5", "#FFED6F")

# Final plots
png("PCA All Final.png", units="in", width=8, height=6, res=900)
ggplot(pca_all, aes(PC1,PC2, col=pca_all$Concrete_Location)) + 
  geom_point(aes(shape=pca_all$Type), size=4, stroke=1, alpha=0.8) +
  xlab(paste0("PC1 (",formatC(pvariance_all$pvariance[1], digits=4),"%)")) + 
  ylab(paste0("PC2 (",formatC(pvariance_all$pvariance[2], digits=4),"%)")) +
  #geom_text(aes(label=pca_all$Individual), hjust=0, vjust=0, size=2) +
  scale_shape_manual(values=c("Freshwater"=16, "Marine"=1)) +
  scale_color_manual(name="Location", values=c("Adriatic"="orange", "Alaska"="darkolivegreen4", "North_Scandinavia"="royalblue3", "East_Rusia"="tomato4", "English_Channel"="chartreuse3", "Finnish_Golf"="purple3", "Iberian_Peninsula"="red", "North_Sea"="turquoise4"),
                     labels=c("Adriatic"="Adriatic Sea", "Alaska"="Alaska", "North_Scandinavia"="North Scandinavia", "East_Rusia"="East Russia", "English_Channel"="English Channel", "Finnish_Golf"="Gulf of Finland", "Iberian_Peninsula"="Iberian Peninsula", "North_Sea"="North Sea")) +
  guides(colour=guide_legend(order=1, title="Geographic Region", override.aes=list(size=3, alpha=1)), shape=guide_legend(order=2, title="Habitat Type",override.aes=list(size=3, alpha=1))) +
  theme(plot.title = element_text(hjust=0.5)) +
  theme_bw()
dev.off()

png("PCA Europe Final.png", units="in", width=8, height=6, res=900)
ggplot(pca_notpacific, aes(PC1,PC2, col=pca_notpacific$Concrete_Location)) + 
  geom_point(aes(shape=pca_notpacific$Type), size=4, stroke=1, alpha=0.8) +
  xlab(paste0("PC1 (",formatC(pvariance_notpacific$pvariance[1], digits=4),"%)")) + 
  ylab(paste0("PC2 (",formatC(pvariance_notpacific$pvariance[2], digits=4),"%)")) +
  #geom_text(aes(label=pca_notpacific$Individual), hjust=0, vjust=0, size=2) +
  scale_shape_manual(values=c("Freshwater"=16, "Marine"=1)) +
  scale_color_manual(name="Location", values=c("Adriatic"="orange", "North_Scandinavia"="royalblue3", "English_Channel"="chartreuse3", "Finnish_Golf"="purple3", "Iberian_Peninsula"="red", "North_Sea"="turquoise4"),
                     labels=c("Adriatic"="Adriatic Sea", "North_Scandinavia"="North Scandinavia", "English_Channel"="English Channel", "Finnish_Golf"="Gulf of Finland", "Iberian_Peninsula"="Iberian Peninsula", "North_Sea"="North Sea")) +
  guides(colour=guide_legend(order=1, title="Geographic Region", override.aes=list(size=3, alpha=1)), shape=guide_legend(order=2, title="Habitat Type",override.aes=list(size=3, alpha=1))) +
  theme(plot.title = element_text(hjust=0.5)) +
  theme_bw()
dev.off()
