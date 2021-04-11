#################################################################
##                        LDna and PCAs                        ##
#################################################################

## Carla Coll Costa, C3. January 2021.
## This script performs LDna analysis and plots and saves automatically the clusters
## of interest of such analysis. It also plots and saves automatically the heterozygosity
## against PC1 of those clusters.
## The expected input for LDna analysis are 012 files generated with vcftools and with
## the -1 values changed by NA values. To plot the heterozygosity against PC1 it is necessary
## to provide .het of the region of the cluster of interest.

#### Prepare the data ####

# Upload the required libraries
library(gdsfmt)
library(LDna)
library(SNPRelate)
library(data.table)
library(ggplot2)
library(tm)
library(reshape2)
library(vcfR)
library(qqman)
library(easyGgplot2)
library(dplyr)
library(viridis)
library(RColorBrewer)
library(gridExtra)
library(ggpubr)
options(expressions = 500000)

# Set the working directory
setwd("~/Desktop/Master's Thesis/Results/14. LDna_Region/Clusters")

# Read the file with the genotypes (row: SNPs, column: individual)
genotypes1 = t(read.table("~/Desktop/Master's Thesis/Results/14. LDna_Region/0_data/All/Genotypes_Filtered_2.3_Chromosome1_Region.012",header=FALSE))
genotypes9 = t(read.table("~/Desktop/Master's Thesis/Results/14. LDna_Region/0_data/All/Genotypes_Filtered_2.3_Chromosome9_Region.012",header=FALSE))
genotypes11 = t(read.table("~/Desktop/Master's Thesis/Results/14. LDna_Region/0_data/All/Genotypes_Filtered_2.3_Chromosome11_Region.012",header=FALSE))
genotypes21 = t(read.table("~/Desktop/Master's Thesis/Results/14. LDna_Region/0_data/All/Genotypes_Filtered_2.3_Chromosome21_Region.012",header=FALSE))

# Delete first row
genotypes1 = genotypes1[-1,]
genotypes9 = genotypes9[-1,]
genotypes11 = genotypes11[-1,]
genotypes21 = genotypes21[-1,]

# Save the genotypes as a matrix
genotypes1 = as.matrix(genotypes1)
genotypes9 = as.matrix(genotypes9)
genotypes11 = as.matrix(genotypes11)
genotypes21 = as.matrix(genotypes21)

# Read the file with the genotype positions
positions1 = read.table("~/Desktop/Master's Thesis/Results/14. LDna_Region/0_data/All/Genotypes_Filtered_2.3_Chromosome1_Region.012.pos", header = F)
positions9 = read.table("~/Desktop/Master's Thesis/Results/14. LDna_Region/0_data/All/Genotypes_Filtered_2.3_Chromosome9_Region.012.pos", header = F)
positions11 = read.table("~/Desktop/Master's Thesis/Results/14. LDna_Region/0_data/All/Genotypes_Filtered_2.3_Chromosome11_Region.012.pos", header = F)
positions21 = read.table("~/Desktop/Master's Thesis/Results/14. LDna_Region/0_data/All/Genotypes_Filtered_2.3_Chromosome21_Region.012.pos", header = F)

# Create a new column with an ID
positions1$LocID = paste0(positions1$V1,"_",positions1$V2)
positions9$LocID = paste0(positions9$V1,"_",positions9$V2)
positions11$LocID = paste0(positions11$V1,"_",positions11$V2)
positions21$LocID = paste0(positions21$V1,"_",positions21$V2)

# Create a vector as large as the number of genotypes
snps1 = 1:nrow(genotypes1)
snps9 = 1:nrow(genotypes9)
snps11 = 1:nrow(genotypes11)
snps21 = 1:nrow(genotypes21)

# Create a new column as large as the snps
positions1$snps = positions1$LocID
positions9$snps = positions9$LocID
positions11$snps = positions11$LocID
positions21$snps = positions21$LocID

# Read the file with the individuals
samples1 = read.table("~/Desktop/Master's Thesis/Results/14. LDna_Region/0_data/All/Genotypes_Filtered_2.3_Chromosome1_Region.012.indv", header = F, stringsAsFactors = F)
samples9 = read.table("~/Desktop/Master's Thesis/Results/14. LDna_Region/0_data/All/Genotypes_Filtered_2.3_Chromosome9_Region.012.indv", header = F, stringsAsFactors = F)
samples11 = read.table("~/Desktop/Master's Thesis/Results/14. LDna_Region/0_data/All/Genotypes_Filtered_2.3_Chromosome11_Region.012.indv", header = F, stringsAsFactors = F)
samples21 = read.table("~/Desktop/Master's Thesis/Results/14. LDna_Region/0_data/All/Genotypes_Filtered_2.3_Chromosome21_Region.012.indv", header = F, stringsAsFactors = F)

# Create a new column with the populations
samples1$population = gsub("[0-9]+","", samples1$V1)
samples1$population = gsub("\\.","", samples1$population)
samples1$population = gsub(".{1}$", "", samples1$population)
samples9$population = gsub("[0-9]+","", samples9$V1)
samples9$population = gsub("\\.","", samples9$population)
samples9$population = gsub(".{1}$", "", samples9$population)
samples11$population = gsub("[0-9]+","", samples11$V1)
samples11$population = gsub("\\.","", samples11$population)
samples11$population = gsub(".{1}$", "", samples11$population)
samples21$population = gsub("[0-9]+","", samples21$V1)
samples21$population = gsub("\\.","", samples21$population)
samples21$population = gsub(".{1}$", "", samples21$population)

# Take the file with the populations
pop_labs = read.table("~/Desktop/Master's Thesis/Results/14. LDna_Region/0_data/All/Populations_230.txt", header = T)

# Order the files by the individuals
pop_labs = pop_labs[with(pop_labs, order(pop_labs$Individual)),]

# Create .gds file
snpgdsCreateGeno("Chromosome 1.gds", genmat = genotypes1, sample.id = pop_labs$Individual, snp.id = positions1$LocID, snpfirstdim=TRUE) ### sample ID should be 1: how many individuals you have
snpgdsCreateGeno("Chromosome 9.gds", genmat = genotypes9, sample.id = pop_labs$Individual, snp.id = positions9$LocID, snpfirstdim=TRUE) ### sample ID should be 1: how many individuals you have
snpgdsCreateGeno("Chromosome 11.gds", genmat = genotypes11, sample.id = pop_labs$Individual, snp.id = positions11$LocID, snpfirstdim=TRUE) ### sample ID should be 1: how many individuals you have
snpgdsCreateGeno("Chromosome 21.gds", genmat = genotypes21, sample.id = pop_labs$Individual, snp.id = positions21$LocID, snpfirstdim=TRUE) ### sample ID should be 1: how many individuals you have

# Save gds file
file1 = snpgdsOpen("Chromosome 1.gds")
file9 = snpgdsOpen("Chromosome 9.gds")
file11 = snpgdsOpen("Chromosome 11.gds")
file21 = snpgdsOpen("Chromosome 21.gds")


#### Calculate LD ####

# Linkage Disequilibrium 
matrixx1 = snpgdsLDMat(file1, method="r",slide=0, verbose=TRUE)
LD_matrix1 = (matrixx1$LD)^2
LD_matrix1[upper.tri(LD_matrix1, diag = TRUE)] <- NA
matrixx9 = snpgdsLDMat(file9, method="r",slide=0, verbose=TRUE)
LD_matrix9 = (matrixx9$LD)^2
LD_matrix9[upper.tri(LD_matrix9, diag = TRUE)] <- NA
matrixx11 = snpgdsLDMat(file11, method="r",slide=0, verbose=TRUE)
LD_matrix11 = (matrixx11$LD)^2
LD_matrix11[upper.tri(LD_matrix11, diag = TRUE)] <- NA
matrixx21 = snpgdsLDMat(file21, method="r",slide=0, verbose=TRUE)
LD_matrix21 = (matrixx21$LD)^2
LD_matrix21[upper.tri(LD_matrix21, diag = TRUE)] <- NA

# Change row names and column names
rownames(LD_matrix1) = positions1$LocID
colnames(LD_matrix1) = positions1$LocID
rownames(LD_matrix9) = positions9$LocID
colnames(LD_matrix9) = positions9$LocID
rownames(LD_matrix11) = positions11$LocID
colnames(LD_matrix11) = positions11$LocID
rownames(LD_matrix21) = positions21$LocID
colnames(LD_matrix21) = positions21$LocID

# Look how many pairwise comparisons we have
length(LD_matrix1)
length(LD_matrix9)
length(LD_matrix11)
length(LD_matrix21)

# How many missing values we have
length(LD_matrix1[is.nan(LD_matrix1)]) 
length(LD_matrix9[is.nan(LD_matrix9)]) 
length(LD_matrix11[is.nan(LD_matrix11)]) 
length(LD_matrix21[is.nan(LD_matrix21)]) 

# Recode the missing values with a 0
LD_matrix1[is.nan(LD_matrix1)] = 0
LD_matrix9[is.nan(LD_matrix9)] = 0
LD_matrix11[is.nan(LD_matrix11)] = 0
LD_matrix21[is.nan(LD_matrix21)] = 0

#### LDna analysis ####

# Run LDna
ldna1 = LDnaRaw(LD_matrix1)
ldna9 = LDnaRaw(LD_matrix9)
ldna11 = LDnaRaw(LD_matrix11)
ldna21 = LDnaRaw(LD_matrix21)

# Clusters of LDna
clusters1 = extractClusters(ldna1, 20, rm.COCs = T, lambda.lim = 5, plot.graph = TRUE)
clusters9 = extractClusters(ldna9, 20, rm.COCs = T, lambda.lim = 5, plot.graph = TRUE)
clusters11 = extractClusters(ldna11, 20, rm.COCs = T, lambda.lim = 5, plot.graph = TRUE)
clusters21 = extractClusters(ldna21, 20, rm.COCs = T, lambda.lim = 5, plot.graph = TRUE)

# LDna summary of clusters
LDna_sum1 = summaryLDna(ldna = ldna1, LD_matrix1, clusters = clusters1)
LDna_sum9 = summaryLDna(ldna = ldna9, LD_matrix9, clusters = clusters9)
LDna_sum11 = summaryLDna(ldna = ldna11, LD_matrix11, clusters = clusters11)
LDna_sum21 = summaryLDna(ldna = ldna21, LD_matrix21, clusters = clusters21)

# Print the clusters
clusters1$clusters
clusters9$clusters
clusters11$clusters
clusters21$clusters

# Read gds file to perform PCA
RV1 = snpgdsPCA(file1)
RV9 = snpgdsPCA(file9)
RV11 = snpgdsPCA(file11)
RV21 = snpgdsPCA(file21)

# Number of clustered loci
clustered_loci1 = length(-1*as.numeric(unlist(clusters1[[1]])))
clustered_loci1
clustered_loci9 = length(-1*as.numeric(unlist(clusters9[[1]])))
clustered_loci9
clustered_loci11 = length(-1*as.numeric(unlist(clusters11[[1]])))
clustered_loci11
clustered_loci21 = length(-1*as.numeric(unlist(clusters21[[1]])))
clustered_loci21

#### Information to plot PCAs ####

# Set the working directory
setwd("~/Desktop/Master's Thesis/Results/14. LDna_Region/Clusters")

# Take the cluster of interest
cluster_interest1 = as.data.frame(clusters1$clusters$`498_0.57`)
cluster_interest9 = as.data.frame(clusters9$clusters$`405_0.72`)
cluster_interest11 = as.data.frame(clusters11$clusters$`127_0.48`)
cluster_interest21 = as.data.frame(clusters21$clusters$`172_0.75`)

# Change column name to SNP_ID
names(cluster_interest1)[1] = "SNP_ID"
names(cluster_interest9)[1] = "SNP_ID"
names(cluster_interest11)[1] = "SNP_ID"
names(cluster_interest21)[1] = "SNP_ID"

# Copy the column
cluster_interest1$CHR = cluster_interest1$SNP_ID
cluster_interest9$CHR = cluster_interest9$SNP_ID
cluster_interest11$CHR = cluster_interest11$SNP_ID
cluster_interest21$CHR = cluster_interest21$SNP_ID

# Get the chromosome number
cluster_interest1$CHR = gsub("_.*","", cluster_interest1$CHR)
cluster_interest1$CHR = gsub("group","", cluster_interest1$CHR)
cluster_interest9$CHR = gsub("_.*","", cluster_interest9$CHR)
cluster_interest9$CHR = gsub("group","", cluster_interest9$CHR)
cluster_interest11$CHR = gsub("_.*","", cluster_interest11$CHR)
cluster_interest11$CHR = gsub("group","", cluster_interest11$CHR)
cluster_interest21$CHR = gsub("_.*","", cluster_interest21$CHR)
cluster_interest21$CHR = gsub("group","", cluster_interest21$CHR)

# Change the roman number to numeric number
cluster_interest1$CHR = as.numeric(as.roman(cluster_interest1$CHR))
cluster_interest9$CHR = as.numeric(as.roman(cluster_interest9$CHR))
cluster_interest11$CHR = as.numeric(as.roman(cluster_interest11$CHR))
cluster_interest21$CHR = as.numeric(as.roman(cluster_interest21$CHR))

# Copy the column
cluster_interest1$BP = cluster_interest1$SNP_ID
cluster_interest9$BP = cluster_interest9$SNP_ID
cluster_interest11$BP = cluster_interest11$SNP_ID
cluster_interest21$BP = cluster_interest21$SNP_ID

# Get the location of the SNPs
cluster_interest1$BP = gsub(".*_","", cluster_interest1$BP)
cluster_interest9$BP = gsub(".*_","", cluster_interest9$BP)
cluster_interest11$BP = gsub(".*_","", cluster_interest11$BP)
cluster_interest21$BP = gsub(".*_","", cluster_interest21$BP)

# PCA information from cluster of interest
PCA1 = snpgdsPCA(file1, snp.id = as.character(clusters1$clusters$`498_0.57`))
PCA9 = snpgdsPCA(file9, snp.id = as.character(clusters9$clusters$`405_0.72`))
PCA11 = snpgdsPCA(file11, snp.id = as.character(clusters11$clusters$`127_0.48`))
PCA21 = snpgdsPCA(file21, snp.id = as.character(clusters21$clusters$`172_0.75`))

# Import data on population type and geographic location
pops_all = read.table("/Users/carlacollcosta/Desktop/Master's Thesis/Lists/Clasification_List_Ind_Pop_Type_Location_230_All(NotRepeatedFilesFiltered).txt", header=T, strip.white = T)

# Convert elements of the list in a data frame
samples1 = data.frame(PCA1$sample.id)
eigenvec1 = data.frame(PCA1$eigenvect)
names(eigenvec1)[1:ncol(eigenvec1)] = paste0("PC", 1:(ncol(eigenvec1)))
varprop1 = data.frame(PCA1$varprop)
samples9 = data.frame(PCA9$sample.id)
eigenvec9 = data.frame(PCA9$eigenvect)
names(eigenvec9)[1:ncol(eigenvec9)] = paste0("PC", 1:(ncol(eigenvec9)))
varprop9 = data.frame(PCA9$varprop)
samples11 = data.frame(PCA11$sample.id)
eigenvec11 = data.frame(PCA11$eigenvect)
names(eigenvec11)[1:ncol(eigenvec11)] = paste0("PC", 1:(ncol(eigenvec11)))
varprop11 = data.frame(PCA11$varprop)
samples21 = data.frame(PCA21$sample.id)
eigenvec21 = data.frame(PCA21$eigenvect)
names(eigenvec21)[1:ncol(eigenvec21)] = paste0("PC", 1:(ncol(eigenvec21)))
varprop21 = data.frame(PCA21$varprop)

# Set the final data frame to plot the PCAs
finalPCA1 = merge(samples1, pops_all, by.x="PCA1.sample.id", by.y="Individual")
names(finalPCA1)[1] = "Individual"
finalPCA1 = data.frame(finalPCA1, eigenvec1)
finalPCA9 = merge(samples9, pops_all, by.x="PCA9.sample.id", by.y="Individual")
names(finalPCA9)[1] = "Individual"
finalPCA9 = data.frame(finalPCA9, eigenvec9)
finalPCA11 = merge(samples11, pops_all, by.x="PCA11.sample.id", by.y="Individual")
names(finalPCA11)[1] = "Individual"
finalPCA11 = data.frame(finalPCA11, eigenvec11)
finalPCA21 = merge(samples21, pops_all, by.x="PCA21.sample.id", by.y="Individual")
names(finalPCA21)[1] = "Individual"
finalPCA21 = data.frame(finalPCA21, eigenvec21)

#### Information to plot Heterozygosity vs PC1 ####

# Heterozygosity of the region
het1 = read.table("~/Desktop/Master's Thesis/Results/14. LDna_Region/0_data/HeterozygosityClusters/Genotypes_Filtered_2.3_Chromosome1_Cluster.het", header = T)
het9 = read.table("~/Desktop/Master's Thesis/Results/14. LDna_Region/0_data/HeterozygosityClusters/Genotypes_Filtered_2.3_Chromosome9_Cluster.het", header = T)
het11 = read.table("~/Desktop/Master's Thesis/Results/14. LDna_Region/0_data/HeterozygosityClusters/Genotypes_Filtered_2.3_Chromosome11_Cluster.het", header = T)
het21 = read.table("~/Desktop/Master's Thesis/Results/14. LDna_Region/0_data/HeterozygosityClusters/Genotypes_Filtered_2.3_Chromosome21_Cluster.het", header = T)

# Calculate observed heterozygosity
het1$O.HET = (het1$N_SITES-het1$O.HOM.)/het1$N_SITES
het9$O.HET = (het9$N_SITES-het9$O.HOM.)/het9$N_SITES
het11$O.HET = (het11$N_SITES-het11$O.HOM.)/het11$N_SITES
het21$O.HET = (het21$N_SITES-het21$O.HOM.)/het21$N_SITES

# Merge dataframes
finalPCA1 = merge(finalPCA1, het1, by.x="Individual", by.y="INDV")
finalPCA9 = merge(finalPCA9, het9, by.x="Individual", by.y="INDV")
finalPCA11 = merge(finalPCA11, het11, by.x="Individual", by.y="INDV")
finalPCA21 = merge(finalPCA21, het21, by.x="Individual", by.y="INDV")


#### Plots ####

# PCAs

plotPCA1 = ggplot(data=finalPCA1, aes(x=PC1, y=PC2, colour=Concrete_Location)) +
  geom_point(aes(shape=Type), size=5, alpha=0.8) +
  ggtitle("Chromosome 1 - Cluster 498_0.57") +
  xlab(paste("PC1",round(varprop1[1,1]*100,2),"%")) +
  ylab(paste("PC2",round(varprop1[2,1]*100,2),"%")) +
  scale_shape_manual(values=c("Freshwater"=16, "Marine"=1)) +
  scale_colour_manual(values=c("Adriatic"="orange", "Alaska"="darkolivegreen4", "North_Scandinavia"="royalblue3", "East_Rusia"="tomato4", "English_Channel"="yellowgreen", "Finnish_Golf"="purple3", "Iberian_Peninsula"="red", "North_Sea"="turquoise4"),
                      labels=c("Adriatic"="Adriatic Sea", "Alaska"="Alaska", "North_Scandinavia"="North Scandinavia", "East_Rusia"="East Russia", "English_Channel"="English Channel", "Finnish_Golf"="Finnish Golf", "Iberian_Peninsula"="Iberian Peninsula", "North_Sea"="North Sea")) +
  guides(colour=guide_legend(order=1, title="Geographic Region", override.aes=list(size=3, alpha=1)), shape=guide_legend(order=2, title="Habitat Type",override.aes=list(size=3, alpha=1))) +
  theme(plot.title = element_text(hjust = 0.5))

plotPCA9 = ggplot(data=finalPCA9, aes(x=PC1, y=PC2, colour=Concrete_Location)) +
  geom_point(aes(shape=Type), size=5, alpha=0.8) +
  ggtitle("Chromosome 9 - Cluster 405_0.72") +
  xlab(paste("PC1",round(varprop9[1,1]*100,2),"%")) +
  ylab(paste("PC2",round(varprop9[2,1]*100,2),"%")) +
  scale_shape_manual(values=c("Freshwater"=16, "Marine"=1)) +
  scale_colour_manual(values=c("Adriatic"="orange", "Alaska"="darkolivegreen4", "North_Scandinavia"="royalblue3", "East_Rusia"="tomato4", "English_Channel"="yellowgreen", "Finnish_Golf"="purple3", "Iberian_Peninsula"="red", "North_Sea"="turquoise4"),
                      labels=c("Adriatic"="Adriatic Sea", "Alaska"="Alaska", "North_Scandinavia"="North Scandinavia", "East_Rusia"="East Russia", "English_Channel"="English Channel", "Finnish_Golf"="Finnish Golf", "Iberian_Peninsula"="Iberian Peninsula", "North_Sea"="North Sea")) +
  guides(colour=guide_legend(order=1, title="Geographic Region", override.aes=list(size=3, alpha=1)), shape=guide_legend(order=2, title="Habitat Type",override.aes=list(size=3, alpha=1))) +
  theme(plot.title = element_text(hjust = 0.5))

plotPCA11 = ggplot(data=finalPCA11, aes(x=PC1, y=PC2, colour=Concrete_Location)) +
  geom_point(aes(shape=Type), size=5, alpha=0.8) +
  ggtitle("Chromosome 11 - Cluster 127_0.48") +
  xlab(paste("PC1",round(varprop11[1,1]*100,2),"%")) +
  ylab(paste("PC2",round(varprop11[2,1]*100,2),"%")) +
  scale_shape_manual(values=c("Freshwater"=16, "Marine"=1)) +
  scale_colour_manual(values=c("Adriatic"="orange", "Alaska"="darkolivegreen4", "North_Scandinavia"="royalblue3", "East_Rusia"="tomato4", "English_Channel"="yellowgreen", "Finnish_Golf"="purple3", "Iberian_Peninsula"="red", "North_Sea"="turquoise4"),
                      labels=c("Adriatic"="Adriatic Sea", "Alaska"="Alaska", "North_Scandinavia"="North Scandinavia", "East_Rusia"="East Russia", "English_Channel"="English Channel", "Finnish_Golf"="Finnish Golf", "Iberian_Peninsula"="Iberian Peninsula", "North_Sea"="North Sea")) +
  guides(colour=guide_legend(order=1, title="Geographic Region", override.aes=list(size=3, alpha=1)), shape=guide_legend(order=2, title="Habitat Type",override.aes=list(size=3, alpha=1))) +
  theme(plot.title = element_text(hjust = 0.5))

plotPCA21 = ggplot(data=finalPCA21, aes(x=PC1, y=PC2, colour=Concrete_Location)) +
  geom_point(aes(shape=Type), size=5, alpha=0.8) +
  ggtitle("Chromosome 21 - Cluster 172_0.75") +
  xlab(paste("PC1",round(varprop21[1,1]*100,2),"%")) +
  ylab(paste("PC2",round(varprop21[2,1]*100,2),"%")) +
  scale_shape_manual(values=c("Freshwater"=16, "Marine"=1)) +
  scale_colour_manual(values=c("Adriatic"="orange", "Alaska"="darkolivegreen4", "North_Scandinavia"="royalblue3", "East_Rusia"="tomato4", "English_Channel"="yellowgreen", "Finnish_Golf"="purple3", "Iberian_Peninsula"="red", "North_Sea"="turquoise4"),
                      labels=c("Adriatic"="Adriatic Sea", "Alaska"="Alaska", "North_Scandinavia"="North Scandinavia", "East_Rusia"="East Russia", "English_Channel"="English Channel", "Finnish_Golf"="Finnish Golf", "Iberian_Peninsula"="Iberian Peninsula", "North_Sea"="North Sea")) +
  guides(colour=guide_legend(order=1, title="Geographic Region", override.aes=list(size=3, alpha=1)), shape=guide_legend(order=2, title="Habitat Type",override.aes=list(size=3, alpha=1))) +
  theme(plot.title = element_text(hjust = 0.5))

# Heterozygosity vs PC1

plothet1 = ggplot(data=finalPCA1, aes(x=PC1, y=`O.HET`, colour=Concrete_Location)) +
  geom_point(aes(shape=Type), size=5, alpha=0.8) +
  ggtitle("Chromosome 1 - Cluster 498_0.57") +
  xlab(paste("PC1",round(varprop1[1,1]*100,2),"%")) +
  ylab(expression(H[O])) +
  scale_shape_manual(values=c("Freshwater"=16, "Marine"=1)) +
  scale_colour_manual(values=c("Adriatic"="orange", "Alaska"="darkolivegreen4", "North_Scandinavia"="royalblue3", "East_Rusia"="tomato4", "English_Channel"="yellowgreen", "Finnish_Golf"="purple3", "Iberian_Peninsula"="red", "North_Sea"="turquoise4"),
                      labels=c("Adriatic"="Adriatic Sea", "Alaska"="Alaska", "North_Scandinavia"="North Scandinavia", "East_Rusia"="East Russia", "English_Channel"="English Channel", "Finnish_Golf"="Finnish Golf", "Iberian_Peninsula"="Iberian Peninsula", "North_Sea"="North Sea")) +
  guides(colour=guide_legend(order=1, title="Geographic Region", override.aes=list(size=3, alpha=1)), shape=guide_legend(order=2, title="Habitat Type",override.aes=list(size=3, alpha=1))) +
  theme(plot.title = element_text(hjust = 0.5))

plothet9 = ggplot(data=finalPCA9, aes(x=PC1, y=`O.HET`, colour=Concrete_Location)) +
  geom_point(aes(shape=Type), size=5, alpha=0.8) +
  ggtitle("Chromosome 9 - Cluster 405_0.72") +
  xlab(paste("PC1",round(varprop9[1,1]*100,2),"%")) +
  ylab(expression(H[O])) +
  scale_shape_manual(values=c("Freshwater"=16, "Marine"=1)) +
  scale_colour_manual(values=c("Adriatic"="orange", "Alaska"="darkolivegreen4", "North_Scandinavia"="royalblue3", "East_Rusia"="tomato4", "English_Channel"="yellowgreen", "Finnish_Golf"="purple3", "Iberian_Peninsula"="red", "North_Sea"="turquoise4"),
                      labels=c("Adriatic"="Adriatic Sea", "Alaska"="Alaska", "North_Scandinavia"="North Scandinavia", "East_Rusia"="East Russia", "English_Channel"="English Channel", "Finnish_Golf"="Finnish Golf", "Iberian_Peninsula"="Iberian Peninsula", "North_Sea"="North Sea")) +
  guides(colour=guide_legend(order=1, title="Geographic Region", override.aes=list(size=3, alpha=1)), shape=guide_legend(order=2, title="Habitat Type",override.aes=list(size=3, alpha=1))) +
  theme(plot.title = element_text(hjust = 0.5))

plothet11 = ggplot(data=finalPCA11, aes(x=PC1, y=`O.HET`, colour=Concrete_Location)) +
  geom_point(aes(shape=Type), size=5, alpha=0.8) +
  ggtitle("Chromosome 11 - Cluster 127_0.48") +
  xlab(paste("PC1",round(varprop11[1,1]*100,2),"%")) +
  ylab(expression(H[O])) +
  scale_shape_manual(values=c("Freshwater"=16, "Marine"=1)) +
  scale_colour_manual(values=c("Adriatic"="orange", "Alaska"="darkolivegreen4", "North_Scandinavia"="royalblue3", "East_Rusia"="tomato4", "English_Channel"="yellowgreen", "Finnish_Golf"="purple3", "Iberian_Peninsula"="red", "North_Sea"="turquoise4"),
                      labels=c("Adriatic"="Adriatic Sea", "Alaska"="Alaska", "North_Scandinavia"="North Scandinavia", "East_Rusia"="East Russia", "English_Channel"="English Channel", "Finnish_Golf"="Finnish Golf", "Iberian_Peninsula"="Iberian Peninsula", "North_Sea"="North Sea")) +
  guides(colour=guide_legend(order=1, title="Geographic Region", override.aes=list(size=3, alpha=1)), shape=guide_legend(order=2, title="Habitat Type",override.aes=list(size=3, alpha=1))) +
  theme(plot.title = element_text(hjust = 0.5))

plothet21 = ggplot(data=finalPCA21, aes(x=PC1, y=`O.HET`, colour=Concrete_Location)) +
  geom_point(aes(shape=Type), size=5, alpha=0.8) +
  ggtitle("Chromosome 21 - Cluster 172_0.75") +
  xlab(paste("PC1",round(varprop21[1,1]*100,2),"%")) +
  ylab(expression(H[O])) +
  scale_shape_manual(values=c("Freshwater"=16, "Marine"=1)) +
  scale_colour_manual(values=c("Adriatic"="orange", "Alaska"="darkolivegreen4", "North_Scandinavia"="royalblue3", "East_Rusia"="tomato4", "English_Channel"="yellowgreen", "Finnish_Golf"="purple3", "Iberian_Peninsula"="red", "North_Sea"="turquoise4"),
                      labels=c("Adriatic"="Adriatic Sea", "Alaska"="Alaska", "North_Scandinavia"="North Scandinavia", "East_Rusia"="East Russia", "English_Channel"="English Channel", "Finnish_Golf"="Finnish Golf", "Iberian_Peninsula"="Iberian Peninsula", "North_Sea"="North Sea")) +
  guides(colour=guide_legend(order=1, title="Geographic Region", override.aes=list(size=3, alpha=1)), shape=guide_legend(order=2, title="Habitat Type",override.aes=list(size=3, alpha=1))) +
  theme(plot.title = element_text(hjust = 0.5))

# Save plots

setwd("~/Desktop/Master's Thesis/Results/14. LDna_Region")

png("Chromosomes 1, 11, 21.png", units="in", width=15, height=8, res=900)
ggarrange(plotPCA1,
          plotPCA11,
          plotPCA21,
          plothet1,
          plothet11,
          plothet21,
          ncol=3, nrow=2, labels="auto", hjust=-1, align="hv", common.legend = TRUE, legend="right")
dev.off()

png("Chromosome 9.png", units="in", width=6.5, height=8, res=900)
ggarrange(plotPCA9,
          plothet9,
          ncol=1, nrow=2, labels="auto", hjust=-1, align="hv", common.legend = TRUE, legend="right")
dev.off()

