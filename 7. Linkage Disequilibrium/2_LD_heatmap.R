##################################################################
##                          LD heatmap                          ##
##################################################################


## Carla Coll Costa, C3. December 2020.
## This script plots LD heatmaps and saves them automatically. It takes as input 012 files
## generated with vcftools and with the -1 values changed by NA.

# Upload the required libraries
library(gdsfmt)
library(LDna)
library(SNPRelate)
library(data.table)
library(ggplot2)
library(tm)
library(reshape2)
options(expressions = 500000)

# Set the working directory
setwd("~/Desktop/Master's Thesis/Results/13. LDheatmap_Region/All")

# Read the file with the genotypes (row: SNPs, column: individual)
genotypes = t(read.table("~/Desktop/Master's Thesis/Results/13. LDheatmap_Region/0_data/All/Genotypes_Filtered_2.3_Chromosome9_Region.012",header=FALSE))
# Delete first row
genotypes = genotypes[-1,]
# Save the genotypes as a matrix
genotypes = as.matrix(genotypes)

# Read the file with the genotype positions
positions = read.table("~/Desktop/Master's Thesis/Results/13. LDheatmap_Region/0_data/All/Genotypes_Filtered_2.3_Chromosome9_Region.012.pos", header = F)

# Create a new column with an ID
positions$LocID = paste0(positions$V1,"_",positions$V2)

# Create a vector as large as the number of genotypes
snps = 1:nrow(genotypes)

# Create a new column as large as the snps
positions$snps = positions$LocID

# Read the file with the individuals
samples = read.table("~/Desktop/Master's Thesis/Results/13. LDheatmap_Region/0_data/All/Genotypes_Filtered_2.3_Chromosome9_Region.012.indv", header = F, stringsAsFactors = F)

# Create a new column with the populations
samples$population = gsub("[0-9]+","", samples$V1)
samples$population = gsub("\\.","", samples$population)
samples$population = gsub(".{1}$", "", samples$population)

# Take the file with the populations
pop_labs = read.table("~/Desktop/Master's Thesis/Lists/Clasification_List_Ind_Pop_Type_Location_230_All(NotRepeatedFilesFiltered).txt", header = T)

# Order the files by the individuals
pop_labs = pop_labs[with(pop_labs, order(pop_labs$Individual)),]

# Create .gds file
snpgdsCreateGeno("Chromosome9.gds", genmat = genotypes, sample.id = pop_labs$Individual, snp.id = positions$LocID, snpfirstdim=TRUE) ### sample ID should be 1: how many individuals you have

# Save gds file
file = snpgdsOpen("Chromosome9.gds")

# Linkage Disequilibrium 
matrixx = snpgdsLDMat(file, method="r",slide=0, verbose=TRUE)
LD_matrix = (matrixx$LD)^2
LD_matrix[upper.tri(LD_matrix, diag = TRUE)] <- NA

# Change row names and column names
rownames(LD_matrix) = positions$LocID
colnames(LD_matrix) = positions$LocID

# Melt the ld matrix
LD_matrix_melted = melt(LD_matrix)

# Plot the LD matrix
png("LD Matrix Chromosome 9.png", units="in", width=10, height=10, res=600)
ggplot(data = LD_matrix_melted, aes(Var1, Var2)) +
  geom_tile(aes(fill=value)) +
  scale_fill_gradientn(colors=c("blue", "yellow", "red"), na.value = "white") +
  ggtitle("LD Matrix Chromosome 9")  +
  coord_fixed() + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())
dev.off()

