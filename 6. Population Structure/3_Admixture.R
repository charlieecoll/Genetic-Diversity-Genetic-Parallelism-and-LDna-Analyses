#################################################################
##                       ADMIXTURE PLOTS                       ##
#################################################################

## Carla Coll Costa, C3. December 2020.
## This script takes data from ADMIXTURE and plots the cross-validation
## error for different Ks and the ancestry proportions for the higher K.
## Plots are saved automatically.


# Upload the libraries required
library(ggplot2)
library(reshape2)

# Set the working directory
setwd("<Working_Directory>")

# Get the file with the cross-validation errors for each K
cverror = read.table("<File_with_CV_Information>")
cverror$K = as.numeric(gsub("[^[:digit:]]+", "", cverror$V1))

# Print the cv values by increasing order
cverror[with(cverror, order(cverror$V4)),]

# Save the cross-validation plot
pdf("Cross-Validation.pdf")
ggplot(data=cverror, aes(x=K, y=V4)) +
  geom_line() +
  geom_point() +
  labs(y="Cross-validation Error")
dev.off()  

## Now we plot the admixture results for different K

# Bam files with the names of the individuals of each of the populations
individuals = read.table("<Input_File.fam>") 
names = individuals$V1

# txt file with the population names and other information
population = read.table("<File_with_Population_Information>", header=TRUE)
population = population[order(population$Individual),]

# Admixture file with ancestry values
values5 = read.table("<Input_File_K5.Q>") 
values6 = read.table("<Input_File_K6.Q>") 
values7 = read.table("<Input_File_K7.Q>") 
values8 = read.table("<Input_File_K8.Q>") 
values9 = read.table("<Input_File_K9.Q>") 
values10 = read.table("Input_File_K10.Q") 
values11 = read.table("Input_File_K11.Q") 
values12 = read.table("Input_File_K12.Q") 
values13 = read.table("Input_File_K13.Q") 
values14 = read.table("Input_File_K14.Q") 
values15 = read.table("Input_File_K15.Q") 

# Rename the rows of the ancestry dataframe with the individual name
rownames(values5) = names
rownames(values6) = names
rownames(values7) = names
rownames(values8) = names
rownames(values9) = names
rownames(values10) = names
rownames(values11) = names
rownames(values12) = names
rownames(values13) = names
rownames(values14) = names
rownames(values15) = names

# Add the population names to the matrix of ancestry
ancestry5 = data.frame(population,values5)
ancestry6 = data.frame(population,values6)
ancestry7 = data.frame(population,values7)
ancestry8 = data.frame(population,values8)
ancestry9 = data.frame(population,values9)
ancestry10 = data.frame(population,values10)
ancestry11 = data.frame(population,values11)
ancestry12 = data.frame(population,values12)
ancestry13 = data.frame(population,values13)
ancestry14 = data.frame(population,values14)
ancestry15 = data.frame(population,values15)

# We melt the matrix to transform it to a dataframe that can be plotted with ggplot2
admixture5 = melt(ancestry5)
admixture6 = melt(ancestry6)
admixture7 = melt(ancestry7)
admixture8 = melt(ancestry8)
admixture9 = melt(ancestry9)
admixture10 = melt(ancestry10)
admixture11 = melt(ancestry11)
admixture12 = melt(ancestry12)
admixture13 = melt(ancestry13)
admixture14 = melt(ancestry14)
admixture15 = melt(ancestry15)

# Vector of colors
colors = c("#7FC97F", "#BEAED4", "#FDC086", "#FFFF99", "#386CB0", "#F0027F", 
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

# Plot the admixture plots
png("Admixture_K5.png", units="in", width=18, height=8, res=900)
ggplot(data=admixture5, aes(x=Individual, y=value, fill=variable)) +
  geom_bar(stat="identity") +
  xlab("Population") +
  ylab("Ancestry") + 
  scale_fill_manual(values=c(colors)) +
  facet_grid(~Population, scales="free_x", switch="x") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.background = element_blank(),
        legend.position = "none",
        strip.placement = "outside", 
        strip.background = element_rect(fill = "white"),
        strip.text.x=element_text(angle=90, size=10))
dev.off()

png("Admixture_K6.png", units="in", width=18, height=8, res=900)
ggplot(data=admixture6, aes(x=Individual, y=value, fill=variable)) +
  geom_bar(stat="identity") +
  xlab("Population") +
  ylab("Ancestry") + 
  scale_fill_manual(values=c(colors)) +
  facet_grid(~Population, scales="free_x", switch="x") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.background = element_blank(),
        legend.position = "none",
        strip.placement = "outside", 
        strip.background = element_rect(fill = "white"),
        strip.text.x=element_text(angle=90, size=10))
dev.off()

png("Admixture_K7.png", units="in", width=18, height=8, res=900)
ggplot(data=admixture7, aes(x=Individual, y=value, fill=variable)) +
  geom_bar(stat="identity") +
  xlab("Population") +
  ylab("Ancestry") + 
  scale_fill_manual(values=c(colors)) +
  facet_grid(~Population, scales="free_x", switch="x") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.background = element_blank(),
        legend.position = "none",
        strip.placement = "outside", 
        strip.background = element_rect(fill = "white"),
        strip.text.x=element_text(angle=90, size=10))
dev.off()

png("Admixture_K8.png", units="in", width=18, height=8, res=900)
ggplot(data=admixture8, aes(x=Individual, y=value, fill=variable)) +
  geom_bar(stat="identity") +
  xlab("Population") +
  ylab("Ancestry") + 
  scale_fill_manual(values=c(colors)) +
  facet_grid(~Population, scales="free_x", switch="x") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.background = element_blank(),
        legend.position = "none",
        strip.placement = "outside", 
        strip.background = element_rect(fill = "white"),
        strip.text.x=element_text(angle=90, size=10))
dev.off()

png("Admixture_K9.png", units="in", width=18, height=8, res=900)
ggplot(data=admixture9, aes(x=Individual, y=value, fill=variable)) +
  geom_bar(stat="identity") +
  xlab("Population") +
  ylab("Ancestry") + 
  scale_fill_manual(values=c(colors)) +
  facet_grid(~Population, scales="free_x", switch="x") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.background = element_blank(),
        legend.position = "none",
        strip.placement = "outside", 
        strip.background = element_rect(fill = "white"),
        strip.text.x=element_text(angle=90, size=10))
dev.off()

png("Admixture_K10.png", units="in", width=18, height=8, res=900)
ggplot(data=admixture10, aes(x=Individual, y=value, fill=variable)) +
  geom_bar(stat="identity") +
  xlab("Population") +
  ylab("Ancestry") + 
  scale_fill_manual(values=c(colors)) +
  facet_grid(~Population, scales="free_x", switch="x") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.background = element_blank(),
        legend.position = "none",
        strip.placement = "outside", 
        strip.background = element_rect(fill = "white"),
        strip.text.x=element_text(angle=90, size=10))
dev.off()

png("Admixture_K11.png", units="in", width=18, height=8, res=900)
ggplot(data=admixture11, aes(x=Individual, y=value, fill=variable)) +
  geom_bar(stat="identity") +
  xlab("Population") +
  ylab("Ancestry") + 
  scale_fill_manual(values=c(colors)) +
  facet_grid(~Population, scales="free_x", switch="x") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.background = element_blank(),
        legend.position = "none",
        strip.placement = "outside", 
        strip.background = element_rect(fill = "white"),
        strip.text.x=element_text(angle=90, size=10))
dev.off()

png("Admixture_K12.png", units="in", width=18, height=8, res=900)
ggplot(data=admixture12, aes(x=Individual, y=value, fill=variable)) +
  geom_bar(stat="identity") +
  xlab("Population") +
  ylab("Ancestry") + 
  scale_fill_manual(values=c(colors)) +
  facet_grid(~Population, scales="free_x", switch="x") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.background = element_blank(),
        legend.position = "none",
        strip.placement = "outside", 
        strip.background = element_rect(fill = "white"),
        strip.text.x=element_text(angle=90, size=10))
dev.off()

png("Admixture_K13.png", units="in", width=18, height=8, res=900)
ggplot(data=admixture13, aes(x=Individual, y=value, fill=variable)) +
  geom_bar(stat="identity") +
  xlab("Population") +
  ylab("Ancestry") + 
  scale_fill_manual(values=c(colors)) +
  facet_grid(~Population, scales="free_x", switch="x") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.background = element_blank(),
        legend.position = "none",
        strip.placement = "outside", 
        strip.background = element_rect(fill = "white"),
        strip.text.x=element_text(angle=90, size=10))
dev.off()

png("Admixture_K14.png", units="in", width=18, height=8, res=900)
ggplot(data=admixture14, aes(x=Individual, y=value, fill=variable)) +
  geom_bar(stat="identity") +
  xlab("Population") +
  ylab("Ancestry") + 
  scale_fill_manual(values=c(colors)) +
  facet_grid(~Population, scales="free_x", switch="x") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.background = element_blank(),
        legend.position = "none",
        strip.placement = "outside", 
        strip.background = element_rect(fill = "white"),
        strip.text.x=element_text(angle=90, size=10))
dev.off()

png("Admixture_K15.png", units="in", width=18, height=8, res=900)
ggplot(data=admixture15, aes(x=Individual, y=value, fill=variable)) +
  geom_bar(stat="identity") +
  xlab("Population") +
  ylab("Ancestry") + 
  scale_fill_manual(values=c(colors)) +
  facet_grid(~Population, scales="free_x", switch="x") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.background = element_blank(),
        legend.position = "none",
        strip.placement = "outside", 
        strip.background = element_rect(fill = "white"),
        strip.text.x=element_text(angle=90, size=10))
dev.off()
