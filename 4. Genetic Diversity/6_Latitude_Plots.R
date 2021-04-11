#################################################################
##             PLOT GENETIC DIVERSITY VS. LATITUDE             ##
#################################################################

## Carla Coll Costa, C3. January 2021.
## This script plots genetic diversity statistics (nucleotide diversity and mean Froh per population) against
## latitude and saves the plots automatically. It expects .txt files with the statistics for each population
## and with the latitude coordinates.

# Upload packages
library(readxl)
library(ggplot2)
library(scales)
library(ggpubr)

# Set working directory
setwd("<Working_Directory>")

# Get the data
pi = read.table("<File_with_Nucleotide_Diversity_of_each_Population>", header=TRUE)
froh = read.table("<File_with_Froh_for_each_Individual>", header=TRUE)
latitude = data.frame(read_excel("<File_with_Information_on_the_Latitude_Coordinates>"))
pop_list = read.table("<File_with_Additional_Information_on_the_Populations>", header=TRUE)

# Order by population name
pi = pi[order(pi$POP),]
froh = froh[order(froh$Population),]
latitude = latitude[order(latitude$POP_ID),]
pop_list = pop_list[order(pop_list$Population),]

# Calculate mean Froh per population
froh = aggregate(Froh_genome~Population, data=froh, FUN=mean)
froh = data.frame(froh$Population, froh$Froh_genome)

# Fuse data frames for pi and froh
data_stats = merge(latitude, pi, by.x="POP_ID", by.y="POP")
data_stats = merge(data_stats, froh, by.x="POP_ID", by.y="froh.Population")
data_stats = merge(data_stats, pop_list, by.x="POP_ID", by.y="Population")

# Delete columns that are not necessary
data_stats$N.x = NULL
data_stats$N.y = NULL
data_stats$nsites = NULL
data_stats$S = NULL
data_stats$theta = NULL
data_stats$TajimaD = NULL
data_stats$Type.x = NULL
data_stats$Individual = NULL

# Remove duplicate rows
data_stats = unique(data_stats)

# Order the geographic location in alphabetical order
data_stats$Concrete_Location = factor(data_stats$Concrete_Location, levels = c("Adriatic", "Alaska", "East_Rusia", "English_Channel", "North_Scandinavia", "Finnish_Golf", "Iberian_Peninsula", "North_Sea"))

# Plots
png("Nucleotide Diversity vs Latitude.png", units="in", width=8, height=6, res=900)
ggplot(data=data_stats, aes(x=GPS_N, y=pi, color=Concrete_Location)) +
  geom_point(size=6, aes(shape=Type.y), stroke=1, alpha=0.8) +
  labs(x="Latitude (GPS North)", y="Nucleotide Diversity (π)", color="Geographic Region") +
  guides(colour=guide_legend(order=1, title="Geographic Region", override.aes=list(size=3, alpha=1)), shape=guide_legend(order=2, title="Habitat Type",override.aes=list(size=3, alpha=1))) +
  scale_shape_manual(values=c("Freshwater"=16, "Marine"=1)) +
  scale_color_manual(values=c("Adriatic"="orange", "Alaska"="darkolivegreen4", "North_Scandinavia"="royalblue3", "East_Rusia"="tomato4", "English_Channel"="chartreuse3", "Finnish_Golf"="purple3", "Iberian_Peninsula"="red", "North_Sea"="turquoise4"),
                     labels=c("Adriatic"="Adriatic Sea", "Alaska"="Alaska", "North_Scandinavia"="Fennoscandia", "East_Rusia"="East Russia", "English_Channel"="English Channel", "Finnish_Golf"="Gulf of Finland", "Iberian_Peninsula"="Iberian Peninsula", "North_Sea"="North Sea")) +
  theme_bw()
dev.off()

png("Froh vs Latitude.png", units="in", width=8, height=6, res=900)
ggplot(data=data_stats, aes(x=GPS_N, y=froh.Froh_genome, color=Concrete_Location)) +
  geom_point(size=6, aes(shape=Type.y), stroke=1, alpha=0.8) +
  labs(x="Latitude (GPS North)", y=expression(F[ROH]), color="Geographic Region") +
  guides(colour=guide_legend(order=1, title="Geographic Region", override.aes=list(size=3, alpha=1)), shape=guide_legend(order=2, title="Habitat Type",override.aes=list(size=3, alpha=1))) +
  scale_shape_manual(values=c("Freshwater"=16, "Marine"=1)) +
  scale_color_manual(values=c("Adriatic"="orange", "Alaska"="darkolivegreen4", "North_Scandinavia"="royalblue3", "East_Rusia"="tomato4", "English_Channel"="chartreuse3", "Finnish_Golf"="purple3", "Iberian_Peninsula"="red", "North_Sea"="turquoise4"),
                     labels=c("Adriatic"="Adriatic Sea", "Alaska"="Alaska", "East_Rusia"="East Russia", "English_Channel"="English Channel", "Finnish_Golf"="Gulf of Finland", "Iberian_Peninsula"="Iberian Peninsula", "North_Sea"="North Sea", "North_Scandinavia"="Fennoscandia")) +
  theme_bw()
dev.off()

# Fused plots
nucdiv_plot = ggplot(data=data_stats, aes(x=GPS_N, y=pi, color=Concrete_Location)) +
  geom_point(size=6, aes(shape=Type.y), stroke=1, alpha=0.8) +
  labs(x="Latitude (GPS North)", y="Nucleotide Diversity (π)", color="Geographic Region") +
  guides(colour=guide_legend(order=1, title="Geographic Region", override.aes=list(size=3, alpha=1)), shape=guide_legend(order=2, title="Habitat Type",override.aes=list(size=3, alpha=1))) +
  scale_shape_manual(values=c("Freshwater"=16, "Marine"=1)) +
  scale_color_manual(values=c("Adriatic"="orange", "Alaska"="darkolivegreen4", "North_Scandinavia"="royalblue3", "East_Rusia"="tomato4", "English_Channel"="chartreuse3", "Finnish_Golf"="purple3", "Iberian_Peninsula"="red", "North_Sea"="turquoise4"),
                     labels=c("Adriatic"="Adriatic Sea", "Alaska"="Alaska", "North_Scandinavia"="Fennoscandia", "East_Rusia"="East Russia", "English_Channel"="English Channel", "Finnish_Golf"="Gulf of Finland", "Iberian_Peninsula"="Iberian Peninsula", "North_Sea"="North Sea")) +
  theme_bw()

froh_plot = ggplot(data=data_stats, aes(x=GPS_N, y=froh.Froh_genome, color=Concrete_Location)) +
  geom_point(size=6, aes(shape=Type.y), stroke=1, alpha=0.8) +
  labs(x="Latitude (GPS North)", y=expression(F[ROH]), color="Geographic Region") +
  guides(colour=guide_legend(order=1, title="Geographic Region", override.aes=list(size=3, alpha=1)), shape=guide_legend(order=2, title="Habitat Type",override.aes=list(size=3, alpha=1))) +
  scale_shape_manual(values=c("Freshwater"=16, "Marine"=1)) +
  scale_color_manual(values=c("Adriatic"="orange", "Alaska"="darkolivegreen4", "North_Scandinavia"="royalblue3", "East_Rusia"="tomato4", "English_Channel"="chartreuse3", "Finnish_Golf"="purple3", "Iberian_Peninsula"="red", "North_Sea"="turquoise4"),
                     labels=c("Adriatic"="Adriatic Sea", "Alaska"="Alaska", "East_Rusia"="East Russia", "English_Channel"="English Channel", "Finnish_Golf"="Gulf of Finland", "Iberian_Peninsula"="Iberian Peninsula", "North_Sea"="North Sea", "North_Scandinavia"="Fennoscandia")) +
  theme_bw()

png("Statistics vs Latitude.png", units="in", width=8, height=4, res=900)
ggarrange(nucdiv_plot,
          froh_plot,
          ncol=2, nrow=1, labels="auto", hjust=-1, align="hv", common.legend = TRUE, legend="bottom")
dev.off()
