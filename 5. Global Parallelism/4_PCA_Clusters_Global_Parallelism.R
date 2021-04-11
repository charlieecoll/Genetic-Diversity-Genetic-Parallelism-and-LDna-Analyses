##################################################################
##                   PCA PARALLELISM CLUSTERS                   ##
##################################################################

## Carla Coll Costa, C3. February 2021.
## This script takes information of covariance matrices calculated with PCAngsd
## and plots PCAs and boxplots of populations against the PC1 axis. The plots are
## saved automatically. It takes several covariance matrices and plots information on
## each of those in a combined plot. Importantly, it also expects information on each 
## of the populations to colour the plots according to different characteristics: habitat
## type and geographic region.

#### Prepare the data ####

# Load the required libraries
library(sf)
library(RcppCNPy)
library(scales)
library(ggplot2)
library(gridExtra)
library(ggpubr)
library(patchwork)

# Set working directory
setwd("<Working_Directory>")

# Import covariance matrices
cov5 = as.matrix(read.table("<Input_Covariance_Matrix>"))
cov6 = as.matrix(read.table("<Input_Covariance_Matrix>"))
cov10 = as.matrix(read.table("<Input_Covariance_Matrix>"))
cov11 = as.matrix(read.table("<Input_Covariance_Matrix>"))
cov12 = as.matrix(read.table("<Input_Covariance_Matrix>"))
cov13 = as.matrix(read.table("<Input_Covariance_Matrix>"))
cov16 = as.matrix(read.table("<Input_Covariance_Matrix>"))
cov18 = as.matrix(read.table("<Input_Covariance_Matrix>"))
cov20 = as.matrix(read.table("<Input_Covariance_Matrix>"))
cov22 = as.matrix(read.table("<Input_Covariance_Matrix>"))
cov25 = as.matrix(read.table("<Input_Covariance_Matrix>"))
cov27 = as.matrix(read.table("<Input_Covariance_Matrix>"))

# Import individuals in the order of the bam files when calculating beagle file
individuals = read.table("<File_with_Individual_Names>")

# Calculate eigen data
eigen_data5 = eigen(cov5)
eigen_data6 = eigen(cov6)
eigen_data10 = eigen(cov10)
eigen_data11 = eigen(cov11)
eigen_data12 = eigen(cov12)
eigen_data13 = eigen(cov13)
eigen_data16 = eigen(cov16)
eigen_data18 = eigen(cov18)
eigen_data20 = eigen(cov20)
eigen_data22 = eigen(cov22)
eigen_data25 = eigen(cov25)
eigen_data27 = eigen(cov27)

# Transform eigenvectors into a data frame
eigenvec5 = as.data.frame(eigen_data5$vectors)
eigenvec6 = as.data.frame(eigen_data6$vectors)
eigenvec10 = as.data.frame(eigen_data10$vectors)
eigenvec11 = as.data.frame(eigen_data11$vectors)
eigenvec12 = as.data.frame(eigen_data12$vectors)
eigenvec13 = as.data.frame(eigen_data13$vectors)
eigenvec16 = as.data.frame(eigen_data16$vectors)
eigenvec18 = as.data.frame(eigen_data18$vectors)
eigenvec20 = as.data.frame(eigen_data20$vectors)
eigenvec22 = as.data.frame(eigen_data22$vectors)
eigenvec25 = as.data.frame(eigen_data25$vectors)
eigenvec27 = as.data.frame(eigen_data27$vectors)

# Add info about the individuals and populations in the data frame
eigenvec5 = cbind(individuals, eigenvec5)
eigenvec6 = cbind(individuals, eigenvec6)
eigenvec10 = cbind(individuals, eigenvec10)
eigenvec11 = cbind(individuals, eigenvec11)
eigenvec12 = cbind(individuals, eigenvec12)
eigenvec13 = cbind(individuals, eigenvec13)
eigenvec16 = cbind(individuals, eigenvec16)
eigenvec18 = cbind(individuals, eigenvec18)
eigenvec20 = cbind(individuals, eigenvec20)
eigenvec22 = cbind(individuals, eigenvec22)
eigenvec25 = cbind(individuals, eigenvec25)
eigenvec27 = cbind(individuals, eigenvec27)

# Change column name
names(eigenvec5)[1] = "Individual"
names(eigenvec6)[1] = "Individual"
names(eigenvec10)[1] = "Individual"
names(eigenvec11)[1] = "Individual"
names(eigenvec12)[1] = "Individual"
names(eigenvec13)[1] = "Individual"
names(eigenvec16)[1] = "Individual"
names(eigenvec18)[1] = "Individual"
names(eigenvec20)[1] = "Individual"
names(eigenvec22)[1] = "Individual"
names(eigenvec25)[1] = "Individual"
names(eigenvec27)[1] = "Individual"

# Reorder by name of individual
eigenvec5 = eigenvec5[order(eigenvec5$Individual),]
eigenvec6 = eigenvec6[order(eigenvec6$Individual),]
eigenvec10 = eigenvec10[order(eigenvec10$Individual),]
eigenvec11 = eigenvec11[order(eigenvec11$Individual),]
eigenvec12 = eigenvec12[order(eigenvec12$Individual),]
eigenvec13 = eigenvec13[order(eigenvec13$Individual),]
eigenvec16 = eigenvec16[order(eigenvec16$Individual),]
eigenvec18 = eigenvec18[order(eigenvec18$Individual),]
eigenvec20 = eigenvec20[order(eigenvec20$Individual),]
eigenvec22 = eigenvec22[order(eigenvec22$Individual),]
eigenvec25 = eigenvec25[order(eigenvec25$Individual),]
eigenvec27 = eigenvec27[order(eigenvec27$Individual),]

# Import information on the individuals
individuals_information = read.table("<File_with_Information_on_Individuals>", header=T)

# Merge the two dataframes
pca5 = merge(individuals_information, eigenvec5, by="Individual")
pca6 = merge(individuals_information, eigenvec6, by="Individual")
pca10 = merge(individuals_information, eigenvec10, by="Individual")
pca11 = merge(individuals_information, eigenvec11, by="Individual")
pca12 = merge(individuals_information, eigenvec12, by="Individual")
pca13 = merge(individuals_information, eigenvec13, by="Individual")
pca16 = merge(individuals_information, eigenvec16, by="Individual")
pca18 = merge(individuals_information, eigenvec18, by="Individual")
pca20 = merge(individuals_information, eigenvec20, by="Individual")
pca22 = merge(individuals_information, eigenvec22, by="Individual")
pca25 = merge(individuals_information, eigenvec25, by="Individual")
pca27 = merge(individuals_information, eigenvec27, by="Individual")


#### Join PCA plots ####

plotpca5 = ggplot(data=pca5, aes(x=V1, y=V2, col=Concrete_Location)) +
  geom_point(aes(shape=Type), size=3, stroke=1, alpha=0.8) +
  ggtitle("Cluster 5") +
  xlab(paste0("PC1", sep=" (", round((eigen_data5$values[1]/sum(eigen_data5$values)*100),2), "%)")) +
  ylab(paste0("PC2", sep=" (", round((eigen_data5$values[2]/sum(eigen_data5$values)*100),2), "%)")) +
  #geom_text(aes(label=pca_notpacific$Individual), hjust=0, vjust=0, size=2) +
  scale_shape_manual(values=c("Freshwater"=16, "Marine"=1)) +
  scale_color_manual(name="Location", values=c("Adriatic_Sea"="orange", "Alaska"="darkolivegreen4", "North_Scandinavia"="royalblue3", "East_Russia"="tomato4", "English_Channel"="chartreuse3", "Finnish_Golf"="purple3", "Iberian_Peninsula"="red", "North_Sea"="turquoise4"),
                     labels=c("Adriatic_Sea"="Adriatic Sea", "Alaska"="Alaska", "North_Scandinavia"="North Scandinavia", "East_Russia"="East Russia", "English_Channel"="English Channel", "Finnish_Golf"="Gulf of Finland", "Iberian_Peninsula"="Iberian Peninsula", "North_Sea"="North Sea")) +
  guides(colour=guide_legend(order=1, title="Geographic Region", override.aes=list(size=3, alpha=1)), shape=guide_legend(order=2, title="Habitat Type",override.aes=list(size=3, alpha=1))) +
  scale_x_reverse() +
  theme_bw() +
  theme(plot.title = element_text(hjust=0.5))

plotpca6 = ggplot(data=pca6, aes(x=V1, y=V2, col=Concrete_Location)) +
  geom_point(aes(shape=Type), size=3, stroke=1, alpha=0.8) +
  ggtitle("Cluster 6") +
  xlab(paste0("PC1", sep=" (", round((eigen_data6$values[1]/sum(eigen_data6$values)*100),2), "%)")) +
  ylab(paste0("PC2", sep=" (", round((eigen_data6$values[2]/sum(eigen_data6$values)*100),2), "%)")) +
  #geom_text(aes(label=pca_notpacific$Individual), hjust=0, vjust=0, size=2) +
  scale_shape_manual(values=c("Freshwater"=16, "Marine"=1)) +
  scale_color_manual(name="Location", values=c("Adriatic_Sea"="orange", "Alaska"="darkolivegreen4", "North_Scandinavia"="royalblue3", "East_Russia"="tomato4", "English_Channel"="chartreuse3", "Finnish_Golf"="purple3", "Iberian_Peninsula"="red", "North_Sea"="turquoise4"),
                     labels=c("Adriatic_Sea"="Adriatic Sea", "Alaska"="Alaska", "North_Scandinavia"="North Scandinavia", "East_Russia"="East Russia", "English_Channel"="English Channel", "Finnish_Golf"="Gulf of Finland", "Iberian_Peninsula"="Iberian Peninsula", "North_Sea"="North Sea")) +
  guides(colour=guide_legend(order=1, title="Geographic Region", override.aes=list(size=3, alpha=1)), shape=guide_legend(order=2, title="Habitat Type",override.aes=list(size=3, alpha=1))) +
  theme_bw() +
  theme(plot.title = element_text(hjust=0.5))

plotpca10 = ggplot(data=pca10, aes(x=V1, y=V2, col=Concrete_Location)) +
  geom_point(aes(shape=Type), size=3, stroke=1, alpha=0.8) +
  ggtitle("Cluster 10") +
  xlab(paste0("PC1", sep=" (", round((eigen_data10$values[1]/sum(eigen_data10$values)*100),2), "%)")) +
  ylab(paste0("PC2", sep=" (", round((eigen_data10$values[2]/sum(eigen_data10$values)*100),2), "%)")) +
  #geom_text(aes(label=pca_notpacific$Individual), hjust=0, vjust=0, size=2) +
  scale_shape_manual(values=c("Freshwater"=16, "Marine"=1)) +
  scale_color_manual(name="Location", values=c("Adriatic_Sea"="orange", "Alaska"="darkolivegreen4", "North_Scandinavia"="royalblue3", "East_Russia"="tomato4", "English_Channel"="chartreuse3", "Finnish_Golf"="purple3", "Iberian_Peninsula"="red", "North_Sea"="turquoise4"),
                     labels=c("Adriatic_Sea"="Adriatic Sea", "Alaska"="Alaska", "North_Scandinavia"="North Scandinavia", "East_Russia"="East Russia", "English_Channel"="English Channel", "Finnish_Golf"="Gulf of Finland", "Iberian_Peninsula"="Iberian Peninsula", "North_Sea"="North Sea")) +
  guides(colour=guide_legend(order=1, title="Geographic Region", override.aes=list(size=3, alpha=1)), shape=guide_legend(order=2, title="Habitat Type",override.aes=list(size=3, alpha=1))) +
  scale_x_reverse() +
  theme_bw() +
  theme(plot.title = element_text(hjust=0.5))

plotpca11 = ggplot(data=pca11, aes(x=V1, y=V2, col=Concrete_Location)) +
  geom_point(aes(shape=Type), size=3, stroke=1, alpha=0.8) +
  ggtitle("Cluster 11") +
  xlab(paste0("PC1", sep=" (", round((eigen_data11$values[1]/sum(eigen_data11$values)*100),2), "%)")) +
  ylab(paste0("PC2", sep=" (", round((eigen_data11$values[2]/sum(eigen_data11$values)*100),2), "%)")) +
  #geom_text(aes(label=pca_notpacific$Individual), hjust=0, vjust=0, size=2) +
  scale_shape_manual(values=c("Freshwater"=16, "Marine"=1)) +
  scale_color_manual(name="Location", values=c("Adriatic_Sea"="orange", "Alaska"="darkolivegreen4", "North_Scandinavia"="royalblue3", "East_Russia"="tomato4", "English_Channel"="chartreuse3", "Finnish_Golf"="purple3", "Iberian_Peninsula"="red", "North_Sea"="turquoise4"),
                     labels=c("Adriatic_Sea"="Adriatic Sea", "Alaska"="Alaska", "North_Scandinavia"="North Scandinavia", "East_Russia"="East Russia", "English_Channel"="English Channel", "Finnish_Golf"="Gulf of Finland", "Iberian_Peninsula"="Iberian Peninsula", "North_Sea"="North Sea")) +
  guides(colour=guide_legend(order=1, title="Geographic Region", override.aes=list(size=3, alpha=1)), shape=guide_legend(order=2, title="Habitat Type",override.aes=list(size=3, alpha=1))) +
  theme_bw() +
  theme(plot.title = element_text(hjust=0.5))

plotpca12 = ggplot(data=pca12, aes(x=V1, y=V2, col=Concrete_Location)) +
  geom_point(aes(shape=Type), size=3, stroke=1, alpha=0.8) +
  ggtitle("Cluster 12") +
  xlab(paste0("PC1", sep=" (", round((eigen_data12$values[1]/sum(eigen_data12$values)*100),2), "%)")) +
  ylab(paste0("PC2", sep=" (", round((eigen_data12$values[2]/sum(eigen_data12$values)*100),2), "%)")) +
  #geom_text(aes(label=pca_notpacific$Individual), hjust=0, vjust=0, size=2) +
  scale_shape_manual(values=c("Freshwater"=16, "Marine"=1)) +
  scale_color_manual(name="Location", values=c("Adriatic_Sea"="orange", "Alaska"="darkolivegreen4", "North_Scandinavia"="royalblue3", "East_Russia"="tomato4", "English_Channel"="chartreuse3", "Finnish_Golf"="purple3", "Iberian_Peninsula"="red", "North_Sea"="turquoise4"),
                     labels=c("Adriatic_Sea"="Adriatic Sea", "Alaska"="Alaska", "North_Scandinavia"="North Scandinavia", "East_Russia"="East Russia", "English_Channel"="English Channel", "Finnish_Golf"="Gulf of Finland", "Iberian_Peninsula"="Iberian Peninsula", "North_Sea"="North Sea")) +
  guides(colour=guide_legend(order=1, title="Geographic Region", override.aes=list(size=3, alpha=1)), shape=guide_legend(order=2, title="Habitat Type",override.aes=list(size=3, alpha=1))) +
  theme_bw() +
  theme(plot.title = element_text(hjust=0.5))

plotpca13 = ggplot(data=pca13, aes(x=V1, y=V2, col=Concrete_Location)) +
  geom_point(aes(shape=Type), size=3, stroke=1, alpha=0.8) +
  ggtitle("Cluster 13") +
  xlab(paste0("PC1", sep=" (", round((eigen_data13$values[1]/sum(eigen_data13$values)*100),2), "%)")) +
  ylab(paste0("PC2", sep=" (", round((eigen_data13$values[2]/sum(eigen_data13$values)*100),2), "%)")) +
  #geom_text(aes(label=pca_notpacific$Individual), hjust=0, vjust=0, size=2) +
  scale_shape_manual(values=c("Freshwater"=16, "Marine"=1)) +
  scale_color_manual(name="Location", values=c("Adriatic_Sea"="orange", "Alaska"="darkolivegreen4", "North_Scandinavia"="royalblue3", "East_Russia"="tomato4", "English_Channel"="chartreuse3", "Finnish_Golf"="purple3", "Iberian_Peninsula"="red", "North_Sea"="turquoise4"),
                     labels=c("Adriatic_Sea"="Adriatic Sea", "Alaska"="Alaska", "North_Scandinavia"="North Scandinavia", "East_Russia"="East Russia", "English_Channel"="English Channel", "Finnish_Golf"="Gulf of Finland", "Iberian_Peninsula"="Iberian Peninsula", "North_Sea"="North Sea")) +
  guides(colour=guide_legend(order=1, title="Geographic Region", override.aes=list(size=3, alpha=1)), shape=guide_legend(order=2, title="Habitat Type",override.aes=list(size=3, alpha=1))) +
  theme_bw() +
  theme(plot.title = element_text(hjust=0.5))

plotpca16 = ggplot(data=pca16, aes(x=V1, y=V2, col=Concrete_Location)) +
  geom_point(aes(shape=Type), size=3, stroke=1, alpha=0.8) +
  ggtitle("Cluster 16") +
  xlab(paste0("PC1", sep=" (", round((eigen_data16$values[1]/sum(eigen_data16$values)*100),2), "%)")) +
  ylab(paste0("PC2", sep=" (", round((eigen_data16$values[2]/sum(eigen_data16$values)*100),2), "%)")) +
  #geom_text(aes(label=pca_notpacific$Individual), hjust=0, vjust=0, size=2) +
  scale_shape_manual(values=c("Freshwater"=16, "Marine"=1)) +
  scale_color_manual(name="Location", values=c("Adriatic_Sea"="orange", "Alaska"="darkolivegreen4", "North_Scandinavia"="royalblue3", "East_Russia"="tomato4", "English_Channel"="chartreuse3", "Finnish_Golf"="purple3", "Iberian_Peninsula"="red", "North_Sea"="turquoise4"),
                     labels=c("Adriatic_Sea"="Adriatic Sea", "Alaska"="Alaska", "North_Scandinavia"="North Scandinavia", "East_Russia"="East Russia", "English_Channel"="English Channel", "Finnish_Golf"="Gulf of Finland", "Iberian_Peninsula"="Iberian Peninsula", "North_Sea"="North Sea")) +
  guides(colour=guide_legend(order=1, title="Geographic Region", override.aes=list(size=3, alpha=1)), shape=guide_legend(order=2, title="Habitat Type",override.aes=list(size=3, alpha=1))) +
  scale_x_reverse() +
  theme_bw() +
  theme(plot.title = element_text(hjust=0.5))

plotpca18 = ggplot(data=pca18, aes(x=V1, y=V2, col=Concrete_Location)) +
  geom_point(aes(shape=Type), size=3, stroke=1, alpha=0.8) +
  ggtitle("Cluster 18") +
  xlab(paste0("PC1", sep=" (", round((eigen_data18$values[1]/sum(eigen_data18$values)*100),2), "%)")) +
  ylab(paste0("PC2", sep=" (", round((eigen_data18$values[2]/sum(eigen_data18$values)*100),2), "%)")) +
  #geom_text(aes(label=pca_notpacific$Individual), hjust=0, vjust=0, size=2) +
  scale_shape_manual(values=c("Freshwater"=16, "Marine"=1)) +
  scale_color_manual(name="Location", values=c("Adriatic_Sea"="orange", "Alaska"="darkolivegreen4", "North_Scandinavia"="royalblue3", "East_Russia"="tomato4", "English_Channel"="chartreuse3", "Finnish_Golf"="purple3", "Iberian_Peninsula"="red", "North_Sea"="turquoise4"),
                     labels=c("Adriatic_Sea"="Adriatic Sea", "Alaska"="Alaska", "North_Scandinavia"="North Scandinavia", "East_Russia"="East Russia", "English_Channel"="English Channel", "Finnish_Golf"="Gulf of Finland", "Iberian_Peninsula"="Iberian Peninsula", "North_Sea"="North Sea")) +
  guides(colour=guide_legend(order=1, title="Geographic Region", override.aes=list(size=3, alpha=1)), shape=guide_legend(order=2, title="Habitat Type",override.aes=list(size=3, alpha=1))) +
  scale_x_reverse() +
  theme_bw() +
  theme(plot.title = element_text(hjust=0.5))

plotpca20 = ggplot(data=pca20, aes(x=V1, y=V2, col=Concrete_Location)) +
  geom_point(aes(shape=Type), size=3, stroke=1, alpha=0.8) +
  ggtitle("Cluster 20") +
  xlab(paste0("PC1", sep=" (", round((eigen_data20$values[1]/sum(eigen_data20$values)*100),2), "%)")) +
  ylab(paste0("PC2", sep=" (", round((eigen_data20$values[2]/sum(eigen_data20$values)*100),2), "%)")) +
  #geom_text(aes(label=pca_notpacific$Individual), hjust=0, vjust=0, size=2) +
  scale_shape_manual(values=c("Freshwater"=16, "Marine"=1)) +
  scale_color_manual(name="Location", values=c("Adriatic_Sea"="orange", "Alaska"="darkolivegreen4", "North_Scandinavia"="royalblue3", "East_Russia"="tomato4", "English_Channel"="chartreuse3", "Finnish_Golf"="purple3", "Iberian_Peninsula"="red", "North_Sea"="turquoise4"),
                     labels=c("Adriatic_Sea"="Adriatic Sea", "Alaska"="Alaska", "North_Scandinavia"="North Scandinavia", "East_Russia"="East Russia", "English_Channel"="English Channel", "Finnish_Golf"="Gulf of Finland", "Iberian_Peninsula"="Iberian Peninsula", "North_Sea"="North Sea")) +
  guides(colour=guide_legend(order=1, title="Geographic Region", override.aes=list(size=3, alpha=1)), shape=guide_legend(order=2, title="Habitat Type",override.aes=list(size=3, alpha=1))) +
  theme_bw() +
  theme(plot.title = element_text(hjust=0.5))

plotpca22 = ggplot(data=pca22, aes(x=V1, y=V2, col=Concrete_Location)) +
  geom_point(aes(shape=Type), size=3, stroke=1, alpha=0.8) +
  ggtitle("Cluster 22") +
  xlab(paste0("PC1", sep=" (", round((eigen_data22$values[1]/sum(eigen_data22$values)*100),2), "%)")) +
  ylab(paste0("PC2", sep=" (", round((eigen_data22$values[2]/sum(eigen_data22$values)*100),2), "%)")) +
  #geom_text(aes(label=pca_notpacific$Individual), hjust=0, vjust=0, size=2) +
  scale_shape_manual(values=c("Freshwater"=16, "Marine"=1)) +
  scale_color_manual(name="Location", values=c("Adriatic_Sea"="orange", "Alaska"="darkolivegreen4", "North_Scandinavia"="royalblue3", "East_Russia"="tomato4", "English_Channel"="chartreuse3", "Finnish_Golf"="purple3", "Iberian_Peninsula"="red", "North_Sea"="turquoise4"),
                     labels=c("Adriatic_Sea"="Adriatic Sea", "Alaska"="Alaska", "North_Scandinavia"="North Scandinavia", "East_Russia"="East Russia", "English_Channel"="English Channel", "Finnish_Golf"="Gulf of Finland", "Iberian_Peninsula"="Iberian Peninsula", "North_Sea"="North Sea")) +
  guides(colour=guide_legend(order=1, title="Geographic Region", override.aes=list(size=3, alpha=1)), shape=guide_legend(order=2, title="Habitat Type",override.aes=list(size=3, alpha=1))) +
  theme_bw() +
  theme(plot.title = element_text(hjust=0.5))

plotpca25 = ggplot(data=pca25, aes(x=V1, y=V2, col=Concrete_Location)) +
  geom_point(aes(shape=Type), size=3, stroke=1, alpha=0.8) +
  ggtitle("Cluster 25") +
  xlab(paste0("PC1", sep=" (", round((eigen_data25$values[1]/sum(eigen_data25$values)*100),2), "%)")) +
  ylab(paste0("PC2", sep=" (", round((eigen_data25$values[2]/sum(eigen_data25$values)*100),2), "%)")) +
  #geom_text(aes(label=pca_notpacific$Individual), hjust=0, vjust=0, size=2) +
  scale_shape_manual(values=c("Freshwater"=16, "Marine"=1)) +
  scale_color_manual(name="Location", values=c("Adriatic_Sea"="orange", "Alaska"="darkolivegreen4", "North_Scandinavia"="royalblue3", "East_Russia"="tomato4", "English_Channel"="chartreuse3", "Finnish_Golf"="purple3", "Iberian_Peninsula"="red", "North_Sea"="turquoise4"),
                     labels=c("Adriatic_Sea"="Adriatic Sea", "Alaska"="Alaska", "North_Scandinavia"="North Scandinavia", "East_Russia"="East Russia", "English_Channel"="English Channel", "Finnish_Golf"="Gulf of Finland", "Iberian_Peninsula"="Iberian Peninsula", "North_Sea"="North Sea")) +
  guides(colour=guide_legend(order=1, title="Geographic Region", override.aes=list(size=3, alpha=1)), shape=guide_legend(order=2, title="Habitat Type",override.aes=list(size=3, alpha=1))) +
  scale_x_reverse() +
  theme_bw() +
  theme(plot.title = element_text(hjust=0.5))

plotpca27 = ggplot(data=pca27, aes(x=V1, y=V2, col=Concrete_Location)) +
  geom_point(aes(shape=Type), size=3, stroke=1, alpha=0.8) +
  ggtitle("Cluster 27") +
  xlab(paste0("PC1", sep=" (", round((eigen_data27$values[1]/sum(eigen_data27$values)*100),2), "%)")) +
  ylab(paste0("PC2", sep=" (", round((eigen_data27$values[2]/sum(eigen_data27$values)*100),2), "%)")) +
  #geom_text(aes(label=pca_notpacific$Individual), hjust=0, vjust=0, size=2) +
  scale_shape_manual(values=c("Freshwater"=16, "Marine"=1)) +
  scale_color_manual(name="Location", values=c("Adriatic_Sea"="orange", "Alaska"="darkolivegreen4", "North_Scandinavia"="royalblue3", "East_Russia"="tomato4", "English_Channel"="chartreuse3", "Finnish_Golf"="purple3", "Iberian_Peninsula"="red", "North_Sea"="turquoise4"),
                     labels=c("Adriatic_Sea"="Adriatic Sea", "Alaska"="Alaska", "North_Scandinavia"="North Scandinavia", "East_Russia"="East Russia", "English_Channel"="English Channel", "Finnish_Golf"="Gulf of Finland", "Iberian_Peninsula"="Iberian Peninsula", "North_Sea"="North Sea")) +
  guides(colour=guide_legend(order=1, title="Geographic Region", override.aes=list(size=3, alpha=1)), shape=guide_legend(order=2, title="Habitat Type",override.aes=list(size=3, alpha=1))) +
  theme_bw() +
  theme(plot.title = element_text(hjust=0.5))

png("PCAs.png", units="in", width=15, height=5, res=900)
  ggarrange(plotpca5,
            plotpca6,
            plotpca10,
            plotpca11,
            plotpca12,
            plotpca13,
            plotpca16,
            plotpca18,
            plotpca20,
            plotpca22,
            plotpca25,
            plotpca27,
            ncol=6, nrow=2, labels="auto", hjust=-1, common.legend = TRUE, legend="right")
dev.off()


#### Join BoxPlots by Population ####

# Copy PCA variables in another variable
bypop5 = pca5
bypop6 = pca6
bypop10 = pca10
bypop11 = pca11
bypop12 = pca12
bypop13 = pca13
bypop16 = pca16
bypop18 = pca18
bypop20 = pca20
bypop22 = pca22
bypop25 = pca25
bypop27 = pca27

# Reorder the factors
bypop5$Population = factor(bypop5$Population, levels=c("RABBIT-SLOUGH", "RUS-ASH-GA", "RUS-AN-GA", "RUS-KHA-GA", "RUS-LEV-GA", "NOR-BAR-GA", "NOR-SBJ-GA", "GAC-RUS-PRI", "GAC-NOR-KRI", "GAC-SWE-FIS", "NOR-MYR-GA", "BOOT-LAKE", "BEAR-PAW-LAKE", "NOR-KVA-GA", "NOR-SKF-GA", "FIN-KEV-GA", "RUS-SLI-GA", "NOR-ORR-GA", "ENG-BUT-GA", "GA-POR-VO", "GA-POR-RA", "GA-POR-LI", "GA-POR-TE", "GA-POR-MI", "GA-POR-SA", "MU", "NN", "ITA-STE", "BOS-NER-GA", "MON-SKA-GA", "M", "NB", "KR")) 
bypop6$Population = factor(bypop6$Population, levels=c("RABBIT-SLOUGH", "RUS-ASH-GA", "RUS-AN-GA", "RUS-KHA-GA", "RUS-LEV-GA", "NOR-BAR-GA", "NOR-SBJ-GA", "GAC-RUS-PRI", "GAC-NOR-KRI", "GAC-SWE-FIS", "NOR-MYR-GA", "BOOT-LAKE", "BEAR-PAW-LAKE", "NOR-KVA-GA", "NOR-SKF-GA", "FIN-KEV-GA", "RUS-SLI-GA", "NOR-ORR-GA", "ENG-BUT-GA", "GA-POR-VO", "GA-POR-RA", "GA-POR-LI", "GA-POR-TE", "GA-POR-MI", "GA-POR-SA", "MU", "NN", "ITA-STE", "BOS-NER-GA", "MON-SKA-GA", "M", "NB", "KR")) 
bypop10$Population = factor(bypop10$Population, levels=c("RABBIT-SLOUGH", "RUS-ASH-GA", "RUS-AN-GA", "RUS-KHA-GA", "RUS-LEV-GA", "NOR-BAR-GA", "NOR-SBJ-GA", "GAC-RUS-PRI", "GAC-NOR-KRI", "GAC-SWE-FIS", "NOR-MYR-GA", "BOOT-LAKE", "BEAR-PAW-LAKE", "NOR-KVA-GA", "NOR-SKF-GA", "FIN-KEV-GA", "RUS-SLI-GA", "NOR-ORR-GA", "ENG-BUT-GA", "GA-POR-VO", "GA-POR-RA", "GA-POR-LI", "GA-POR-TE", "GA-POR-MI", "GA-POR-SA", "MU", "NN", "ITA-STE", "BOS-NER-GA", "MON-SKA-GA", "M", "NB", "KR")) 
bypop11$Population = factor(bypop11$Population, levels=c("RABBIT-SLOUGH", "RUS-ASH-GA", "RUS-AN-GA", "RUS-KHA-GA", "RUS-LEV-GA", "NOR-BAR-GA", "NOR-SBJ-GA", "GAC-RUS-PRI", "GAC-NOR-KRI", "GAC-SWE-FIS", "NOR-MYR-GA", "BOOT-LAKE", "BEAR-PAW-LAKE", "NOR-KVA-GA", "NOR-SKF-GA", "FIN-KEV-GA", "RUS-SLI-GA", "NOR-ORR-GA", "ENG-BUT-GA", "GA-POR-VO", "GA-POR-RA", "GA-POR-LI", "GA-POR-TE", "GA-POR-MI", "GA-POR-SA", "MU", "NN", "ITA-STE", "BOS-NER-GA", "MON-SKA-GA", "M", "NB", "KR")) 
bypop12$Population = factor(bypop12$Population, levels=c("RABBIT-SLOUGH", "RUS-ASH-GA", "RUS-AN-GA", "RUS-KHA-GA", "RUS-LEV-GA", "NOR-BAR-GA", "NOR-SBJ-GA", "GAC-RUS-PRI", "GAC-NOR-KRI", "GAC-SWE-FIS", "NOR-MYR-GA", "BOOT-LAKE", "BEAR-PAW-LAKE", "NOR-KVA-GA", "NOR-SKF-GA", "FIN-KEV-GA", "RUS-SLI-GA", "NOR-ORR-GA", "ENG-BUT-GA", "GA-POR-VO", "GA-POR-RA", "GA-POR-LI", "GA-POR-TE", "GA-POR-MI", "GA-POR-SA", "MU", "NN", "ITA-STE", "BOS-NER-GA", "MON-SKA-GA", "M", "NB", "KR")) 
bypop13$Population = factor(bypop13$Population, levels=c("RABBIT-SLOUGH", "RUS-ASH-GA", "RUS-AN-GA", "RUS-KHA-GA", "RUS-LEV-GA", "NOR-BAR-GA", "NOR-SBJ-GA", "GAC-RUS-PRI", "GAC-NOR-KRI", "GAC-SWE-FIS", "NOR-MYR-GA", "BOOT-LAKE", "BEAR-PAW-LAKE", "NOR-KVA-GA", "NOR-SKF-GA", "FIN-KEV-GA", "RUS-SLI-GA", "NOR-ORR-GA", "ENG-BUT-GA", "GA-POR-VO", "GA-POR-RA", "GA-POR-LI", "GA-POR-TE", "GA-POR-MI", "GA-POR-SA", "MU", "NN", "ITA-STE", "BOS-NER-GA", "MON-SKA-GA", "M", "NB", "KR")) 
bypop16$Population = factor(bypop16$Population, levels=c("RABBIT-SLOUGH", "RUS-ASH-GA", "RUS-AN-GA", "RUS-KHA-GA", "RUS-LEV-GA", "NOR-BAR-GA", "NOR-SBJ-GA", "GAC-RUS-PRI", "GAC-NOR-KRI", "GAC-SWE-FIS", "NOR-MYR-GA", "BOOT-LAKE", "BEAR-PAW-LAKE", "NOR-KVA-GA", "NOR-SKF-GA", "FIN-KEV-GA", "RUS-SLI-GA", "NOR-ORR-GA", "ENG-BUT-GA", "GA-POR-VO", "GA-POR-RA", "GA-POR-LI", "GA-POR-TE", "GA-POR-MI", "GA-POR-SA", "MU", "NN", "ITA-STE", "BOS-NER-GA", "MON-SKA-GA", "M", "NB", "KR")) 
bypop18$Population = factor(bypop18$Population, levels=c("RABBIT-SLOUGH", "RUS-ASH-GA", "RUS-AN-GA", "RUS-KHA-GA", "RUS-LEV-GA", "NOR-BAR-GA", "NOR-SBJ-GA", "GAC-RUS-PRI", "GAC-NOR-KRI", "GAC-SWE-FIS", "NOR-MYR-GA", "BOOT-LAKE", "BEAR-PAW-LAKE", "NOR-KVA-GA", "NOR-SKF-GA", "FIN-KEV-GA", "RUS-SLI-GA", "NOR-ORR-GA", "ENG-BUT-GA", "GA-POR-VO", "GA-POR-RA", "GA-POR-LI", "GA-POR-TE", "GA-POR-MI", "GA-POR-SA", "MU", "NN", "ITA-STE", "BOS-NER-GA", "MON-SKA-GA", "M", "NB", "KR")) 
bypop20$Population = factor(bypop20$Population, levels=c("RABBIT-SLOUGH", "RUS-ASH-GA", "RUS-AN-GA", "RUS-KHA-GA", "RUS-LEV-GA", "NOR-BAR-GA", "NOR-SBJ-GA", "GAC-RUS-PRI", "GAC-NOR-KRI", "GAC-SWE-FIS", "NOR-MYR-GA", "BOOT-LAKE", "BEAR-PAW-LAKE", "NOR-KVA-GA", "NOR-SKF-GA", "FIN-KEV-GA", "RUS-SLI-GA", "NOR-ORR-GA", "ENG-BUT-GA", "GA-POR-VO", "GA-POR-RA", "GA-POR-LI", "GA-POR-TE", "GA-POR-MI", "GA-POR-SA", "MU", "NN", "ITA-STE", "BOS-NER-GA", "MON-SKA-GA", "M", "NB", "KR")) 
bypop22$Population = factor(bypop22$Population, levels=c("RABBIT-SLOUGH", "RUS-ASH-GA", "RUS-AN-GA", "RUS-KHA-GA", "RUS-LEV-GA", "NOR-BAR-GA", "NOR-SBJ-GA", "GAC-RUS-PRI", "GAC-NOR-KRI", "GAC-SWE-FIS", "NOR-MYR-GA", "BOOT-LAKE", "BEAR-PAW-LAKE", "NOR-KVA-GA", "NOR-SKF-GA", "FIN-KEV-GA", "RUS-SLI-GA", "NOR-ORR-GA", "ENG-BUT-GA", "GA-POR-VO", "GA-POR-RA", "GA-POR-LI", "GA-POR-TE", "GA-POR-MI", "GA-POR-SA", "MU", "NN", "ITA-STE", "BOS-NER-GA", "MON-SKA-GA", "M", "NB", "KR")) 
bypop25$Population = factor(bypop25$Population, levels=c("RABBIT-SLOUGH", "RUS-ASH-GA", "RUS-AN-GA", "RUS-KHA-GA", "RUS-LEV-GA", "NOR-BAR-GA", "NOR-SBJ-GA", "GAC-RUS-PRI", "GAC-NOR-KRI", "GAC-SWE-FIS", "NOR-MYR-GA", "BOOT-LAKE", "BEAR-PAW-LAKE", "NOR-KVA-GA", "NOR-SKF-GA", "FIN-KEV-GA", "RUS-SLI-GA", "NOR-ORR-GA", "ENG-BUT-GA", "GA-POR-VO", "GA-POR-RA", "GA-POR-LI", "GA-POR-TE", "GA-POR-MI", "GA-POR-SA", "MU", "NN", "ITA-STE", "BOS-NER-GA", "MON-SKA-GA", "M", "NB", "KR")) 
bypop27$Population = factor(bypop27$Population, levels=c("RABBIT-SLOUGH", "RUS-ASH-GA", "RUS-AN-GA", "RUS-KHA-GA", "RUS-LEV-GA", "NOR-BAR-GA", "NOR-SBJ-GA", "GAC-RUS-PRI", "GAC-NOR-KRI", "GAC-SWE-FIS", "NOR-MYR-GA", "BOOT-LAKE", "BEAR-PAW-LAKE", "NOR-KVA-GA", "NOR-SKF-GA", "FIN-KEV-GA", "RUS-SLI-GA", "NOR-ORR-GA", "ENG-BUT-GA", "GA-POR-VO", "GA-POR-RA", "GA-POR-LI", "GA-POR-TE", "GA-POR-MI", "GA-POR-SA", "MU", "NN", "ITA-STE", "BOS-NER-GA", "MON-SKA-GA", "M", "NB", "KR")) 


boxplotpop5 = ggplot(data=bypop5, aes(x=Population, y=V1, fill=Concrete_Location)) +
  geom_boxplot() +
  #geom_boxplot(width=0.1, color="grey") +
  ggtitle("Cluster 5") +
  xlab("Population") +
  ylab(paste0("PC1", sep=" (", round((eigen_data5$values[1]/sum(eigen_data5$values)*100),2), "%)")) +
  scale_x_discrete(limits = rev(levels(bypop5$Population))) +
  scale_y_reverse() +
  scale_fill_manual(name="Geographic Region", values=c("Adriatic_Sea"="orange", "Alaska"="darkolivegreen4", "North_Scandinavia"="royalblue3", "East_Russia"="tomato4", "English_Channel"="chartreuse3", "Finnish_Golf"="purple3", "Iberian_Peninsula"="red", "North_Sea"="turquoise4"),
                    labels=c("Adriatic_Sea"="Adriatic Sea", "Alaska"="Alaska", "North_Scandinavia"="North Scandinavia", "East_Russia"="East Russia", "English_Channel"="English Channel", "Finnish_Golf"="Gulf of Finland", "Iberian_Peninsula"="Iberian Peninsula", "North_Sea"="North Sea")) +
  coord_flip() +
  theme_bw() +
  theme(plot.title = element_text(hjust=0.5),
        axis.text.y=element_text(color=c("skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "navy", "navy", "navy", "navy", "navy", "navy", "navy", "navy", "navy", "navy", "navy")))

boxplotpop6 = ggplot(data=bypop6, aes(x=Population, y=V1, fill=Concrete_Location)) +
  geom_boxplot() +
  #geom_boxplot(width=0.1, color="grey") +
  ggtitle("Cluster 6") +
  xlab("Population") +
  ylab(paste0("PC1", sep=" (", round((eigen_data6$values[1]/sum(eigen_data6$values)*100),2), "%)")) +
  scale_x_discrete(limits = rev(levels(bypop6$Population))) +
  scale_fill_manual(name="Geographic Region", values=c("Adriatic_Sea"="orange", "Alaska"="darkolivegreen4", "North_Scandinavia"="royalblue3", "East_Russia"="tomato4", "English_Channel"="chartreuse3", "Finnish_Golf"="purple3", "Iberian_Peninsula"="red", "North_Sea"="turquoise4"),
                    labels=c("Adriatic_Sea"="Adriatic Sea", "Alaska"="Alaska", "North_Scandinavia"="North Scandinavia", "East_Russia"="East Russia", "English_Channel"="English Channel", "Finnish_Golf"="Gulf of Finland", "Iberian_Peninsula"="Iberian Peninsula", "North_Sea"="North Sea")) +
  coord_flip() +
  theme_bw() +
  theme(plot.title = element_text(hjust=0.5),
        axis.text.y=element_text(color=c("skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "navy", "navy", "navy", "navy", "navy", "navy", "navy", "navy", "navy", "navy", "navy")))

boxplotpop10 = ggplot(data=bypop10, aes(x=Population, y=V1, fill=Concrete_Location)) +
  geom_boxplot() +
  #geom_boxplot(width=0.1, color="grey") +
  ggtitle("Cluster 10") +
  xlab("Population") +
  ylab(paste0("PC1", sep=" (", round((eigen_data10$values[1]/sum(eigen_data10$values)*100),2), "%)")) +
  scale_x_discrete(limits = rev(levels(bypop10$Population))) +
  scale_y_reverse() +
  scale_fill_manual(name="Geographic Region", values=c("Adriatic_Sea"="orange", "Alaska"="darkolivegreen4", "North_Scandinavia"="royalblue3", "East_Russia"="tomato4", "English_Channel"="chartreuse3", "Finnish_Golf"="purple3", "Iberian_Peninsula"="red", "North_Sea"="turquoise4"),
                    labels=c("Adriatic_Sea"="Adriatic Sea", "Alaska"="Alaska", "North_Scandinavia"="North Scandinavia", "East_Russia"="East Russia", "English_Channel"="English Channel", "Finnish_Golf"="Gulf of Finland", "Iberian_Peninsula"="Iberian Peninsula", "North_Sea"="North Sea")) +
  coord_flip() +
  theme_bw() +
  theme(plot.title = element_text(hjust=0.5),
        axis.text.y=element_text(color=c("skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "navy", "navy", "navy", "navy", "navy", "navy", "navy", "navy", "navy", "navy", "navy")))

boxplotpop11 = ggplot(data=bypop11, aes(x=Population, y=V1, fill=Concrete_Location)) +
  geom_boxplot() +
  #geom_boxplot(width=0.1, color="grey") +
  ggtitle("Cluster 11") +
  xlab("Population") +
  ylab(paste0("PC1", sep=" (", round((eigen_data11$values[1]/sum(eigen_data11$values)*100),2), "%)")) +
  scale_x_discrete(limits = rev(levels(bypop11$Population))) +
  scale_fill_manual(name="Geographic Region", values=c("Adriatic_Sea"="orange", "Alaska"="darkolivegreen4", "North_Scandinavia"="royalblue3", "East_Russia"="tomato4", "English_Channel"="chartreuse3", "Finnish_Golf"="purple3", "Iberian_Peninsula"="red", "North_Sea"="turquoise4"),
                    labels=c("Adriatic_Sea"="Adriatic Sea", "Alaska"="Alaska", "North_Scandinavia"="North Scandinavia", "East_Russia"="East Russia", "English_Channel"="English Channel", "Finnish_Golf"="Gulf of Finland", "Iberian_Peninsula"="Iberian Peninsula", "North_Sea"="North Sea")) +
  coord_flip() +
  theme_bw() +
  theme(plot.title = element_text(hjust=0.5),
        axis.text.y=element_text(color=c("skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "navy", "navy", "navy", "navy", "navy", "navy", "navy", "navy", "navy", "navy", "navy")))

boxplotpop12 = ggplot(data=bypop12, aes(x=Population, y=V1, fill=Concrete_Location)) +
  geom_boxplot() +
  #geom_boxplot(width=0.1, color="grey") +
  ggtitle("Cluster 12") +
  xlab("Population") +
  ylab(paste0("PC1", sep=" (", round((eigen_data12$values[1]/sum(eigen_data12$values)*100),2), "%)")) +
  scale_x_discrete(limits = rev(levels(bypop12$Population))) +
  scale_fill_manual(name="Geographic Region", values=c("Adriatic_Sea"="orange", "Alaska"="darkolivegreen4", "North_Scandinavia"="royalblue3", "East_Russia"="tomato4", "English_Channel"="chartreuse3", "Finnish_Golf"="purple3", "Iberian_Peninsula"="red", "North_Sea"="turquoise4"),
                    labels=c("Adriatic_Sea"="Adriatic Sea", "Alaska"="Alaska", "North_Scandinavia"="North Scandinavia", "East_Russia"="East Russia", "English_Channel"="English Channel", "Finnish_Golf"="Gulf of Finland", "Iberian_Peninsula"="Iberian Peninsula", "North_Sea"="North Sea")) +
  coord_flip() +
  theme_bw() +
  theme(plot.title = element_text(hjust=0.5),
        axis.text.y=element_text(color=c("skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "navy", "navy", "navy", "navy", "navy", "navy", "navy", "navy", "navy", "navy", "navy")))

boxplotpop13 = ggplot(data=bypop13, aes(x=Population, y=V1, fill=Concrete_Location)) +
  geom_boxplot() +
  #geom_boxplot(width=0.1, color="grey") +
  ggtitle("Cluster 13") +
  xlab("Population") +
  ylab(paste0("PC1", sep=" (", round((eigen_data13$values[1]/sum(eigen_data13$values)*100),2), "%)")) +
  scale_x_discrete(limits = rev(levels(bypop13$Population))) +
  scale_fill_manual(name="Geographic Region", values=c("Adriatic_Sea"="orange", "Alaska"="darkolivegreen4", "North_Scandinavia"="royalblue3", "East_Russia"="tomato4", "English_Channel"="chartreuse3", "Finnish_Golf"="purple3", "Iberian_Peninsula"="red", "North_Sea"="turquoise4"),
                    labels=c("Adriatic_Sea"="Adriatic Sea", "Alaska"="Alaska", "North_Scandinavia"="North Scandinavia", "East_Russia"="East Russia", "English_Channel"="English Channel", "Finnish_Golf"="Gulf of Finland", "Iberian_Peninsula"="Iberian Peninsula", "North_Sea"="North Sea")) +
  coord_flip() +
  theme_bw() +
  theme(plot.title = element_text(hjust=0.5),
        axis.text.y=element_text(color=c("skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "navy", "navy", "navy", "navy", "navy", "navy", "navy", "navy", "navy", "navy", "navy")))

boxplotpop16 = ggplot(data=bypop16, aes(x=Population, y=V1, fill=Concrete_Location)) +
  geom_boxplot() +
  #geom_boxplot(width=0.1, color="grey") +
  ggtitle("Cluster 16") +
  xlab("Population") +
  ylab(paste0("PC1", sep=" (", round((eigen_data16$values[1]/sum(eigen_data16$values)*100),2), "%)")) +
  scale_x_discrete(limits = rev(levels(bypop16$Population))) +
  scale_y_reverse() +
  scale_fill_manual(name="Geographic Region", values=c("Adriatic_Sea"="orange", "Alaska"="darkolivegreen4", "North_Scandinavia"="royalblue3", "East_Russia"="tomato4", "English_Channel"="chartreuse3", "Finnish_Golf"="purple3", "Iberian_Peninsula"="red", "North_Sea"="turquoise4"),
                    labels=c("Adriatic_Sea"="Adriatic Sea", "Alaska"="Alaska", "North_Scandinavia"="North Scandinavia", "East_Russia"="East Russia", "English_Channel"="English Channel", "Finnish_Golf"="Gulf of Finland", "Iberian_Peninsula"="Iberian Peninsula", "North_Sea"="North Sea")) +
  coord_flip() +
  theme_bw() +
  theme(plot.title = element_text(hjust=0.5),
        axis.text.y=element_text(color=c("skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "navy", "navy", "navy", "navy", "navy", "navy", "navy", "navy", "navy", "navy", "navy")))

boxplotpop18 = ggplot(data=bypop18, aes(x=Population, y=V1, fill=Concrete_Location)) +
  geom_boxplot() +
  #geom_boxplot(width=0.1, color="grey") +
  ggtitle("Cluster 18") +
  xlab("Population") +
  ylab(paste0("PC1", sep=" (", round((eigen_data18$values[1]/sum(eigen_data18$values)*100),2), "%)")) +
  scale_x_discrete(limits = rev(levels(bypop18$Population))) +
  scale_y_reverse() +
  scale_fill_manual(name="Geographic Region", values=c("Adriatic_Sea"="orange", "Alaska"="darkolivegreen4", "North_Scandinavia"="royalblue3", "East_Russia"="tomato4", "English_Channel"="chartreuse3", "Finnish_Golf"="purple3", "Iberian_Peninsula"="red", "North_Sea"="turquoise4"),
                    labels=c("Adriatic_Sea"="Adriatic Sea", "Alaska"="Alaska", "North_Scandinavia"="North Scandinavia", "East_Russia"="East Russia", "English_Channel"="English Channel", "Finnish_Golf"="Gulf of Finland", "Iberian_Peninsula"="Iberian Peninsula", "North_Sea"="North Sea")) +
  coord_flip() +
  theme_bw() +
  theme(plot.title = element_text(hjust=0.5),
        axis.text.y=element_text(color=c("skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "navy", "navy", "navy", "navy", "navy", "navy", "navy", "navy", "navy", "navy", "navy")))

boxplotpop20 = ggplot(data=bypop20, aes(x=Population, y=V1, fill=Concrete_Location)) +
  geom_boxplot() +
  #geom_boxplot(width=0.1, color="grey") +
  ggtitle("Cluster 20") +
  xlab("Population") +
  ylab(paste0("PC1", sep=" (", round((eigen_data20$values[1]/sum(eigen_data20$values)*100),2), "%)")) +
  scale_x_discrete(limits = rev(levels(bypop20$Population))) +
  scale_fill_manual(name="Geographic Region", values=c("Adriatic_Sea"="orange", "Alaska"="darkolivegreen4", "North_Scandinavia"="royalblue3", "East_Russia"="tomato4", "English_Channel"="chartreuse3", "Finnish_Golf"="purple3", "Iberian_Peninsula"="red", "North_Sea"="turquoise4"),
                    labels=c("Adriatic_Sea"="Adriatic Sea", "Alaska"="Alaska", "North_Scandinavia"="North Scandinavia", "East_Russia"="East Russia", "English_Channel"="English Channel", "Finnish_Golf"="Gulf of Finland", "Iberian_Peninsula"="Iberian Peninsula", "North_Sea"="North Sea")) +
  coord_flip() +
  theme_bw() +
  theme(plot.title = element_text(hjust=0.5),
        axis.text.y=element_text(color=c("skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "navy", "navy", "navy", "navy", "navy", "navy", "navy", "navy", "navy", "navy", "navy")))

boxplotpop22 = ggplot(data=bypop22, aes(x=Population, y=V1, fill=Concrete_Location)) +
  geom_boxplot() +
  #geom_boxplot(width=0.1, color="grey") +
  ggtitle("Cluster 22") +
  xlab("Population") +
  ylab(paste0("PC1", sep=" (", round((eigen_data22$values[1]/sum(eigen_data22$values)*100),2), "%)")) +
  scale_x_discrete(limits = rev(levels(bypop22$Population))) +
  scale_fill_manual(name="Geographic Region", values=c("Adriatic_Sea"="orange", "Alaska"="darkolivegreen4", "North_Scandinavia"="royalblue3", "East_Russia"="tomato4", "English_Channel"="chartreuse3", "Finnish_Golf"="purple3", "Iberian_Peninsula"="red", "North_Sea"="turquoise4"),
                    labels=c("Adriatic_Sea"="Adriatic Sea", "Alaska"="Alaska", "North_Scandinavia"="North Scandinavia", "East_Russia"="East Russia", "English_Channel"="English Channel", "Finnish_Golf"="Gulf of Finland", "Iberian_Peninsula"="Iberian Peninsula", "North_Sea"="North Sea")) +
  coord_flip() +
  theme_bw() +
  theme(plot.title = element_text(hjust=0.5),
        axis.text.y=element_text(color=c("skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "navy", "navy", "navy", "navy", "navy", "navy", "navy", "navy", "navy", "navy", "navy")))

boxplotpop25 = ggplot(data=bypop25, aes(x=Population, y=V1, fill=Concrete_Location)) +
  geom_boxplot() +
  #geom_boxplot(width=0.1, color="grey") +
  ggtitle("Cluster 25") +
  xlab("Population") +
  ylab(paste0("PC1", sep=" (", round((eigen_data25$values[1]/sum(eigen_data25$values)*100),2), "%)")) +
  scale_x_discrete(limits = rev(levels(bypop25$Population))) +
  scale_y_reverse() +
  scale_fill_manual(name="Geographic Region", values=c("Adriatic_Sea"="orange", "Alaska"="darkolivegreen4", "North_Scandinavia"="royalblue3", "East_Russia"="tomato4", "English_Channel"="chartreuse3", "Finnish_Golf"="purple3", "Iberian_Peninsula"="red", "North_Sea"="turquoise4"),
                    labels=c("Adriatic_Sea"="Adriatic Sea", "Alaska"="Alaska", "North_Scandinavia"="North Scandinavia", "East_Russia"="East Russia", "English_Channel"="English Channel", "Finnish_Golf"="Gulf of Finland", "Iberian_Peninsula"="Iberian Peninsula", "North_Sea"="North Sea")) +
  coord_flip() +
  theme_bw() +
  theme(plot.title = element_text(hjust=0.5),
        axis.text.y=element_text(color=c("skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "navy", "navy", "navy", "navy", "navy", "navy", "navy", "navy", "navy", "navy", "navy")))

boxplotpop27 = ggplot(data=bypop27, aes(x=Population, y=V1, fill=Concrete_Location)) +
  geom_boxplot() +
  #geom_boxplot(width=0.1, color="grey") +
  ggtitle("Cluster 27") +
  xlab("Population") +
  ylab(paste0("PC1", sep=" (", round((eigen_data27$values[1]/sum(eigen_data27$values)*100),2), "%)")) +
  scale_x_discrete(limits = rev(levels(bypop27$Population))) +
  scale_fill_manual(name="Geographic Region", values=c("Adriatic_Sea"="orange", "Alaska"="darkolivegreen4", "North_Scandinavia"="royalblue3", "East_Russia"="tomato4", "English_Channel"="chartreuse3", "Finnish_Golf"="purple3", "Iberian_Peninsula"="red", "North_Sea"="turquoise4"),
                    labels=c("Adriatic_Sea"="Adriatic Sea", "Alaska"="Alaska", "North_Scandinavia"="North Scandinavia", "East_Russia"="East Russia", "English_Channel"="English Channel", "Finnish_Golf"="Gulf of Finland", "Iberian_Peninsula"="Iberian Peninsula", "North_Sea"="North Sea")) +
  coord_flip() +
  theme_bw() +
  theme(plot.title = element_text(hjust=0.5),
        axis.text.y=element_text(color=c("skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "navy", "navy", "navy", "navy", "navy", "navy", "navy", "navy", "navy", "navy", "navy")))

png("BoxPlots Population Multiple Labels 4x3.png", units="in", width=15, height=20, res=900)
ggarrange(boxplotpop5,
          boxplotpop6 + rremove("ylab"),
          boxplotpop10 + rremove("ylab"),
          boxplotpop11,
          boxplotpop12 + rremove("ylab"),
          boxplotpop13 + rremove("ylab"),
          boxplotpop16,
          boxplotpop18 + rremove("ylab"),
          boxplotpop20 + rremove("ylab"),
          boxplotpop22,
          boxplotpop25 + rremove("ylab"),
          boxplotpop27 + rremove("ylab"),
          ncol=3, nrow=4, labels="auto", hjust=-1, common.legend = TRUE, legend="right")
dev.off()

png("BoxPlots Population Label 4x3.png", units="in", width=15, height=20, res=900)
ggarrange(boxplotpop5,
          boxplotpop6 + theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.title.y=element_blank()),
          boxplotpop10 + theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.title.y=element_blank()),
          boxplotpop11,
          boxplotpop12 + theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.title.y=element_blank()),
          boxplotpop13 + theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.title.y=element_blank()),
          boxplotpop16,
          boxplotpop18 + theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.title.y=element_blank()),
          boxplotpop20 + theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.title.y=element_blank()),
          boxplotpop22,
          boxplotpop25 + theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.title.y=element_blank()),
          boxplotpop27 + theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.title.y=element_blank()),
          ncol=3, nrow=4, labels="auto", hjust=-1, align="hv", common.legend = TRUE, legend="right")
dev.off()


#### Join Geom_Point and BoxPlot by Population ####

# Copy PCA variables in another variable
bypop5 = pca5
bypop6 = pca6
bypop10 = pca10
bypop11 = pca11
bypop12 = pca12
bypop13 = pca13
bypop16 = pca16
bypop18 = pca18
bypop20 = pca20
bypop22 = pca22
bypop25 = pca25
bypop27 = pca27

# Reorder the factors
bypop5$Population = factor(bypop5$Population, levels=c("RABBIT-SLOUGH", "RUS-ASH-GA", "RUS-AN-GA", "RUS-KHA-GA", "RUS-LEV-GA", "NOR-BAR-GA", "NOR-SBJ-GA", "GAC-RUS-PRI", "GAC-NOR-KRI", "GAC-SWE-FIS", "NOR-MYR-GA", "BOOT-LAKE", "BEAR-PAW-LAKE", "NOR-KVA-GA", "NOR-SKF-GA", "FIN-KEV-GA", "RUS-SLI-GA", "NOR-ORR-GA", "ENG-BUT-GA", "GA-POR-VO", "GA-POR-RA", "GA-POR-LI", "GA-POR-TE", "GA-POR-MI", "GA-POR-SA", "MU", "NN", "ITA-STE", "BOS-NER-GA", "MON-SKA-GA", "M", "NB", "KR")) 
bypop6$Population = factor(bypop6$Population, levels=c("RABBIT-SLOUGH", "RUS-ASH-GA", "RUS-AN-GA", "RUS-KHA-GA", "RUS-LEV-GA", "NOR-BAR-GA", "NOR-SBJ-GA", "GAC-RUS-PRI", "GAC-NOR-KRI", "GAC-SWE-FIS", "NOR-MYR-GA", "BOOT-LAKE", "BEAR-PAW-LAKE", "NOR-KVA-GA", "NOR-SKF-GA", "FIN-KEV-GA", "RUS-SLI-GA", "NOR-ORR-GA", "ENG-BUT-GA", "GA-POR-VO", "GA-POR-RA", "GA-POR-LI", "GA-POR-TE", "GA-POR-MI", "GA-POR-SA", "MU", "NN", "ITA-STE", "BOS-NER-GA", "MON-SKA-GA", "M", "NB", "KR")) 
bypop10$Population = factor(bypop10$Population, levels=c("RABBIT-SLOUGH", "RUS-ASH-GA", "RUS-AN-GA", "RUS-KHA-GA", "RUS-LEV-GA", "NOR-BAR-GA", "NOR-SBJ-GA", "GAC-RUS-PRI", "GAC-NOR-KRI", "GAC-SWE-FIS", "NOR-MYR-GA", "BOOT-LAKE", "BEAR-PAW-LAKE", "NOR-KVA-GA", "NOR-SKF-GA", "FIN-KEV-GA", "RUS-SLI-GA", "NOR-ORR-GA", "ENG-BUT-GA", "GA-POR-VO", "GA-POR-RA", "GA-POR-LI", "GA-POR-TE", "GA-POR-MI", "GA-POR-SA", "MU", "NN", "ITA-STE", "BOS-NER-GA", "MON-SKA-GA", "M", "NB", "KR")) 
bypop11$Population = factor(bypop11$Population, levels=c("RABBIT-SLOUGH", "RUS-ASH-GA", "RUS-AN-GA", "RUS-KHA-GA", "RUS-LEV-GA", "NOR-BAR-GA", "NOR-SBJ-GA", "GAC-RUS-PRI", "GAC-NOR-KRI", "GAC-SWE-FIS", "NOR-MYR-GA", "BOOT-LAKE", "BEAR-PAW-LAKE", "NOR-KVA-GA", "NOR-SKF-GA", "FIN-KEV-GA", "RUS-SLI-GA", "NOR-ORR-GA", "ENG-BUT-GA", "GA-POR-VO", "GA-POR-RA", "GA-POR-LI", "GA-POR-TE", "GA-POR-MI", "GA-POR-SA", "MU", "NN", "ITA-STE", "BOS-NER-GA", "MON-SKA-GA", "M", "NB", "KR")) 
bypop12$Population = factor(bypop12$Population, levels=c("RABBIT-SLOUGH", "RUS-ASH-GA", "RUS-AN-GA", "RUS-KHA-GA", "RUS-LEV-GA", "NOR-BAR-GA", "NOR-SBJ-GA", "GAC-RUS-PRI", "GAC-NOR-KRI", "GAC-SWE-FIS", "NOR-MYR-GA", "BOOT-LAKE", "BEAR-PAW-LAKE", "NOR-KVA-GA", "NOR-SKF-GA", "FIN-KEV-GA", "RUS-SLI-GA", "NOR-ORR-GA", "ENG-BUT-GA", "GA-POR-VO", "GA-POR-RA", "GA-POR-LI", "GA-POR-TE", "GA-POR-MI", "GA-POR-SA", "MU", "NN", "ITA-STE", "BOS-NER-GA", "MON-SKA-GA", "M", "NB", "KR")) 
bypop13$Population = factor(bypop13$Population, levels=c("RABBIT-SLOUGH", "RUS-ASH-GA", "RUS-AN-GA", "RUS-KHA-GA", "RUS-LEV-GA", "NOR-BAR-GA", "NOR-SBJ-GA", "GAC-RUS-PRI", "GAC-NOR-KRI", "GAC-SWE-FIS", "NOR-MYR-GA", "BOOT-LAKE", "BEAR-PAW-LAKE", "NOR-KVA-GA", "NOR-SKF-GA", "FIN-KEV-GA", "RUS-SLI-GA", "NOR-ORR-GA", "ENG-BUT-GA", "GA-POR-VO", "GA-POR-RA", "GA-POR-LI", "GA-POR-TE", "GA-POR-MI", "GA-POR-SA", "MU", "NN", "ITA-STE", "BOS-NER-GA", "MON-SKA-GA", "M", "NB", "KR")) 
bypop16$Population = factor(bypop16$Population, levels=c("RABBIT-SLOUGH", "RUS-ASH-GA", "RUS-AN-GA", "RUS-KHA-GA", "RUS-LEV-GA", "NOR-BAR-GA", "NOR-SBJ-GA", "GAC-RUS-PRI", "GAC-NOR-KRI", "GAC-SWE-FIS", "NOR-MYR-GA", "BOOT-LAKE", "BEAR-PAW-LAKE", "NOR-KVA-GA", "NOR-SKF-GA", "FIN-KEV-GA", "RUS-SLI-GA", "NOR-ORR-GA", "ENG-BUT-GA", "GA-POR-VO", "GA-POR-RA", "GA-POR-LI", "GA-POR-TE", "GA-POR-MI", "GA-POR-SA", "MU", "NN", "ITA-STE", "BOS-NER-GA", "MON-SKA-GA", "M", "NB", "KR")) 
bypop18$Population = factor(bypop18$Population, levels=c("RABBIT-SLOUGH", "RUS-ASH-GA", "RUS-AN-GA", "RUS-KHA-GA", "RUS-LEV-GA", "NOR-BAR-GA", "NOR-SBJ-GA", "GAC-RUS-PRI", "GAC-NOR-KRI", "GAC-SWE-FIS", "NOR-MYR-GA", "BOOT-LAKE", "BEAR-PAW-LAKE", "NOR-KVA-GA", "NOR-SKF-GA", "FIN-KEV-GA", "RUS-SLI-GA", "NOR-ORR-GA", "ENG-BUT-GA", "GA-POR-VO", "GA-POR-RA", "GA-POR-LI", "GA-POR-TE", "GA-POR-MI", "GA-POR-SA", "MU", "NN", "ITA-STE", "BOS-NER-GA", "MON-SKA-GA", "M", "NB", "KR")) 
bypop20$Population = factor(bypop20$Population, levels=c("RABBIT-SLOUGH", "RUS-ASH-GA", "RUS-AN-GA", "RUS-KHA-GA", "RUS-LEV-GA", "NOR-BAR-GA", "NOR-SBJ-GA", "GAC-RUS-PRI", "GAC-NOR-KRI", "GAC-SWE-FIS", "NOR-MYR-GA", "BOOT-LAKE", "BEAR-PAW-LAKE", "NOR-KVA-GA", "NOR-SKF-GA", "FIN-KEV-GA", "RUS-SLI-GA", "NOR-ORR-GA", "ENG-BUT-GA", "GA-POR-VO", "GA-POR-RA", "GA-POR-LI", "GA-POR-TE", "GA-POR-MI", "GA-POR-SA", "MU", "NN", "ITA-STE", "BOS-NER-GA", "MON-SKA-GA", "M", "NB", "KR")) 
bypop22$Population = factor(bypop22$Population, levels=c("RABBIT-SLOUGH", "RUS-ASH-GA", "RUS-AN-GA", "RUS-KHA-GA", "RUS-LEV-GA", "NOR-BAR-GA", "NOR-SBJ-GA", "GAC-RUS-PRI", "GAC-NOR-KRI", "GAC-SWE-FIS", "NOR-MYR-GA", "BOOT-LAKE", "BEAR-PAW-LAKE", "NOR-KVA-GA", "NOR-SKF-GA", "FIN-KEV-GA", "RUS-SLI-GA", "NOR-ORR-GA", "ENG-BUT-GA", "GA-POR-VO", "GA-POR-RA", "GA-POR-LI", "GA-POR-TE", "GA-POR-MI", "GA-POR-SA", "MU", "NN", "ITA-STE", "BOS-NER-GA", "MON-SKA-GA", "M", "NB", "KR")) 
bypop25$Population = factor(bypop25$Population, levels=c("RABBIT-SLOUGH", "RUS-ASH-GA", "RUS-AN-GA", "RUS-KHA-GA", "RUS-LEV-GA", "NOR-BAR-GA", "NOR-SBJ-GA", "GAC-RUS-PRI", "GAC-NOR-KRI", "GAC-SWE-FIS", "NOR-MYR-GA", "BOOT-LAKE", "BEAR-PAW-LAKE", "NOR-KVA-GA", "NOR-SKF-GA", "FIN-KEV-GA", "RUS-SLI-GA", "NOR-ORR-GA", "ENG-BUT-GA", "GA-POR-VO", "GA-POR-RA", "GA-POR-LI", "GA-POR-TE", "GA-POR-MI", "GA-POR-SA", "MU", "NN", "ITA-STE", "BOS-NER-GA", "MON-SKA-GA", "M", "NB", "KR")) 
bypop27$Population = factor(bypop27$Population, levels=c("RABBIT-SLOUGH", "RUS-ASH-GA", "RUS-AN-GA", "RUS-KHA-GA", "RUS-LEV-GA", "NOR-BAR-GA", "NOR-SBJ-GA", "GAC-RUS-PRI", "GAC-NOR-KRI", "GAC-SWE-FIS", "NOR-MYR-GA", "BOOT-LAKE", "BEAR-PAW-LAKE", "NOR-KVA-GA", "NOR-SKF-GA", "FIN-KEV-GA", "RUS-SLI-GA", "NOR-ORR-GA", "ENG-BUT-GA", "GA-POR-VO", "GA-POR-RA", "GA-POR-LI", "GA-POR-TE", "GA-POR-MI", "GA-POR-SA", "MU", "NN", "ITA-STE", "BOS-NER-GA", "MON-SKA-GA", "M", "NB", "KR")) 


pointpop5 = ggplot(data=bypop5, aes(x=Population, y=V1, color=Concrete_Location)) +
  geom_boxplot(show_guide=FALSE, outlier.shape=NA) + 
  geom_point(aes(shape=Type),alpha=0.8) +
  ggtitle("Cluster 5") +
  xlab("Population") +
  ylab(paste0("PC1", sep=" (", round((eigen_data5$values[1]/sum(eigen_data5$values)*100),2), "%)")) +
  scale_x_discrete(limits = rev(levels(bypop5$Population))) +
  scale_y_reverse() +
  scale_color_manual(name="Geographic Region", values=c("Adriatic_Sea"="orange", "Alaska"="darkolivegreen4", "North_Scandinavia"="royalblue3", "East_Russia"="tomato4", "English_Channel"="chartreuse3", "Finnish_Golf"="purple3", "Iberian_Peninsula"="red", "North_Sea"="turquoise4"),
                     labels=c("Adriatic_Sea"="Adriatic Sea", "Alaska"="Alaska", "North_Scandinavia"="North Scandinavia", "East_Russia"="East Russia", "English_Channel"="English Channel", "Finnish_Golf"="Gulf of Finland", "Iberian_Peninsula"="Iberian Peninsula", "North_Sea"="North Sea")) +
  scale_shape_manual(values=c("Freshwater"=16, "Marine"=1)) +
  coord_flip() +
  guides(colour=guide_legend(order=1, title="Geographic Region", override.aes=list(size=3, alpha=1)), shape=guide_legend(order=2, title="Habitat Type",override.aes=list(size=3, alpha=1))) +
  theme_bw() +
  theme(plot.title = element_text(hjust=0.5),
        axis.text.y=element_text(color=c("skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "navy", "navy", "navy", "navy", "navy", "navy", "navy", "navy", "navy", "navy", "navy")))

pointpop6 = ggplot(data=bypop6, aes(x=Population, y=V1, color=Concrete_Location)) +
  geom_boxplot(show_guide=FALSE, outlier.shape=NA) + 
  geom_point(aes(shape=Type),alpha=0.8) +
  ggtitle("Cluster 6") +
  xlab("Population") +
  ylab(paste0("PC1", sep=" (", round((eigen_data6$values[1]/sum(eigen_data6$values)*100),2), "%)")) +
  scale_x_discrete(limits = rev(levels(bypop6$Population))) +
  scale_color_manual(name="Geographic Region", values=c("Adriatic_Sea"="orange", "Alaska"="darkolivegreen4", "North_Scandinavia"="royalblue3", "East_Russia"="tomato4", "English_Channel"="chartreuse3", "Finnish_Golf"="purple3", "Iberian_Peninsula"="red", "North_Sea"="turquoise4"),
                     labels=c("Adriatic_Sea"="Adriatic Sea", "Alaska"="Alaska", "North_Scandinavia"="North Scandinavia", "East_Russia"="East Russia", "English_Channel"="English Channel", "Finnish_Golf"="Gulf of Finland", "Iberian_Peninsula"="Iberian Peninsula", "North_Sea"="North Sea")) +
  scale_shape_manual(values=c("Freshwater"=16, "Marine"=1)) +
  coord_flip() +
  guides(colour=guide_legend(order=1, title="Geographic Region", override.aes=list(size=3, alpha=1)), shape=guide_legend(order=2, title="Habitat Type",override.aes=list(size=3, alpha=1))) +
  theme_bw() +
  theme(plot.title = element_text(hjust=0.5),
        axis.text.y=element_text(color=c("skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "navy", "navy", "navy", "navy", "navy", "navy", "navy", "navy", "navy", "navy", "navy")))

pointpop10 = ggplot(data=bypop10, aes(x=Population, y=V1, color=Concrete_Location)) +
  geom_boxplot(show_guide=FALSE, outlier.shape=NA) + 
  geom_point(aes(shape=Type),alpha=0.8) +
  ggtitle("Cluster 10") +
  xlab("Population") +
  ylab(paste0("PC1", sep=" (", round((eigen_data10$values[1]/sum(eigen_data10$values)*100),2), "%)")) +
  scale_x_discrete(limits = rev(levels(bypop10$Population))) +
  scale_y_reverse() +
  scale_color_manual(name="Geographic Region", values=c("Adriatic_Sea"="orange", "Alaska"="darkolivegreen4", "North_Scandinavia"="royalblue3", "East_Russia"="tomato4", "English_Channel"="chartreuse3", "Finnish_Golf"="purple3", "Iberian_Peninsula"="red", "North_Sea"="turquoise4"),
                     labels=c("Adriatic_Sea"="Adriatic Sea", "Alaska"="Alaska", "North_Scandinavia"="North Scandinavia", "East_Russia"="East Russia", "English_Channel"="English Channel", "Finnish_Golf"="Gulf of Finland", "Iberian_Peninsula"="Iberian Peninsula", "North_Sea"="North Sea")) +
  scale_shape_manual(values=c("Freshwater"=16, "Marine"=1)) +
  coord_flip() +
  guides(colour=guide_legend(order=1, title="Geographic Region", override.aes=list(size=3, alpha=1)), shape=guide_legend(order=2, title="Habitat Type",override.aes=list(size=3, alpha=1))) +
  theme_bw() +
  theme(plot.title = element_text(hjust=0.5),
        axis.text.y=element_text(color=c("skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "navy", "navy", "navy", "navy", "navy", "navy", "navy", "navy", "navy", "navy", "navy")))

pointpop11 = ggplot(data=bypop11, aes(x=Population, y=V1, color=Concrete_Location)) +
  geom_boxplot(show_guide=FALSE, outlier.shape=NA) + 
  geom_point(aes(shape=Type),alpha=0.8) +
  ggtitle("Cluster 11") +
  xlab("Population") +
  ylab(paste0("PC1", sep=" (", round((eigen_data11$values[1]/sum(eigen_data11$values)*100),2), "%)")) +
  scale_x_discrete(limits = rev(levels(bypop11$Population))) +
  scale_color_manual(name="Geographic Region", values=c("Adriatic_Sea"="orange", "Alaska"="darkolivegreen4", "North_Scandinavia"="royalblue3", "East_Russia"="tomato4", "English_Channel"="chartreuse3", "Finnish_Golf"="purple3", "Iberian_Peninsula"="red", "North_Sea"="turquoise4"),
                     labels=c("Adriatic_Sea"="Adriatic Sea", "Alaska"="Alaska", "North_Scandinavia"="North Scandinavia", "East_Russia"="East Russia", "English_Channel"="English Channel", "Finnish_Golf"="Gulf of Finland", "Iberian_Peninsula"="Iberian Peninsula", "North_Sea"="North Sea")) +
  scale_shape_manual(values=c("Freshwater"=16, "Marine"=1)) +
  coord_flip() +
  guides(colour=guide_legend(order=1, title="Geographic Region", override.aes=list(size=3, alpha=1)), shape=guide_legend(order=2, title="Habitat Type",override.aes=list(size=3, alpha=1))) +
  theme_bw() +
  theme(plot.title = element_text(hjust=0.5),
        axis.text.y=element_text(color=c("skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "navy", "navy", "navy", "navy", "navy", "navy", "navy", "navy", "navy", "navy", "navy")))

pointpop12 = ggplot(data=bypop12, aes(x=Population, y=V1, color=Concrete_Location)) +
  geom_boxplot(show_guide=FALSE, outlier.shape=NA) + 
  geom_point(aes(shape=Type),alpha=0.8) +
  ggtitle("Cluster 12") +
  xlab("Population") +
  ylab(paste0("PC1", sep=" (", round((eigen_data12$values[1]/sum(eigen_data12$values)*100),2), "%)")) +
  scale_x_discrete(limits = rev(levels(bypop12$Population))) +
  scale_color_manual(name="Geographic Region", values=c("Adriatic_Sea"="orange", "Alaska"="darkolivegreen4", "North_Scandinavia"="royalblue3", "East_Russia"="tomato4", "English_Channel"="chartreuse3", "Finnish_Golf"="purple3", "Iberian_Peninsula"="red", "North_Sea"="turquoise4"),
                     labels=c("Adriatic_Sea"="Adriatic Sea", "Alaska"="Alaska", "North_Scandinavia"="North Scandinavia", "East_Russia"="East Russia", "English_Channel"="English Channel", "Finnish_Golf"="Gulf of Finland", "Iberian_Peninsula"="Iberian Peninsula", "North_Sea"="North Sea")) +
  scale_shape_manual(values=c("Freshwater"=16, "Marine"=1)) +
  coord_flip() +
  guides(colour=guide_legend(order=1, title="Geographic Region", override.aes=list(size=3, alpha=1)), shape=guide_legend(order=2, title="Habitat Type",override.aes=list(size=3, alpha=1))) +
  theme_bw() +
  theme(plot.title = element_text(hjust=0.5),
        axis.text.y=element_text(color=c("skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "navy", "navy", "navy", "navy", "navy", "navy", "navy", "navy", "navy", "navy", "navy")))

pointpop13 = ggplot(data=bypop13, aes(x=Population, y=V1, color=Concrete_Location)) +
  geom_boxplot(show_guide=FALSE, outlier.shape=NA) + 
  geom_point(aes(shape=Type),alpha=0.8) +
  ggtitle("Cluster 13") +
  xlab("Population") +
  ylab(paste0("PC1", sep=" (", round((eigen_data13$values[1]/sum(eigen_data13$values)*100),2), "%)")) +
  scale_x_discrete(limits = rev(levels(bypop13$Population))) +
  scale_color_manual(name="Geographic Region", values=c("Adriatic_Sea"="orange", "Alaska"="darkolivegreen4", "North_Scandinavia"="royalblue3", "East_Russia"="tomato4", "English_Channel"="chartreuse3", "Finnish_Golf"="purple3", "Iberian_Peninsula"="red", "North_Sea"="turquoise4"),
                     labels=c("Adriatic_Sea"="Adriatic Sea", "Alaska"="Alaska", "North_Scandinavia"="North Scandinavia", "East_Russia"="East Russia", "English_Channel"="English Channel", "Finnish_Golf"="Gulf of Finland", "Iberian_Peninsula"="Iberian Peninsula", "North_Sea"="North Sea")) +
  scale_shape_manual(values=c("Freshwater"=16, "Marine"=1)) +
  coord_flip() +
  guides(colour=guide_legend(order=1, title="Geographic Region", override.aes=list(size=3, alpha=1)), shape=guide_legend(order=2, title="Habitat Type",override.aes=list(size=3, alpha=1))) +
  theme_bw() +
  theme(plot.title = element_text(hjust=0.5),
        axis.text.y=element_text(color=c("skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "navy", "navy", "navy", "navy", "navy", "navy", "navy", "navy", "navy", "navy", "navy")))

pointpop16 = ggplot(data=bypop16, aes(x=Population, y=V1, color=Concrete_Location)) +
  geom_boxplot(show_guide=FALSE, outlier.shape=NA) + 
  geom_point(aes(shape=Type),alpha=0.8) +
  ggtitle("Cluster 16") +
  xlab("Population") +
  ylab(paste0("PC1", sep=" (", round((eigen_data16$values[1]/sum(eigen_data16$values)*100),2), "%)")) +
  scale_x_discrete(limits = rev(levels(bypop16$Population))) +
  scale_y_reverse() +
  scale_color_manual(name="Geographic Region", values=c("Adriatic_Sea"="orange", "Alaska"="darkolivegreen4", "North_Scandinavia"="royalblue3", "East_Russia"="tomato4", "English_Channel"="chartreuse3", "Finnish_Golf"="purple3", "Iberian_Peninsula"="red", "North_Sea"="turquoise4"),
                     labels=c("Adriatic_Sea"="Adriatic Sea", "Alaska"="Alaska", "North_Scandinavia"="North Scandinavia", "East_Russia"="East Russia", "English_Channel"="English Channel", "Finnish_Golf"="Gulf of Finland", "Iberian_Peninsula"="Iberian Peninsula", "North_Sea"="North Sea")) +
  scale_shape_manual(values=c("Freshwater"=16, "Marine"=1)) +
  coord_flip() +
  guides(colour=guide_legend(order=1, title="Geographic Region", override.aes=list(size=3, alpha=1)), shape=guide_legend(order=2, title="Habitat Type",override.aes=list(size=3, alpha=1))) +
  theme_bw() +
  theme(plot.title = element_text(hjust=0.5),
        axis.text.y=element_text(color=c("skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "navy", "navy", "navy", "navy", "navy", "navy", "navy", "navy", "navy", "navy", "navy")))

pointpop18 = ggplot(data=bypop18, aes(x=Population, y=V1, color=Concrete_Location)) +
  geom_boxplot(show_guide=FALSE, outlier.shape=NA) + 
  geom_point(aes(shape=Type),alpha=0.8) +
  ggtitle("Cluster 18") +
  xlab("Population") +
  ylab(paste0("PC1", sep=" (", round((eigen_data18$values[1]/sum(eigen_data18$values)*100),2), "%)")) +
  scale_x_discrete(limits = rev(levels(bypop18$Population))) +
  scale_y_reverse() +
  scale_color_manual(name="Geographic Region", values=c("Adriatic_Sea"="orange", "Alaska"="darkolivegreen4", "North_Scandinavia"="royalblue3", "East_Russia"="tomato4", "English_Channel"="chartreuse3", "Finnish_Golf"="purple3", "Iberian_Peninsula"="red", "North_Sea"="turquoise4"),
                     labels=c("Adriatic_Sea"="Adriatic Sea", "Alaska"="Alaska", "North_Scandinavia"="North Scandinavia", "East_Russia"="East Russia", "English_Channel"="English Channel", "Finnish_Golf"="Gulf of Finland", "Iberian_Peninsula"="Iberian Peninsula", "North_Sea"="North Sea")) +
  scale_shape_manual(values=c("Freshwater"=16, "Marine"=1)) +
  coord_flip() +
  guides(colour=guide_legend(order=1, title="Geographic Region", override.aes=list(size=3, alpha=1)), shape=guide_legend(order=2, title="Habitat Type",override.aes=list(size=3, alpha=1))) +
  theme_bw() +
  theme(plot.title = element_text(hjust=0.5),
        axis.text.y=element_text(color=c("skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "navy", "navy", "navy", "navy", "navy", "navy", "navy", "navy", "navy", "navy", "navy")))

pointpop20 = ggplot(data=bypop20, aes(x=Population, y=V1, color=Concrete_Location)) +
  geom_boxplot(show_guide=FALSE, outlier.shape=NA) + 
  geom_point(aes(shape=Type),alpha=0.8) +
  ggtitle("Cluster 20") +
  xlab("Population") +
  ylab(paste0("PC1", sep=" (", round((eigen_data20$values[1]/sum(eigen_data20$values)*100),2), "%)")) +
  scale_x_discrete(limits = rev(levels(bypop20$Population))) +
  scale_color_manual(name="Geographic Region", values=c("Adriatic_Sea"="orange", "Alaska"="darkolivegreen4", "North_Scandinavia"="royalblue3", "East_Russia"="tomato4", "English_Channel"="chartreuse3", "Finnish_Golf"="purple3", "Iberian_Peninsula"="red", "North_Sea"="turquoise4"),
                     labels=c("Adriatic_Sea"="Adriatic Sea", "Alaska"="Alaska", "North_Scandinavia"="North Scandinavia", "East_Russia"="East Russia", "English_Channel"="English Channel", "Finnish_Golf"="Gulf of Finland", "Iberian_Peninsula"="Iberian Peninsula", "North_Sea"="North Sea")) +
  scale_shape_manual(values=c("Freshwater"=16, "Marine"=1)) +
  coord_flip() +
  guides(colour=guide_legend(order=1, title="Geographic Region", override.aes=list(size=3, alpha=1)), shape=guide_legend(order=2, title="Habitat Type",override.aes=list(size=3, alpha=1))) +
  theme_bw() +
  theme(plot.title = element_text(hjust=0.5),
        axis.text.y=element_text(color=c("skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "navy", "navy", "navy", "navy", "navy", "navy", "navy", "navy", "navy", "navy", "navy")))

pointpop22 = ggplot(data=bypop22, aes(x=Population, y=V1, color=Concrete_Location)) +
  geom_boxplot(show_guide=FALSE, outlier.shape=NA) + 
  geom_point(aes(shape=Type),alpha=0.8) +
  ggtitle("Cluster 22") +
  xlab("Population") +
  ylab(paste0("PC1", sep=" (", round((eigen_data22$values[1]/sum(eigen_data22$values)*100),2), "%)")) +
  scale_x_discrete(limits = rev(levels(bypop22$Population))) +
  scale_color_manual(name="Geographic Region", values=c("Adriatic_Sea"="orange", "Alaska"="darkolivegreen4", "North_Scandinavia"="royalblue3", "East_Russia"="tomato4", "English_Channel"="chartreuse3", "Finnish_Golf"="purple3", "Iberian_Peninsula"="red", "North_Sea"="turquoise4"),
                     labels=c("Adriatic_Sea"="Adriatic Sea", "Alaska"="Alaska", "North_Scandinavia"="North Scandinavia", "East_Russia"="East Russia", "English_Channel"="English Channel", "Finnish_Golf"="Gulf of Finland", "Iberian_Peninsula"="Iberian Peninsula", "North_Sea"="North Sea")) +
  scale_shape_manual(values=c("Freshwater"=16, "Marine"=1)) +
  coord_flip() +
  guides(colour=guide_legend(order=1, title="Geographic Region", override.aes=list(size=3, alpha=1)), shape=guide_legend(order=2, title="Habitat Type",override.aes=list(size=3, alpha=1))) +
  theme_bw() +
  theme(plot.title = element_text(hjust=0.5),
        axis.text.y=element_text(color=c("skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "navy", "navy", "navy", "navy", "navy", "navy", "navy", "navy", "navy", "navy", "navy")))

pointpop25 = ggplot(data=bypop25, aes(x=Population, y=V1, color=Concrete_Location)) +
  geom_boxplot(show_guide=FALSE, outlier.shape=NA) + 
  geom_point(aes(shape=Type),alpha=0.8) +
  ggtitle("Cluster 25") +
  xlab("Population") +
  ylab(paste0("PC1", sep=" (", round((eigen_data25$values[1]/sum(eigen_data25$values)*100),2), "%)")) +
  scale_x_discrete(limits = rev(levels(bypop25$Population))) +
  scale_y_reverse() +
  scale_color_manual(name="Geographic Region", values=c("Adriatic_Sea"="orange", "Alaska"="darkolivegreen4", "North_Scandinavia"="royalblue3", "East_Russia"="tomato4", "English_Channel"="chartreuse3", "Finnish_Golf"="purple3", "Iberian_Peninsula"="red", "North_Sea"="turquoise4"),
                     labels=c("Adriatic_Sea"="Adriatic Sea", "Alaska"="Alaska", "North_Scandinavia"="North Scandinavia", "East_Russia"="East Russia", "English_Channel"="English Channel", "Finnish_Golf"="Gulf of Finland", "Iberian_Peninsula"="Iberian Peninsula", "North_Sea"="North Sea")) +
  scale_shape_manual(values=c("Freshwater"=16, "Marine"=1)) +
  coord_flip() +
  guides(colour=guide_legend(order=1, title="Geographic Region", override.aes=list(size=3, alpha=1)), shape=guide_legend(order=2, title="Habitat Type",override.aes=list(size=3, alpha=1))) +
  theme_bw() +
  theme(plot.title = element_text(hjust=0.5),
        axis.text.y=element_text(color=c("skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "navy", "navy", "navy", "navy", "navy", "navy", "navy", "navy", "navy", "navy", "navy")))

pointpop27 = ggplot(data=bypop27, aes(x=Population, y=V1, color=Concrete_Location)) +
  geom_boxplot(show_guide=FALSE, outlier.shape=NA) + 
  geom_point(aes(shape=Type),alpha=0.8) +
  ggtitle("Cluster 27") +
  xlab("Population") +
  ylab(paste0("PC1", sep=" (", round((eigen_data27$values[1]/sum(eigen_data27$values)*100),2), "%)")) +
  scale_x_discrete(limits = rev(levels(bypop27$Population))) +
  scale_color_manual(name="Geographic Region", values=c("Adriatic_Sea"="orange", "Alaska"="darkolivegreen4", "North_Scandinavia"="royalblue3", "East_Russia"="tomato4", "English_Channel"="chartreuse3", "Finnish_Golf"="purple3", "Iberian_Peninsula"="red", "North_Sea"="turquoise4"),
                     labels=c("Adriatic_Sea"="Adriatic Sea", "Alaska"="Alaska", "North_Scandinavia"="North Scandinavia", "East_Russia"="East Russia", "English_Channel"="English Channel", "Finnish_Golf"="Gulf of Finland", "Iberian_Peninsula"="Iberian Peninsula", "North_Sea"="North Sea")) +
  scale_shape_manual(values=c("Freshwater"=16, "Marine"=1)) +
  coord_flip() +
  guides(colour=guide_legend(order=1, title="Geographic Region", override.aes=list(size=3, alpha=1)), shape=guide_legend(order=2, title="Habitat Type",override.aes=list(size=3, alpha=1))) +
  theme_bw() +
  theme(plot.title = element_text(hjust=0.5),
        axis.text.y=element_text(color=c("skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "skyblue1", "navy", "navy", "navy", "navy", "navy", "navy", "navy", "navy", "navy", "navy", "navy")))


png("PointPlots Population Multiple Labels 4x3.png", units="in", width=15, height=20, res=900)
ggarrange(pointpop5,
          pointpop6 + rremove("ylab"),
          pointpop10 + rremove("ylab"),
          pointpop11,
          pointpop12 + rremove("ylab"),
          pointpop13 + rremove("ylab"),
          pointpop16,
          pointpop18 + rremove("ylab"),
          pointpop20 + rremove("ylab"),
          pointpop22,
          pointpop25 + rremove("ylab"),
          pointpop27 + rremove("ylab"),
          ncol=3, nrow=4, labels="auto", hjust=-1, common.legend = TRUE, legend="right")
dev.off()

png("PointPlots Population Label 4x3.png", units="in", width=15, height=20, res=900)
ggarrange(pointpop5,
          pointpop6 + theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.title.y=element_blank()),
          pointpop10 + theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.title.y=element_blank()),
          pointpop11,
          pointpop12 + theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.title.y=element_blank()),
          pointpop13 + theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.title.y=element_blank()),
          pointpop16,
          pointpop18 + theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.title.y=element_blank()),
          pointpop20 + theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.title.y=element_blank()),
          pointpop22,
          pointpop25 + theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.title.y=element_blank()),
          pointpop27 + theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.title.y=element_blank()),
          ncol=3, nrow=4, labels="auto", hjust=-1, align="hv", common.legend = TRUE, legend="right")
dev.off()
