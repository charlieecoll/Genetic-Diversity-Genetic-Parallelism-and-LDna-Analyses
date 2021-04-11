##################################################################
##                          Fst Matrix                          ##
##################################################################

## Carla Coll Costa, C3. January 2021.
## This script takes a .gen file and calculates several statistics from this file. 
## Among these statistics, there is pairwise Fst for each pair of populations,
## which is plotted and saved automatically.


# Upload the required libraries
library(diveRsity)
library(reshape2)
library(ggplot2)
library(gdata)
library(easyGgplot2)

# Set the working directory
setwd("<Working_Directory>")

# Calculate Gst (Nei), G'st (Hedrick), Theta (Weir & Cockerman), D (Jost) and Fstatistics (Weir & Cockerman)
diffCalc(infile="<Input_GEN_file>", outfile="<Output_Folder_Name>", fst=TRUE, pairwise=TRUE)

# Read the Fst pairwise matrix and add the upper triangle values
fst_pairwise = read.table("<Input_File_with_Pairwise_Infromation_Calculated_with_diffCalc>", header=TRUE, sep="\t", row.names=1)
values = lowerTriangle(fst_pairwise, diag=FALSE)
upperTriangle(fst_pairwise, diag=FALSE, byrow=TRUE) = values

# Get the population names
names = row.names(fst_pairwise)

# Set the column names as the population names
colnames(fst_pairwise) = names

# Set the data as a matrix
fst_pairwise = as.matrix(fst_pairwise)

# Melt the matrix so it can be plotted with ggplot2
fst_pairwise_melted = melt(data=fst_pairwise)

# Order by geographic region and habitat type
fst_pairwise_melted$Var1 = factor(fst_pairwise_melted$Var1, levels=c("BEAR-PAW-LAKE","BOOT-LAKE","RABBIT-SLOUGH","RUS-AN-GA","RUS-ASH-GA","RUS-KHA-GA","FIN-KEV-GA","NOR-SKF-GA","NOR-KVA-GA","NOR-SBJ-GA","RUS-LEV-GA","RUS-SLI-GA","GAC-RUS-PRI","NOR-ORR-GA","NOR-MYR-GA","GA-POR-LI","GA-POR-MI","GA-POR-RA","GA-POR-SA","GA-POR-TE","GA-POR-VO","BOS-NER-GA","KR","M","MON-SKA-GA","MU","NB","NN"))
fst_pairwise_melted$Var2 = factor(fst_pairwise_melted$Var2, levels=c("BEAR-PAW-LAKE","BOOT-LAKE","RABBIT-SLOUGH","RUS-AN-GA","RUS-ASH-GA","RUS-KHA-GA","FIN-KEV-GA","NOR-SKF-GA","NOR-KVA-GA","NOR-SBJ-GA","RUS-LEV-GA","RUS-SLI-GA","GAC-RUS-PRI","NOR-ORR-GA","NOR-MYR-GA","GA-POR-LI","GA-POR-MI","GA-POR-RA","GA-POR-SA","GA-POR-TE","GA-POR-VO","BOS-NER-GA","KR","M","MON-SKA-GA","MU","NB","NN"))

# Final Plot
png("Fst Matrix Final.png", units="in", width=8, height=8, res=900)
ggplot(data=fst_pairwise_melted, aes(x=Var1, y=Var2, fill=value)) +
  geom_tile(color="black") +
  scale_fill_gradientn(colours=c("blue", "yellow", "red"), na.value="white") +
  scale_x_discrete(labels=c("BEAR-PAW-LAKE","BOOT-LAKE","RABBIT-SLOUGH*","RUS-AN-GA*","RUS-ASH-GA*","RUS-KHA-GA*","FIN-KEV-GA","NOR-SKF-GA","NOR-KVA-GA","NOR-SBJ-GA*","RUS-LEV-GA*","RUS-SLI-GA","GAC-RUS-PRI*","NOR-ORR-GA","NOR-MYR-GA*","GA-POR-LI","GA-POR-MI","GA-POR-RA","GA-POR-SA","GA-POR-TE","GA-POR-VO","BOS-NER-GA","KR","M","MON-SKA-GA","MU","NB","NN")) +
  scale_y_discrete(labels=c("BEAR-PAW-LAKE","BOOT-LAKE","RABBIT-SLOUGH*","RUS-AN-GA*","RUS-ASH-GA*","RUS-KHA-GA*","FIN-KEV-GA","NOR-SKF-GA","NOR-KVA-GA","NOR-SBJ-GA*","RUS-LEV-GA*","RUS-SLI-GA","GAC-RUS-PRI*","NOR-ORR-GA","NOR-MYR-GA*","GA-POR-LI","GA-POR-MI","GA-POR-RA","GA-POR-SA","GA-POR-TE","GA-POR-VO","BOS-NER-GA","KR","M","MON-SKA-GA","MU","NB","NN")) +
  labs(x="Population", y="Population") +
  theme(axis.text.x=element_text(angle=90, vjust=0.5),
        axis.text=element_text(color=c("darkolivegreen4","darkolivegreen4","darkolivegreen4","tomato4","tomato4","tomato4","royalblue3","royalblue3","royalblue3","royalblue3","royalblue3","purple3","purple3","turquoise4","turquoise4","red","red","red","red","red","red","orange","orange","orange","orange","orange","orange","orange")),
        plot.title=element_text(hjust=0.5),
        panel.background=element_blank())
dev.off()