#################################################################
##                             ROH                             ##
#################################################################

## Carla Coll Costa, C3. January 2021.
## This script takes the files on the ROH data calculated in BCFTOOLS with the script 
## 5_ROH. It also needs a map file and a genotype file for the ROH data calculated. 


# Upload the required packaged
library(detectRUNS)
library(stringr)
library(Eagle)
library(easyGgplot2)
library(plyr)

# Set the working directory
setwd("<Working_Directory>")

# Take the data
data = readExternalRuns(inputFile="<Input_File_with_ROH_Data>", program="BCFtools")

# Remove the scaffolds
data = data[- grep("scaffold", data$chrom),]

# Get the map path file
map = file.path("<Input_Map_File>")

# Get the genotype file
genotype=file.path("~/Desktop/Master's Thesis/Results/17. ROH/0_data/Genotype_File.ped")

# We need to add the population to the group column
data$group = data$id
data$group = sub("-[^-]+$", "", data$group)

# We substitute the chromosome name to an appropriate name
data$chrom = str_replace_all(data$chrom,"groupI\\b","1")
data$chrom = str_replace_all(data$chrom,"groupII\\b","2")
data$chrom = str_replace_all(data$chrom,"groupIII\\b","3")
data$chrom = str_replace_all(data$chrom,"groupIV\\b","4")
data$chrom = str_replace_all(data$chrom,"groupV\\b","5")
data$chrom = str_replace_all(data$chrom,"groupVI\\b","6")
data$chrom = str_replace_all(data$chrom,"groupVII\\b","7")
data$chrom = str_replace_all(data$chrom,"groupVIII\\b","8")
data$chrom = str_replace_all(data$chrom,"groupIX\\b","9")
data$chrom = str_replace_all(data$chrom,"groupX\\b","10")
data$chrom = str_replace_all(data$chrom,"groupXI\\b","11")
data$chrom = str_replace_all(data$chrom,"groupXII\\b","12")
data$chrom = str_replace_all(data$chrom,"groupXIII\\b","13")
data$chrom = str_replace_all(data$chrom,"groupXIV\\b","14")
data$chrom = str_replace_all(data$chrom,"groupXV\\b","15")
data$chrom = str_replace_all(data$chrom,"groupXVI\\b","16")
data$chrom = str_replace_all(data$chrom,"groupXVII\\b","17")
data$chrom = str_replace_all(data$chrom,"groupXVIII\\b","18")
data$chrom = str_replace_all(data$chrom,"groupXIX\\b","19")
data$chrom = str_replace_all(data$chrom,"groupXX\\b","20")
data$chrom = str_replace_all(data$chrom,"groupXXI\\b","21")

# Get information on the populations and individuals
pop_labs = read.table("<File_with_Information_on_Population_and_Individuals>", header=TRUE)

#### Plot NROH against SROH ####

# Sum ROH genome
sum_ROH_genome = ddply(data,.(id),summarize,sum=sum(lengthBps)/10^6)

# Sum ROH for sample
count_ROH_genome = count(data,"id")
sum_ROH_genome=merge(sum_ROH_genome,count_ROH_genome,by='id')
sum_ROH_genome=merge(sum_ROH_genome,unique(data[,c("id","group")]),by='id')

# Merge other information to plot latter
sum_ROH_genome = merge(sum_ROH_genome, pop_labs, by.x="id", by.y="Individual")

# Delete duplicated column
sum_ROH_genome$group = NULL

# Final plot
png("NROHvsSROH.png", units="in", width=8, height=6, res=900)
ggplot(data=sum_ROH_genome, aes(x=sum, y=freq, colour=Concrete_Location)) +
  geom_point(size=2, aes(shape=Type), stroke=1, alpha=0.8) +
  geom_text(data=subset(sum_ROH_genome, sum > 150 | freq > 200),
            aes(sum,freq,label=Population), size=3,
            nudge_x=rep((5), times=3), 
            nudge_y=rep((5), times=3)) + 
  xlab("Sum of ROH in Mbps") +
  ylab ("Number of ROH for Individual") +
  scale_shape_manual(values=c("Freshwater"=16, "Marine"=1)) +
  scale_color_manual(values=c("Adriatic"="orange", "Alaska"="darkolivegreen4", "North_Scandinavia"="royalblue3", "East_Rusia"="tomato4", "English_Channel"="chartreuse3", "Finnish_Golf"="purple3", "Iberian_Peninsula"="red", "North_Sea"="turquoise4"),
                     labels=c("Adriatic"="Adriatic Sea", "Alaska"="Alaska", "East_Rusia"="East Russia", "English_Channel"="English Channel", "Finnish_Golf"="Finnish Golf", "Iberian_Peninsula"="Iberian Peninsula", "North_Sea"="North Sea", "North_Scandinavia"="North Scandinavia")) +
  guides(colour=guide_legend(order=1, title="Geographic Region", override.aes=list(size=3, alpha=1)), shape=guide_legend(order=2, title="Habitat Type",override.aes=list(size=3, alpha=1)))+
  theme_bw()
dev.off()
