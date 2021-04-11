##################################################################
##                   AVERAGE COVERAGE BOXPLOT                   ##
##################################################################

## Carla Coll Costa, C3. September 2020.
## This script plots the average coverage or depth of sequencing for different populations
## in the same plot. It expects a tab separated file with 3 columns:
##  1. Population name
##  2. Individual name
##  3. Average depth of each individual
## The plot generated is a boxplot with the average coverage per population. 


# Upload the libraries required
library(ggplot2)
library(dplyr)
library("readxl")

# Upload the data
data = read.delim("<Tab_File>", header=FALSE)
type = data.frame(read_excel("<Population_Information>"))

# Plot the average coverage
ggplot(data=data, aes(x=data$V1, y=data$V3)) +
  geom_boxplot() +
  labs(x="Population", y="Average Coverage") +
  ggtitle ("Average Coverage per Population") +
  theme (axis.text.x=element_text(angle=90),
         plot.title=element_text(hjust=0.5))

