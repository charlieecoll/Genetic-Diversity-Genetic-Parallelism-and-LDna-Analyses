##################################################################
##                   Comparison of Statistics                   ##
##################################################################

## Carla Coll Costa, C3. January 2021.
## This script compares nucleotide diversity calculated from the population 1D SFS
## across samples, grouping them in different geographical regions. It tests if the
## differences between these groups are statistically significant. The input data consists
## on tab separated files or excel files in which each line corresponds to one population
## and the columns have different information.
## The scripts tests the assumptions of the ANOVA test and then performs such test and a
## Tukey-Kramer post-hoc analysis to see which pairs of comparisons are significantly different.


# Upload the required packages
library(ggplot2)

# Set the directory
setwd("<Working_Directory>")

# Load the data on nucleotide diversity
div_stats = read.table("<File_with_Summary_Statistics>", header=TRUE)

# Sort dataframe by population name
div_stats = div_stats[order(div_stats$POP),]

# Change name of population column
colnames(div_stats)[1] = "Population"

# load data with information on each population
pop_list = read.table("<File_with_Information_on_the_Populations>", header=TRUE)

# Fuse dataframes
div_stats = merge(div_stats, pop_list, by="Population")

# Remove groups with less than three observations: English_Channel
div_stats_concrete_location = div_stats[!(div_stats$Concrete_Location=="English_Channel"),]

# Check the assumptions of the ANOVA test for the data on nucleotide diversity. These assumptions are:
  # Normality: Each sample is taken from a normally distributed population.
  # Homoscedasticity: Variance Equality.
  # Sample Independence.
  # Continuous dependent variable.
# We assume sample independence since a ranzomized design was used when collecting the samples.
# The dependent variable, nucleotide diversity, is continuous.

# Test normality
shapiro.test(div_stats_concrete_location$pi)

# Test variance equality
bartlett.test(pi~Concrete_Location, data=div_stats_concrete_location)

# ANOVA to test nucleotide diversity across geographic regions
anovapi_concrete_location = aov(data=div_stats_concrete_location, formula=pi~Concrete_Location)
summary(anovapi_concrete_location)

# Optional: Check assumptions of ANOVA test with plots
plot(anovapi_concrete_location)

# Post-hoc analysis with a tukey-kramer test
tukeypi = TukeyHSD(x=anovapi_concrete_location, conf.level=0.95)
tukeypi