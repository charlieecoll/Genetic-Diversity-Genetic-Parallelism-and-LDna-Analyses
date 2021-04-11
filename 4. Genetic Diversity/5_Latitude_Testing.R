#################################################################
##             TEST GENETIC DIVERSITY VS. LATITUDE             ##
#################################################################

## Carla Coll Costa, C3. January 2021.
## This script has two parts: the first part tests genetic diversity statistics (nucleotide diversity and Froh) 
## against latitude and habitat type at a population level using linear models to see if the variables are 
## good predictors of the statistics. It expects .txt files with the data for each of the populations.
## The second part tests genetic diversity statistics (nucleotide diversity and Froh) against latitude and 
## habitat type at an individual level using population as a random factor in the  linear models.
## It expects .txt files with the data for each of the populations.

#### PART 1 ####
#### Prepare the Data ####

# Upload packages
library(readxl)
library(ggplot2)
library(scales)
library(lme4)
library(lmtest)
library(olsrr)
library(lattice)
library(ICC)
library(rsq)

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
data_stats$Colour_Concrete_Location = NULL
data_stats$Colour_Location = NULL
data_stats$Colour_Type = NULL

# Remove duplicate rows
data_stats = unique(data_stats)

#### Can we test Habitat Type and Latitude Simultaneously? ####

# ANOVA one-way: Check for correlation between habitat type and latitude (Significance Test)
correlation = aov(data_stats$GPS_N~data_stats$Type.y)
summary(correlation)

# Should the correlation worry us? We calculate the intraclass correlation coefficient (Effect Size)
ICCbare(x=data_stats$Type.y, y=data_stats$GPS_N)

#### Nucleotide Diversity LINEAR MODEL: Test Habitat Type and Latitude Simultaneously ####

# We apply a linear model using habitat type and latitude as predictors
linearmodel1 = lm(pi~GPS_N+Type.y, data=data_stats)
anova(linearmodel1)
summary(linearmodel1)

# We check if the assumptions of generalized linear models are met:
# 1. Linearity of residuals.
# 2. Independence of residuals.
# 3. Normal distribution of residuals.
# 4. Equal variance of residuals.

# Plot diagnostic plots:
par(mfrow = c(2, 2), mar=c(5.1, 4.1, 4.1, 2.1))
plot(linearmodel1)
hist(resid(linearmodel1))
# Residuals vs Fitted: Linearity of residuals and homoscedasticity.
# Normal Q-Q: Normality of residuals.
# Scale-Location: Homoscedasticity
# Residuals vs Leverage: Identify influential cases

# Residual histogram
par(mfrow = c(1, 1))
ols_plot_resid_hist(linearmodel1)

# Check for normal distribution of residuals with different tests
ols_test_normality(linearmodel1)
shapiro.test(resid(linearmodel1))

# Check for equal variance of residuals with Breusch-Pagan test
bptest(linearmodel1)


#### Mean Froh LINEAR MODEL: Test Habitat Type and Latitude Simultaneously ####

# We apply a linear model using habitat type and latitude as predictors
linearmodel2 = lm(froh.Froh_genome~GPS_N+Type.y, data=data_stats)
anova(linearmodel2)
summary(linearmodel2)

# We check if the assumptions of generalized linear models are met:
# 1. Linearity of residuals.
# 2. Independence of residuals.
# 3. Normal distribution of residuals.
# 4. Equal variance of residuals.

# Plot diagnostic plots:
par(mfrow = c(2, 2), mar=c(5.1, 4.1, 4.1, 2.1))
plot(linearmodel2)
# Residuals vs Fitted: Linearity of residuals and homoscedasticity.
# Normal Q-Q: Normality of residuals.
# Scale-Location: Homoscedasticity
# Residuals vs Leverage: Identify influential cases

# Residual histogram
par(mfrow = c(1, 1))
ols_plot_resid_hist(linearmodel2)

# Check for normal distribution of residuals with different tests
ols_test_normality(linearmodel2)
shapiro.test(resid(linearmodel2))

# Check for equal variance of residuals with Breusch-Pagan test
bptest(linearmodel4)

## TRANSFORM THE DATA ##

# LOG #
log1 = log(data_stats$froh.Froh_genome)
linearmodellog1 = lm(log1~GPS_N+Type.y, data=data_stats)
summary(linearmodellog1)
# Check for the normality assumption with the Shapiro-Wilk test
ols_test_normality(linearmodellog1)
shapiro.test(resid(linearmodellog1))

# Check for variance equality with the Barlett test
bptest(linearmodellog2)

# GAMMA #
linearmodelgamma = glm(froh.Froh_genome~GPS_N+Type.y, data=data_stats, family=Gamma(link=log))
summary(linearmodelgamma)
anova(linearmodelgamma, test="F")

par(mfrow = c(2, 2), mar=c(5.1, 4.1, 4.1, 2.1))
plot(linearmodelgamma)

# Check for the normality assumption with the Shapiro-Wilk test
shapiro.test(resid(linearmodelgamma))

# Check for variance equality with the Barlett test
bptest(linearmodelgamma)

# Calculate R^2 adjusted values
rsq(linearmodelgamma,adj=TRUE,type="v")

#### PART 2 ####
#### Prepare the Data ####

# Upload packages
library(readxl)
library(ggplot2)
library(scales)
library(lme4)
library(lmtest)
library(olsrr)
library(lattice)
library(car)
library(ICC)

# Set working directory
setwd("<Working_Directory>")

# Get the data
het = read.table("<File_with_Individual_Heterozygosity>", header=TRUE)
froh = read.table("<File_with_Froh_for_each_Individual>", header=TRUE)
latitude = data.frame(read_excel("<File_with_Individual_Latitude_Coordinates>"))

# Order by population name
het = het[order(het$individuals),]
froh = froh[order(froh$id),]
latitude = latitude[order(latitude$POP_ID),]

# Delete not needed columns
froh$sum = NULL
froh$Location = NULL
froh$Colour_Location = NULL
froh$Colour_Concrete_Location = NULL
froh$Colour_Type = NULL

# Fuse data frames for heterozygosity and froh
data_stats = merge(het, froh, by.x="individuals", by.y="id", all=TRUE)

# Delete not needed columns
data_stats$Population.y = NULL
data_stats$Type.y = NULL
data_stats$Concrete_Location.y = NULL

# Fuse data frames to get latitude information
data_stats = merge(data_stats, latitude, by.x="Population.x", by.y="POP_ID", all=TRUE)

# Delete not needed columns
data_stats$Type = NULL
data_stats$N = NULL

#### Can we test Habitat Type and Latitude Simultaneously? ####

# ANOVA one-way: Check for correlation between habitat type and latitude (Significance Test)
correlation = aov(data_stats$GPS_N~data_stats$Type.x)
summary(correlation)

# Should the correlation worry us? We calculate the intraclass correlation (Effect Size)
ICCbare(x=data_stats$Type.x, y=data_stats$GPS_N)

#### Test Heterozygosity: Habitat Type and Latitude Simultaneously ####

# We apply a linear models with and without random factor and compare them
lm3 = lm(heterozygosities~GPS_N+Type.x, data=data_stats)
lmm1 = lmer(heterozygosities~GPS_N+Type.x+(1|Population.x), data=data_stats)
anova(lmm1, lm3)

# Look at the summary
summary(lmm1)

# Here we don't obtain the p-values, we need to calculate them
betas = fixef(lmm1) # Take the estimates of the slopes from the summary
se = sqrt(diag(vcov(lmm1))) # Calculate the SE
pval = 2*pnorm(-abs(betas/se)) # We divide the slope by the SE to know how many SD away from 0 are we
output = cbind(betas, se, pval)
print(output, digist=3)

# We check if the assumptions of generalized linear mixed models are met:
# 1. Linearity of residuals.
# 2. Independence of residuals.
# 3. Normal distribution of residuals.
# 4. Equal variance of residuals.

# Plot diagnostic plots:
residuals_lmm1 = resid(lmm1)
predicted_lmm1= fitted.values(lmm1)
par(mfrow = c(1, 1))
hist(residuals_lmm1)
qqmath(lmm1, id=0.05)
plot(lmm1)
abline(a=0, b=0)

# Check for normal distribution of residuals with different tests
shapiro.test(resid(lmm1))

# Check for equal variance of residuals with Bartlett test
bartlett.test(resid(lmm1)~data_stats$GPS_N)

#### Test Froh: Habitat Type and Latitude Simultaneously ####

# We apply a linear models with and without random factor and compare them
lm4 = lm(Froh_genome~GPS_N+Type.x, data=data_stats)
lmm2 = lmer(Froh_genome~GPS_N+Type.x+(1|Population.x), data=data_stats)
anova(lmm2, lm4)

# Look at the summary
summary(lmm2)

# Here we don't obtain the p-values, we need to calculate them
betas = fixef(lmm2) # Take the estimates of the slopes from the summary
se = sqrt(diag(vcov(lmm2))) # Calculate the SE
pval = 2*pnorm(-abs(betas/se)) # We divide the slope by the SE to know how many SD away from 0 are we
output = cbind(betas, se, pval)
print(output, digist=3)

# We check if the assumptions of generalized linear mixed models are met:
# 1. Linearity of residuals.
# 2. Independence of residuals.
# 3. Normal distribution of residuals.
# 4. Equal variance of residuals.

# Plot diagnostic plots:
residuals_lmm2 = resid(lmm2)
predicted_lmm2 = fitted.values(lmm2)
par(mfrow = c(1, 1))
hist(residuals_lmm2)
qqmath(lmm2, id=0.05)
plot(lmm2)
abline(a=0, b=0)

# Check for normal distribution of residuals with different tests
shapiro.test(resid(lmm2))

# Check for equal variance of residuals with Bartlett test
bartlett.test(resid(lmm2)~data_stats$GPS_N)

