# Genetic-Diversity-Genetic-Parallelism-and-LDna-Analyses

This repository contains the pipeline used for an analysis of genetic diversity, genetic parallelism and LD of three-spine stickleback populations using RAD data (.sh and .R files). It can be used to conduct similar analysis on other species. Folders and files are ordered in the needed order to conduct the analysis.

Preparing the Data: Scripts to remove the PCR duplicates and map and align the raw data to a reference genome.
Genotype Likelihoods: Scripts to get the site frequency spectrum (SFS) of the individuals and populations (1D and 2D) in case of low coverage of the data.
Hard-Call Data: Scripts to hard-call the SNP variants.
Genetic Diversity: Scripts to conduct genetic diversity analyses on the data (nucleotide diversity, Watterson's theta, Tajima's D, runs of homozygosity and inbreeding based on runs of homozygosity (Froh)).
Global Parallelism: Scripts to conduct analysis of global parallelism through PCA using genotype likelihoods.
Population Structure: Scripts to conduct genetic differentiation (Fst), admixture and PCA analyses.
Linkage Disequilibrium (LD): Scripts to conduct studies of LD by plotting LD heatmaps and performing LDna analysis.
