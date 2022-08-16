library(tidyverse)
library(microbiome)
library(microbiomeutilities)
library(phyloseq)
library(ggplot2)
library(network)
library(igraph)
library(intergraph)
library(GGally)
library(ggnet)
library(ggnetwork)
library(dplyr) # data handling
library(SpiecEasi) # Network analysis for sparse compositional data  
library(devtools)

ps_3_clean

ps_D<-subset_samples(ps_3_clean, sample_data(ps_3_clean)$Caste_composite == "T1 colony Worker")
ps_D_200<-subset_taxa(ps_D, taxa_sums(ps_D)>200)
 
otu_D<-otu_table(ps_D_200)

set.seed(123)
net_D<-spiec.easi(otu_D, "mb", lambda.min.ratio = 1e-3, nlambda = 300)
net_D_beta<-symBeta(getOptBeta(net_D))
net_D_refit <- getRefit(net_D)