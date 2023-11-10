library(vegan)
library(tidyverse)
library(phyloseq)
library(magrittr)
library(here)

#' Read in phyloseq object
phy.clean <- readRDS(here("output", "compiled", "phy.clean.rds"))

#' Remove low abundance OTUs by prevalence and Hellinger transformation
(phy.clean %<>% phyloseq::filter_taxa(., function(x){sum(x>0) > 0.01*nsamples(.)}, TRUE)) 
phy.clean <-  transform_sample_counts(phy.clean, function(x) (x / sum(x)))

otu_test <- phy.clean %>% otu_table %>% data.frame

#' Subsetting data
phy.w = subset_samples(phy.clean, Region=="West")
phy.e = subset_samples(phy.clean, Region=="East")

phy.corv = subset_samples(phy.clean, Garden=="Corvallis")
phy.herm = subset_samples(phy.clean, Garden=="Hermiston")

phy.corv_w = subset_samples(phy.clean, Garden=="Corvallis" & Region=="West" )
phy.herm_e = subset_samples(phy.clean, Garden=="Hermiston" & Region=="East")
phy.corv_e = subset_samples(phy.clean, Garden=="Corvallis" & Region=="East")
phy.herm_w = subset_samples(phy.clean, Garden=="Hermiston" & Region=="West")

adonis_phy <- function(phy.in, x){
  X <- phy.in %>% sample_data %>% data.frame %>% dplyr::select(all_of=x) %>% .[,1]
  Y <- phy.in %>% phyloseq::distance("bray")
  adonis2(Y~X, permutations=9999, parallel=parallel::detectCores())
}


#### perMANOVA functions (Devin function from EricoidChronosequence) ####
# .[,1] a df to vector

# 1 predictor 
adonis_phy <- function(phy.in, x){
  X <- phy.in %>% sample_data %>% data.frame %>% dplyr::select(all_of=x) %>% .[,1]
  Y <- phy.in %>% phyloseq::distance("bray")
  adonis2(Y~X, permutations=9999, parallel=parallel::detectCores())
}

# interaction
adonis_phy_interaction <- function(phy.in,x,z){
  X <- phy.in %>% sample_data %>% data.frame %>% select(all_of=x) %>% .[,1]
  Z <- phy.in %>% sample_data %>% data.frame %>% select(all_of=z) %>% .[,1]
  Y <- phy.in %>% otu_table %>% data.frame %>% 
    philentropy::distance('jensen-shannon', est.prob = "empirical", use.row.names = T, as.dist.obj = T) %>% 
    sqrt
  adonis2(Y ~ X*Z, permutations=9999, parallel=parallel::detectCores())
}

#### Block effect ####
phy.corv_block <- adonis_phy(phy.corv,"Block") # 29.5% 
phy.herm_block <- adonis_phy(phy.herm,"Block") # 25.2%

#### Genotype, Garden, and GxE effects ####
phy_all <- adonis_phy_interaction(phy.clean, "Genotype", "Garden")
phy_all
# Genotype = .16
# Garden = .21
# GxE = .07

#### Genotype effect within garden ####
phy.corv_adj <- adonis_phy(phy.corv,"Genotype") # 29.5% 
phy.herm_adj <- adonis_phy(phy.herm,"Genotype") # 26.8%

#### Genotype effect within garden and ecotype ####
phy.corv_w_adj <- adonis_phy(phy.corv_w,"Genotype") # 18.1%
phy.herm_e_adj <- adonis_phy(phy.herm_e,"Genotype") # 19.9%
phy.corv_e_adj<- adonis_phy(phy.corv_e,"Genotype") # 29.8%
phy.herm_w_adj <- adonis_phy(phy.herm_w,"Genotype") # 21.4%
