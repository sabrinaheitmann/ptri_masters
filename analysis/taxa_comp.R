#10/18/22 - to identify fungal community composition at different taxonomic levels (beginning of MS results)

library(vegan)
library(tidyverse)
library(phyloseq)
library(magrittr)
library(here)
library(ggthemes)
library(mvabund)
library(foreach)
library(MASS)
library(ggtext)
library(ggnewscale)
library(ggvegan)
library(cowplot)
library(scico)
library(patchwork)

source("code/ewg_functions.R")

phy <- readRDS(here("output", "compiled", "phy.clean.rds"))
#(phy %<>% phyloseq::filter_taxa(., function(x){sum(x>0) > 0.01*nsamples(.)}, TRUE)) 

otu_table <- phy %>% otu_table %>% data.frame



#PHYLUM
phylum <- phy %>% tax_glom("Phylum") #This collapses everything to the phylum level first.
tax<- tax_table(phylum) %>% data.frame #create tax table
phylum_otu <- otu_table(phylum)%>% data.frame #create OTU table

phylum_otu %<>% t() %>% data.frame #flip matrix and convert back to data frame
phylum_otu$OTUSums <- rowSums(phylum_otu) #get sums of read counts for each OTU

otunames <- row.names(phylum_otu)#extract names
Sums <- phylum_otu$OTUSums #extract sums

df <- data.frame(otunames, Sums) #create data frame

df %<>%
  mutate(sum = sum(Sums)) %>%
  mutate(percent = (Sums/sum)*100)

df$tax <- tax$Phylum
df

#CLASS
class <- phy %>% tax_glom("Class") #This collapses everything to the class level first.
tax<- tax_table(class) %>% data.frame #create tax table
class_otu <- otu_table(class)%>% data.frame #create OTU table

class_otu %<>% t() %>% data.frame #flip matrix and convert back to data frame
class_otu$OTUSums <- rowSums(class_otu) #get sums of read counts for each OTU

otunames <- row.names(class_otu)#extract names
Sums <- class_otu$OTUSums #extract sums

df <- data.frame(otunames, Sums) #create data frame

df %<>%
  mutate(sum = sum(Sums)) %>%
  mutate(percent = (Sums/sum)*100)

df$tax <- tax$Class
df

#OTU
genus <- phy %>% tax_glom("Genus") #This collapses everything to the genus level first.
tax<- tax_table(genus) %>% data.frame #create tax table
genus_otu <- otu_table(genus)%>% data.frame #create OTU table

genus_otu %<>% t() %>% data.frame #flip matrix and convert back to data frame
genus_otu$OTUSums <- rowSums(genus_otu) #get sums of read counts for each OTU

otunames <- row.names(genus_otu)#extract names
Sums <- genus_otu$OTUSums #extract sums

df <- data.frame(otunames, Sums) #create data frame

df %<>%
  mutate(sum = sum(Sums)) %>%
  mutate(percent = (Sums/sum)*100)

df$tax <- tax$Genus
df

phy.clean <- readRDS(here("output", "compiled", "phy.clean.rds"))
#clean_tax <- phy.clean %>% tax_table %>% data.frame
phy.clean %<>% phyloseq::filter_taxa(., function(x){sum(x>0) > 0.01*nsamples(.)}, TRUE)
phy.clean <-  transform_sample_counts(phy.clean, function(x) (x / sum(x)))

phy.rust = subset_taxa(phy.clean, Genus=="g__Melampsora")
phy.mycosph = subset_taxa(phy.clean, Genus=="g__Mycosphaerella")

phy.clean_otu <- phy.clean %>% otu_table %>% data.frame

phy.clean_otu %<>%
  rownames_to_column() %>%
  dplyr::select(F.1, F.2, rowname) %>%
  mutate(add = F.1 + F.2)
#mutate(over1 = ifelse(add > 1, "TRUE", add))

meta_data <- phy.clean %>% sample_data %>% data.frame %>% rownames_to_column() %>% dplyr::select(Region, Garden, rowname)

otu_abund <- left_join(phy.clean_otu, meta_data, by = "rowname")

otu_abund <- otu_abund %>%
  dplyr::select(Garden, Region, F.1, F.2) %>%
  group_by(Garden,Region) %>%
  summarise(across(everything(), ~ mean(.x, na.rm = TRUE)))


##### By garden #####

#PHYLUM
phy.garden <- phy %>% merge_samples("Garden")
phylum.garden <- phy.garden %>% tax_glom("Phylum") #This collapses everything to the phylum level first.

tax<- tax_table(phylum.garden) %>% data.frame #create tax table
phylum_otu <- phylum.garden %>% otu_table %>% data.frame #create OTU table
phylum_otu %<>% t() %>% data.frame #flip matrix and convert back to data frame

phylum_otu %<>% rownames_to_column("OTU")
df.long <- pivot_longer(phylum_otu, cols = 2:3, names_to = "Garden", values_to = "count")

df.long %<>%
  group_by(OTU) %>%
  mutate(sum = sum(count)) %>%
  mutate(percent = (count/sum)*100)

tax %<>% rownames_to_column("OTU") %>% dplyr::select(Phylum, OTU)
phylum <- left_join(tax, df.long)

#Classes
#phy.garden <- phy %>% merge_samples("Garden")

class.garden <- phy.garden %>% tax_glom("Class") #This collapses everything to the class level first.

tax<- tax_table(class.garden) %>% data.frame #create tax table
class_otu <- class.garden %>% otu_table %>% data.frame #create OTU table
class_otu %<>% t() %>% data.frame #flip matrix and convert back to data frame

class_otu %<>% rownames_to_column("OTU")
df.long <- pivot_longer(class_otu, cols = 2:3, names_to = "Garden", values_to = "count")

df.long %<>%
  #group_by(OTU) %>%
  mutate(sum = sum(count)) %>%
  mutate(percent = (count/sum)*100)

tax %<>% rownames_to_column("OTU") %>% dplyr::select(Class, OTU)
class <- left_join(tax, df.long)

#### GENUS ####

phy.corv = subset_samples(phy, Garden=="Corvallis")
phy.herm = subset_samples(phy, Garden=="Hermiston")

genus <- phy.herm %>% tax_glom("Genus") #This collapses everything to the Genus level first.
tax<- tax_table(genus) %>% data.frame #create tax table
genus_otu <- otu_table(genus)%>% data.frame #create OTU table

genus_otu %<>% t() %>% data.frame #flip matrix and convert back to data frame
genus_otu$OTUSums <- rowSums(genus_otu) #get sums of read counts for each OTU

otunames <- row.names(genus_otu)#extract names
Sums <- genus_otu$OTUSums #extract sums

df <- data.frame(otunames, Sums) #create data frame

df %<>%
  mutate(sum = sum(Sums)) %>%
  mutate(percent = (Sums/sum)*100)

df$tax <- tax$genus
df
