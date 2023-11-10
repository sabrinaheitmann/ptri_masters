library(vegan)
library(tidyverse)
library(phyloseq)
library(magrittr)
library(here)
library(ggthemes)
library(foreach)
library(MASS)
library(ggtext)
library(ggnewscale)
library(ggvegan)
library(cowplot)
library(scico)
library(patchwork)

pal.region <- c("#fb832d", "#016392")

#' Read in phyloseq object
phy.clean <- readRDS(here("output", "compiled", "phy.clean.rds"))

#' Remove low abundance OTUs by prevalence
(phy.clean %<>% phyloseq::filter_taxa(., function(x){sum(x>0) > 0.01*nsamples(.)}, TRUE)) 

#' Rel abund
phy.clean <-  transform_sample_counts(phy.clean, function(x) (x / sum(x)))

#' Edit row names
meta_data <- phy.clean %>% sample_data %>% data.frame

meta_gr <- meta_data %>%  
  mutate(garden_genotype = paste0(Garden, "_" , Genotype)) %>%
  mutate(garden_population = paste0(Garden, "_" , Population)) %>%
  mutate(garden_region = paste0(Garden, "_" , Region)) 

meta_gr %<>% 
  rownames_to_column() %>%
  mutate(SampID=substring(Sample_ID,10)) %>%
  mutate(SampID = rowname) %>%
  column_to_rownames("SampID") %>%
  sample_data

sample_data(phy.clean) <- meta_gr
meta_data <- phy.clean %>% sample_data %>% data.frame
otu <- phy.clean %>% otu_table() %>% data.frame

#' Make garden_genotype (i.e. genotype within garden), garden_population, garden_region columns
df.tidy <- phy.clean %>%
  otu_table %>% data.frame %>% 
  bind_cols(sample_data(phy.clean) %>% data.frame %>% 
              dplyr::select(Garden, Region, Genotype, Population, rowname)) %>%
  pivot_longer(all_of(taxa_names(phy.clean)),names_to="Taxa",values_to="abundance") %>%
  mutate(ID=paste(rowname,Genotype),
         garden_genotype = paste(Garden, Genotype, sep = "_"),
         garden_region = paste(Garden, Region, sep = "_"),
         garden_population = paste(Garden, Population, sep = "_"),
         Region=ifelse(Region=="East","E","W") %>%
           factor(levels=c("W","E")))

#' Add mean abundance to garden_region, garden_genotype, and garden_population
df.1c <- df.tidy %>% 
  group_by(Region,Genotype,Garden, garden_region, garden_genotype, garden_population, Population, rowname,Taxa) %>%
  summarise(abundance=mean(abundance)) %>%
  ungroup %>%
  pivot_wider(names_from = Taxa, values_from=abundance)
set.seed(99)
JSDdist <- df.1c %>% select_if(is.numeric) %>%
  #wisconsin %>% 
  as.matrix %>% 
  #philentropy::JSD(est.prob = "empirical") %>% sqrt
  vegdist("bray")

# fit dbRDA constrained by treatment, conditioned on genotype
(dbrda.treatment <- dbrda(JSDdist~Region*Garden, data=df.1c))

# get variance explained by the ordination axes
dbrda.treatment.r2 <- (eigenvals(dbrda.treatment)/sum(eigenvals(dbrda.treatment))) %>% multiply_by(100) %>% round(1)

# ordination 
fortify(dbrda.treatment,display="wa") %>%
  filter(score=="sites") %>% 
  bind_cols(df.1c) %>%
  group_by(garden_region) %>%
  summarise(n=n(),x=mean(dbRDA1),y=mean(dbRDA2),sd1=sd(dbRDA1),sd2=sd(dbRDA2)) %>%
  left_join(sample_data(phy.clean) %>%
              data.frame %>%
              dplyr::select(Region, Garden, garden_region), 
            by = "garden_region" ) %>%
  ggplot(aes(x=x,y=y,fill = Region, shape = Garden)) +
  geom_errorbar(aes(ymin=y-sd2,ymax=y+sd2),width=0, color = "darkgrey")+
  geom_errorbarh(aes(xmin=x-sd1,xmax=x+sd1),height=0, color = "darkgrey")+
  geom_point(color = 'black', size = 6) +
  #coord_fixed(eigenvals(dbrda.treatment)[2]/eigenvals(dbrda.treatment)[1])+
  labs(x=paste0("dbRDA1 [",dbrda.treatment.r2[1],"%]"),
       y=paste0("dbRDA2 [",dbrda.treatment.r2[2],"%]")) +
  #scale_shape_manual("Garden", values = c(16, 17, 0, 1), labels = c('Corvallis', 'Hermiston')) +
  #guides(colour = guide_legend(override.aes = list(size=5))) +
  theme_few()+
  theme(legend.position = 'bottom', 
    legend.margin = margin(0,.05,.05,.05, unit="cm"),
    legend.box = 'vertical', 
    legend.text = element_text(size = 20),
    legend.title = element_text(size = 20, face = "bold"),
    axis.title = element_text(size=18),
    axis.text = element_text(size=15),
    panel.border = element_rect(fill=NA, colour = "black", size=1.5)) +
  scale_shape_manual("Garden", values = c(21L, 23L), labels = c('Corvallis', 'Hermiston')) +
  scale_fill_manual("Host Ecotype", values = pal.region, labels = c("East", "West")) +
  guides(fill = guide_legend(override.aes = list(shape = 21)))

ggsave("output/figs/figure2.png")
