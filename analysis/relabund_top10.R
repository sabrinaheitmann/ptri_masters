#CLR after top 10

library(bipartite)
library(vegan)
library(tidyverse)
library(phyloseq)
library(magrittr)
library(here)
library(lme4)
library(broom)
library(emmeans)
library(ggthemes)
library(microbiome)

pal.region <- c("#fb832d", "#016392")

source("code/ewg_functions.R")
#' Read in phyloseq object
phy.clean <- readRDS(here("output", "compiled", "phy.clean.rds"))

#' Remove low abundance OTUs by prevalence
(phy.clean %<>% phyloseq::filter_taxa(., function(x){sum(x>0) > 0.01*nsamples(.)}, TRUE)) 
#' relativize
phy.clean <-  transform_sample_counts(phy.clean, function(x) (x / sum(x)))

tb <- psmelt(phy.clean) %>%
  as_tibble

#' find top 10 OTUs by relative abundace
tb %>%
  group_by(OTU) %>%
  summarise(Mean=mean(Abundance)) %>%
  unique() %>%
  top_n(10) %>%
  arrange(desc(Mean)) #F.1, 2, 4, 3,6, 7, 5, 8, 9, 10 were top otus

#' Center log ratio of top 10 dominant OTUs by relative abundance
phy.clean <-  microbiome::transform(phy.clean, "clr")
tb <- psmelt(phy.clean) %>%
  as_tibble


#' Select top 10 OTUs
re_mod <- tb %>%
  filter(OTU == "F.1" | OTU == "F.2" | OTU == "F.4" | OTU == "F.3" | OTU == "F.6" | OTU == "F.7" | OTU == "F.5" | OTU == "F.8" | OTU == "F.9" | OTU == "F.10" ) 

re_mod %>% dplyr::select(OTU, Genus) %>% unique()

out <- split(re_mod, f = re_mod$OTU)

#' Find 95% confidence intervals for each garden-ecotype combination
mylmer <- function(df){
  fit1 <- lm(Abundance ~ Garden*Genotype, data = df)
  emm1 <- emmeans(fit1, ~Garden*Genotype)
  emmdf1 = as.data.frame(emm1)
  
  region <- test %>% dplyr::select(Region, Genotype) %>% unique()
  
  emmdf1 <- left_join(emmdf1, region)
  
  emmdf1 %<>% 
    group_by(Garden, Region) %>%
    summarise(emmean = mean(emmean),
              SE = mean(SE),
              lower.CL = mean(lower.CL),
              upper.CL = mean(upper.CL)) %>%
    unite(Garden_Region, Garden, Region, remove = FALSE) 
}

list <- out %>% lapply(mylmer)

fig_df <- reshape2::melt(list)

fig_df <- pivot_wider(fig_df, names_from = variable)

fig_df %<>% 
  mutate(across(L1, factor, levels=c("F.1","F.2","F.4", "F.3","F.6","F.7","F.5","F.8","F.9", "F.10"))) %>%
  mutate(across(Garden_Region, factor, levels=c("Corvallis_West", "Corvallis_East", "Hermiston_West", "Hermiston_East")))

new <- c("Cladosporium", "Melampsora", "Alternaria.1", "Didymellaceae", "Ramularia", "Alternaria.2", "Tilletiopsis", "Neocamarosporium", "Phaeosphaeria", "Epicoccum")

names(new) <- c("F.1","F.2","F.4", "F.3","F.6","F.7","F.5","F.8","F.9", "F.10")

#### Plot Figure S1 ####
ggplot(fig_df, aes(x = Garden_Region, y = emmean)) +
  geom_point(aes(fill = Region, shape = Garden), size = 3) +
  geom_linerange(aes(ymin=lower.CL, ymax=upper.CL), width=.2,
                 position=position_dodge(0.05)) +
  facet_wrap(~L1, nrow = 5, labeller = labeller(L1 = new))+
  scale_colour_manual(values = c("#fb832d", "#016392")) +
  #guides(fill = guide_legend(override.aes = list(size=0.8)))+
  theme_few()+
  coord_flip() +
  #scale_fill_manual("Host genotype",values=pal.ecotype)+
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text = element_text(size=14,face="italic"),
        axis.title.x = element_text(size=14),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12),
        legend.margin=margin(6,6,6,12),
        plot.tag = element_text(angle=90,hjust=0),
        plot.tag.position = c(0.848,0.115)) +
  scale_shape_manual("Garden", values = c(21L, 23L), labels = c('Corvallis', 'Hermiston')) +
  scale_fill_manual("Host Ecotype", values = c("#fb832d", "#016392"), labels = c("East", "West")) +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  labs(y = "Centered log-ratio of relative abundance")

ggsave("output/figs/relabund_v2.tiff")

#### Genotype-by-garden effect ANOVA ####
anova <- function(df){
  fit1 <- aov(Abundance ~ Garden*Genotype, data = df)
  summary(fit1) 
  fit1df <- tidy(fit1)
}

anova_list <- out %>% lapply(anova)

anova_df <- reshape2::melt(anova_list)

anova_df <- pivot_wider(anova_df, names_from = variable)

#anova_df %<>% filter(term != "Residuals")
#anova_df %>%
#  filter(term == "Garden:Genotype") %>%
#  filter(p.value < 0.05)

write_csv(anova_df, "output/tables/top_otu_anova_tableS1.csv")
