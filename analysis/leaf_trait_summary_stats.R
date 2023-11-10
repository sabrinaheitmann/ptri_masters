#' mean +- std error is made into a table in EWG/output/leaf_trait_summary_stats.xlx

library(magrittr)
library(tidyverse)
library(here)
library(ggthemes)
library(psych)

traits <- read.csv(here("data", "gardenData", "trait_plus_meta_csvs", "trait_wmeta_allblocks.csv"))

traits %<>%
  dplyr::select(Garden, Region, Genotype, thick, SLA, LDMC, stomata_length, stomata_count, percent_C_weight,percent_N_weight)

#' Rename traits
traits %<>% 
  rename(LT = thick) %>%
  rename(SL = stomata_length) %>%
  rename(SD = stomata_count) %>%
  rename(N = percent_N_weight) %>%
  rename(C = percent_C_weight)

#' Tidy columns 
new <- traits %>% pivot_longer(
  cols = LT:N,
  names_to = "trait",
  values_to = "measure"
) %>%
  drop_na()

#### INTERACTION PLOT ####

#' Standard error function
se <- function(x) sqrt(var(x) / length(x))

#' Find mean and standard error of traits for each ecotype within each garden
means <- new %>% group_by(Garden, Region,trait) %>%  
  summarise(se = se(measure),
            mean = mean(measure))

#' Plot seperately for each leaf trait
pal.region <- c("#fb832d", "#016392")

ggplot(means, aes(x=Garden, y=mean, group=Region, colour = Region)) + 
  scale_colour_manual(values = pal.region) +
  facet_wrap(~trait, scales = "free", ncol = 2) +
  geom_line() +
  geom_point()+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2,
                position=position_dodge(0.05)) + 
  labs(color='Host Ecotype')  +
  theme_bw() +
  theme(axis.title.x=element_blank(),
        axis.title.y = element_blank(),
        legend.title=element_text(size=20),
        legend.text =element_text(size=18),
        strip.text = element_text(size=20),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15))
  
ggsave("output/figs/figure4.png", width = 8, height = 6, units = "in")
