library(magrittr)
library(tidyverse)
library(here)
library(ggthemes)

#trait genotype effects from leaf_trait_data_analysis.R)
geno_effect <- read_csv(here("output", "csvs", "ew_geno_effect.csv"))

geno_effect %<>% select(-1)

#microbiome genotype effects from permanovas.R
#these are the genotype effects in my masters thesis, not the updated genotype effects reflected in permanovas.R
microbiome_geno_wcorv <- c(.176, "Corvallis", "microbiome", "West", 1)
microbiome_geno_ecorv <- c(.282, "Corvallis", "microbiome", "East", 1)
microbiome_geno_wherm <- c(.192, "Hermiston", "microbiome", "West", 1)
microbiome_geno_eherm <- c(.199, "Hermiston", "microbiome", "East", 1)

geno_effect <- rbind(geno_effect, microbiome_geno_wcorv, microbiome_geno_ecorv, microbiome_geno_wherm, microbiome_geno_eherm)
geno_effect$geno_effect <- as.numeric(geno_effect$geno_effect)

#adds new column with asterisks to represent significant genotype effect
geno_effect <- geno_effect %>%
  dplyr::mutate(
    label = case_when(
      significance == 1 ~ "*",
      significance == 0 ~ "",
      TRUE ~ NA_character_
    )
  )

#### Genotype effect plot ####

#Leaf trait abbrevs
leaftraits <- c("SL", "SD", "LT", "SLA", "LDMC", "Microbiome")
#adds significance asterisks to figure based on bar position
dodger = position_dodge(width = 0.9)
#color palette 
pal.region <- c("#016392", "#fb832d")
#orders facets from West to East
geno_effect$ecotype = factor(geno_effect$ecotype, levels=c('West','East'))

ggplot(geno_effect, aes(x = geno_effect, y = reorder(measurement, -geno_effect), fill = garden)) +
  geom_bar(position=dodger,stat="identity") +
  coord_flip() + 
  facet_wrap(~ecotype) +
  scale_fill_manual(values=pal.region) + 
  theme_clean() +
  scale_y_discrete(labels = leaftraits) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.title.x = element_blank()) +
  xlab("R\u00b2 genotype") +
  geom_text(aes(label=label, group = garden), position = dodger) 

ggsave("output/figs/figureS2.png")
