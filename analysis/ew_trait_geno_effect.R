library(magrittr)
library(tidyverse)
library(here)
library(ggthemes)

pal.region <- c("#fb832d", "#016392")

geno_effect <- read_csv(here("newphy", "ew_geno_effect.csv"))

geno_effect %<>% select(-1)
#from spring_anlyses/permanovas.R
microbiome_geno_wcorv <- c(.072, "Corvallis", "microbiome", "West", 1)
microbiome_geno_ecorv <- c(.155, "Corvallis", "microbiome", "East", 1)
microbiome_geno_wherm <- c(.035, "Hermiston", "microbiome", "West", 1)
microbiome_geno_eherm <- c(.076, "Hermiston", "microbiome", "East", 1)

geno_effect <- rbind(geno_effect, microbiome_geno_wcorv, microbiome_geno_ecorv, microbiome_geno_wherm, microbiome_geno_eherm)
geno_effect$geno_effect <- as.numeric(geno_effect$geno_effect)

#geno_effect %<>% filter(significance == 1)

geno_effect <- geno_effect %>%
  dplyr::mutate(
    label = case_when(
      significance == 1 ~ "*",
      significance == 0 ~ "",
      TRUE ~ NA_character_
    )
  )

leaftraits <- c("Stomata Length", "Stomata Density", "Leaf thickness", "Specific leaf \narea (SLA)", "Leaf Dry Matter \nContent (LDMC)", "Leaf \nMicrobiome")

#adds significance asterisks to figure based on bar position
dodger = position_dodge(width = 0.9)
ggplot(geno_effect, aes(x = geno_effect, y = reorder(measurement, -geno_effect), fill = ecotype)) +
  geom_bar(position=dodger,stat="identity") +
  coord_flip() + 
  facet_wrap(~garden) +
  scale_fill_manual(values=pal.region) + 
  theme_clean() +
  scale_y_discrete(labels = leaftraits) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.title.x = element_blank()) +
  xlab("Genotype effect") +
  geom_text(aes(label=label, group = ecotype), position = dodger) 
