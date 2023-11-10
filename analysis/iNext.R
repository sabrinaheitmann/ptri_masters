library(vegan)
library(tidyverse)
library(phyloseq)
library(magrittr)
library(here)
library(lme4)
library(ggthemes)
library(ggpubr)
library(iNEXT)

##### UNHASH TO CHANGE iNEXT PARAMETERS (OUTPUT FOR CODE IS IN CSV BELOW ) ####
# # remove empty OTUs from phyloseq object
# sweepOTUs <- function(phy){
#   prune_taxa(taxa_sums(phy) > 0,phy)
# }
# 
# phy.clean <- readRDS(here("output", "compiled", "phy.clean.rds")) %>%
#   sweepOTUs
# 
# # Get estimated diversity and richness estimates from phyloseq object with iNEXT
#   iNEXT.phy <- function(phy.in,level, use.cores=parallel::detectCores()){
#   require(foreach)
#   require(doMC)
#   registerDoMC(cores=use.cores)
#   otuTab <- phy.in %>% otu_table %>% t %>% data.frame 
#   out <- foreach(i=1:ncol(otuTab), .combine=bind_rows) %dopar% {
#     estimateD(otuTab[,i],level=level) %>% 
#       tibble %>%
#       mutate(SampID=colnames(otuTab)[i])
#   }
#   out %>% 
#     filter(order!=2) %>%
#     transmute(SampID=SampID,
#               Metric=case_when(order==0 ~ 'Species richness',
#                                order==1 ~ 'Shannon diversity'),
#               Estimate=qD,
#               LCL=qD.LCL,
#               UCL=qD.UCL) %>%
#     left_join(sample_data(phy.in) %>% data.frame %>%
#                 rownames_to_column("SampID"))
# }
# 
# div_dat <- iNEXT.phy(phy.clean,level=1000)
# 
# rich <- div_dat %>% filter(Metric == "Species richness")
# write_csv(rich, "output/csvs/inext_rich.csv")
# div <- div_dat %>% filter(Metric == "Shannon diversity")
# write_csv(div, "output/csvs/inext_div.csv")

#' Format rich and div data 
rich <- read_csv("output/csvs/inext_rich.csv")
rich = mutate(rich, Region_Garden = paste(Region, "_", Garden, sep = ""))
div <- read_csv("output/csvs/inext_div.csv")
div = mutate(div, Region_Garden = paste(Region, "_", Garden, sep = ""))

#' Calculate mean diversity and richness
mean_div <- div %>%
  group_by(Garden, Genotype, Region) %>%
  summarise("Shannon Diversity" = mean(Estimate))

mean_rich <- rich %>%
  group_by(Garden, Genotype, Region) %>%
  summarise(Richness = mean(Estimate))

#' Combine dfs
both <- left_join(mean_div,mean_rich)

both_long <- both %>% pivot_longer(cols = c("Shannon Diversity", Richness),
             names_to = "metric",
             values_to = "measure")


#### FIGURE 3 - Plot richness and diversity boxplots with averaged values by genotype_garden ####
pal.region <- c("#fb832d", "#016392")

p <- ggplot(both_long, aes(x = Garden, y = measure, color = Region)) +
  geom_boxplot(alpha = 0.1, width=0.75) +
  geom_point(position = position_jitterdodge(jitter.width=0.15), size = .5) +
  facet_wrap(~metric, scales = "free") +
  scale_color_manual(values = pal.region) +
  theme_clean() +
  theme(axis.title.x=element_blank(),
        axis.ticks.x=element_blank(), 
        axis.title.y=element_blank(),
        axis.ticks.y=element_blank(), 
        axis.text.x = element_text(size = 18),
        strip.text = element_text(size=20),
        legend.title = element_text(size=20),
        legend.text = element_text(size = 20),
        axis.text.y = element_text(size = 18)) +
  #ylab("Genotype mean diversity values") +
  labs(color='Host \necotype') 

ggsave('output/figs/figure3.png')

##### FIGURE 3 - Ecotype x E ANOVA #####

#' Fungal rich ANOVA (Region*Garden)
rich.aov<- lm(Estimate ~ Region*Garden, data = rich)
summary(rich.aov)$r.squared
GxR_int_rich <- aov(rich.aov) %>% broom::tidy()
GxR_int_rich

#' rich tukey
aov(rich.aov) %>% TukeyHSD

#' Fungal diversity ANOVA (Genotype*Garden)
#div.aov<- lm(Estimate ~ Genotype*Garden, data = div)
#summary(div.aov)$r.squared
#GxR_int_div <- aov(div.aov) %>% broom::tidy()

#' div tukey
#GxR_int_div
#aov(div.aov) %>% TukeyHSD

##### MERGE DIVERSITY DATA WITH LEAF TRAIT DATA #####

leaf_trait_data <- read.csv(here("data", "gardenData", "trait_plus_meta_csvs", "trait_wmeta_allblocks.csv"))

leaf_trait_data %<>%
  dplyr::select(Sample_ID, Garden, Genotype, Region, Population, thick, SLA, LDMC, stomata_length, stomata_count, percent_C_weight,percent_N_weight)

leaf_trait_data %<>% 
  rename(LT = thick) %>%
  rename(SL = stomata_length) %>%
  rename(SD = stomata_count) %>%
  rename(N = percent_N_weight) %>%
  rename(C = percent_C_weight)

ltd_avg <- leaf_trait_data %>%
  dplyr::select(Garden, Genotype, LT:N) %>%
  group_by(Garden, Genotype) %>%
  summarise(across(everything(), ~ mean(.x, na.rm = TRUE))) %>%
  na.omit()

#' Add ltd metadata (Region, Population)
ltd_avg <- left_join(ltd_avg, both, by = c("Genotype", "Garden")) %>% unique()

##### FIGURE 7 - PLOT DIVERSITY MEASURES AGAINST STOMATAL LENGTH #####
trait_diversity_plots <- function(trait, measure, trait_name, measure_name){
  ggplot(ltd_avg, aes(trait, measure)) +
    geom_point(aes(colour=Region), size = 2.5) +
    #geom_smooth(method="lm", aes(color=Region), se = F, alpha = .3) +
    facet_wrap(~Garden) +
    scale_color_manual(values = pal.region) +
    theme(strip.text = element_text(size=20),
          axis.title = element_text(size=18),
          legend.title=element_text(size=20),
          legend.text =element_text(size=18)) +
     # stat_cor(
     #   aes(label = paste(..rr.label..,
     #                     if_else(readr::parse_number(..p.label..) < 0.001, 
     #                             "p<0.001", ..p.label..), sep = "~`,`~")), size = 4, vjust = 1, hjust = -.2) +
    labs(color='Host Ecotype')  +
    theme_hc() +
    xlab(trait_name) + 
    ylab(measure_name)
}

SL_a <- trait_diversity_plots(ltd_avg$SL, ltd_avg$`Shannon Diversity`, "Stomatal length (pixels)", "Shannon Diversity")
SL_b <- trait_diversity_plots(ltd_avg$SL, ltd_avg$Richness, "Stomatal length (pixels)", "Richness")
ggarrange(SL_a + rremove("xlab"), SL_b, nrow = 2, common.legend = TRUE, legend = "right", align = "hv")
ggsave("output/figs/figure7.png")
