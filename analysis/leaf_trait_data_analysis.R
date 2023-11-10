library(magrittr)
library(tidyverse)
library(here)
library(ggthemes)
library(partR2)
library(ggstatsplot)
library(broom)
library(data.table)
#source("code/ewg_functions.R")

#' Upload full metadata and leaf trait trees (Corvallis: B2, B6, B9, B12 & Hermiston: B4, B7, B11, B12)
all_tree_meta_data <- read_csv(here("data", "gardenData", "all_tree_meta_data.csv"))
meta <- read_csv(here("data", "gardenData", "Merged_meta&sampledata.csv"))
trait_meta <- read_csv(here("data", "gardenData", "leaftraits", "leafmeta.csv"))

#' Editing all_tree_meta_data to only include trait blocks
all_tree_meta_data = subset(all_tree_meta_data, Garden == "Corvallis" & Block == "2" | Garden == "Corvallis" & Block == "6" | Garden == "Corvallis" & Block == "9" | Garden == "Corvallis" & Block == "12" | Garden == "Hermiston" & Block == "4" | Garden == "Hermiston" & Block == "7" | Garden == "Hermiston" & Block == "11" | Garden == "Hermiston" & Block == "12")

all_tree_meta_data %<>% dplyr::select(-c(Longitude, Latitude, Comment, Pop.1, Number, Collection))

all_tree_meta_data <- all_tree_meta_data[-c(1)]

#' Upload leaf trait data
wet_thick <- read.csv(here("data", "gardenData", "leaftraits", "wet_thick_2020.csv"), header = TRUE, na.strings=c(""," ","NA"))
dry <- read_csv(here("data", "gardenData", "leaftraits", "dry_2020.csv"))
area <- read_csv(here("data", "gardenData", "leaftraits", "area_2020.csv"))
stomata <- read_csv(here("data", "gardenData", "leaftraits", "stomata_2020_data.csv"))
stomata_2021 <- read.csv(here("data", "gardenData", "leaftraits", "stomata_measurements_2021.csv"))
cn_data <- read_csv(here("data", "gardenData", "leaftraits", "elemental_analysis_full_norepeats.csv"))

#### LEAF TRAIT TIDYING ####

#' Find any mistypes
summary(area) # none
summary(wet_thick) # none
summary(dry)

#' Obvious mistypes were removed
dry$dry_L4[dry$dry_L4==0] <- NA
dry$dry_L2[dry$dry_L2==0.005] <- NA

#### Leaf wet weight and thickness tidying #### 
#' Combine wet and thick duplicated genotypes (this is due to sampling tagged leaves on one day and untagged leaves the next day) 
#' Edited the column names and removed any rows that had N/A's for all leaves

#' Separate column, row, and block indicators into separate columns, while also saving as ID
wet_thick %<>%
  separate(Tree_ID, c("C", "R", "Block")) %>%
  dplyr::select(-Garden)

wet_thick = mutate(wet_thick,
                   ID = paste(C, R, Block, sep = ""))

#' Remove trees that do not have any data associated with them
wet_thick <- wet_thick[!rowSums(is.na(wet_thick)) == 12, ]

#' These are all the samples that measurements were taken on less than 4 leaves then new untagged leaves were added to also measure
#' Remove them from merged datatable
n_occur <- data.frame(table(wet_thick$ID))
duplicate <- n_occur[n_occur$Freq > 1,]
removed <- wet_thick[!(wet_thick$ID %in% duplicate$Var1 ),]
added_untagged <- wet_thick[wet_thick$ID %in% n_occur$Var1[n_occur$Freq > 1],]

#' This removes all of the duplicates that weren't exactly the same when combined
setDT(added_untagged)[, lapply(.SD, na.omit), by = ID]
error <- setDT(added_untagged)[, code := NULL][, lapply(.SD, na.omit), by = ID]   
error.keep <- error[duplicated(error), ]

#' Combine the corrected duplicates in merged datatable
wet_thick <- rbind(removed, error.keep)

wet_thick %>% count(No_tag) #108 samples included 1 or more untagged leaves
wet_thick %<>% dplyr::select(-No_tag)

#' Pivot longer
wet_thick %<>% 
  dplyr::select(-c(C, R, Block)) %>%
  pivot_longer(- ID, 
               names_to = "measurement", 
               values_to = "obs") %>%
  separate(measurement, c("measurement", "leaf", "side"))


#### Leaf dry weight tidying ####
dry %<>% 
  dplyr::select(-(Garden)) %>%
  pivot_longer(- ID, 
               names_to = "measurement", 
               values_to = "obs") %>%
  separate(measurement, c("measurement", "leaf"))

#### Leaf area tidying ####
area %<>%
  separate(ID, c("C","R", "Block", "Garden", "leaf", "tif"), sep = "([.])") %>%
  unite("ID", "C", "R", "Block", sep = "") %>%
  dplyr::select(-c(Garden, tif))

area %<>% 
  pivot_longer(c(area, length), 
               names_to = "measurement", 
               values_to = "obs") 

#### Stomata length and number tidying ####
stomata_2021 %<>%
  separate(Label, c("Garden", "C", "R", "Block", "Image", "Numbers"), sep = "([.])") %>%
  unite("ID", "C", "R", "Block", sep = "") %>%
  dplyr::select(-c(Numbers, Garden))

stomata_2020 <- stomata %>%
  separate(Label, c("Garden", "C", "R", "Block", "Side", "Image", "tif"), sep = "([.])") %>%
  unite("ID", "C", "R", "Block", sep = "") %>%
  dplyr::select(-c(tif, Side, Garden))

stomata_full <- rbind(stomata_2020, stomata_2021)

#' There were 3 mistypes (greater than 80 so delete)
summary(stomata_full)
stomata_full <- stomata_full[stomata_full$Length < 80, ]

#' Found mean stomata length/ID
stomata_mean <- stomata_full %>%
  group_by(ID, Image) %>%
  summarise(stomata_length = mean(Length)) %>%
  group_by(ID) %>%
  summarise(stomata_length = mean(stomata_length))

#' Found mean stomata count/ID
stomata_count <- stomata_full %>%
  group_by(ID, Image) %>%
  count(name = "stomata_count") %>%
  group_by(ID) %>%
  summarise(stomata_count = mean(stomata_count))

stomata_joined <- left_join(stomata_mean, stomata_count, by.x = "ID")

#### C and N data tidying ####
cn_data %<>%
  mutate(CtoN = percent_C_weight/percent_N_weight) %>%
  #dplyr::select(-Notes, -percent_C_weight, -percent_N_weight)
  dplyr::select(-Notes)

cn_key2 <- all_tree_meta_data %>%
  dplyr::select(ID, Genotype, Garden)

cn_data <- left_join(cn_data, cn_key2)
cn_data %<>% dplyr::select(-c(Genotype, Garden)) 


#### BIND ALL LEAF TRAIT DATA ####

#' bind area,  wet weight,  dry weight, thickness
traits2 <- bind_rows(area, wet_thick, dry)

#' Find trait mean for each sample and pivot wide 
mean_traits2 <- traits2 %>%
  group_by(ID,measurement) %>%
  summarise(mean_bysample = mean(obs, na.rm = TRUE))

mean_traits2 <- pivot_wider(mean_traits2, names_from = measurement, values_from = mean_bysample)

#' Change units then calculate useful leaf trait measurements
mean_traits2 %<>%
  mutate(area = area*10) %>% # cm^2 to mm^2
  mutate(dry = dry*1000) # g to mg

mean_traits2 %<>%
  mutate(SLA = area/dry) %>% # mm^2 mg-1
  mutate(LMA = dry/area) %>%
  mutate(LDMC = dry/wet) # mg g-1

mean_traits2 %<>%
  dplyr::select(-c(wet, length, dry))

#' join stomata and CN trait values to leaf trait data
mean_traits2 <- left_join(mean_traits2, stomata_joined)
mean_traits2 <- left_join(mean_traits2, cn_data, by = "ID")

#' join with metadata
mean_traits2 <- left_join(mean_traits2, all_tree_meta_data) 

#' remove sample labels if didn't correspond to metadata
mean_traits2 <- mean_traits2[!is.na(mean_traits2$Garden), ]

#'  Output of all trait data trees + meta data
write.csv(mean_traits2, "data/gardenData/trait_plus_meta_csvs/all_trait_wmeta.csv", row.names = FALSE)

mean_traits2 = mutate(mean_traits2,
                      ID_Garden = paste(ID, Garden, sep = "_"))

#' Output of trait data only corresponding with microbiome data (all blocks) 
meta = mutate(meta,
              ID_Garden = paste("C", C, "R", R, "B", Block, "_", Garden, sep = ""))

all_traits <- mean_traits2 %>% ungroup() %>% dplyr::select(thick:CtoN, ID_Garden)
trait_wmeta_allblocks <- left_join(meta, all_traits, by = "ID_Garden") 

write.csv(trait_wmeta_allblocks, "data/gardenData/trait_plus_meta_csvs/trait_wmeta_allblocks.csv", row.names = FALSE)

#' Output of trait data only corresponding with microbiome data (only trait blocks) 
trait_meta = mutate(trait_meta,
                    ID_Garden = paste("C", C, "R", R, "B", Block, "_", Garden, sep = ""))

trait_wmeta_traitblocks <- left_join(trait_meta, all_traits, by = "ID_Garden") 

write.csv(trait_wmeta_traitblocks, "data/gardenData/trait_plus_meta_csvs/trait_wmeta_traitblocks.csv", row.names = FALSE)


#' Genotype averages of leaf traits calculated
#' NAs were ignored - if NA's were not ignored there would be only 40/118 stomata values)
genotype_mean_trait <- mean_traits2 %>%
  group_by(Garden, Genotype) %>%
  mutate(across(area:percent_C_weight, mean, na.rm = TRUE)) %>%
  dplyr::select(Garden, Region, Genotype, thick, SLA, LDMC, CtoN, stomata_length, area, LMA, percent_C_weight,percent_N_weight, stomata_count) %>%
  unique()

#' Almost all genotypes are replicated 4 times
#' However, the number of leaf traits associated with those replicates may vary
genotype_count <- mean_traits2 %>%
  group_by(Garden, Genotype) %>%
  count()

#' 21/118 stomata have 0 or 1 rep/genotype/garden
test <- mean_traits2 %>% 
  group_by(Garden, Genotype) %>% 
  summarise(is_na = sum(is.na(stomata_count)))

stomata_test <- left_join(genotype_count, test)

stomata_test <- stomata_test %>%
  mutate(reps = n - is_na) 

sum(stomata_test$reps == 0 | stomata_test$reps == 1)

#' Reshape data
genotype_mean_trait <- genotype_mean_trait[!is.na(genotype_mean_trait$Garden), ]
trait_join <- genotype_mean_trait %>% reshape2::melt() %>% drop_na()

#### CODE FOR TRAIT ANOVA FIG ####

#' Use the G,E,GxE values from output to input into trait_anova_fig semi_r2 values

mean_traits2 = mutate(mean_traits2,
                      Block_Garden = paste(Block, "_", Garden, sep = ""))

mean_traits2 = mutate(mean_traits2,
                      Region_Garden = paste(Region, "_", Garden, sep = ""))

thick_df <- subset(mean_traits2, !is.na(thick))
LDMC_df <- subset(mean_traits2, !is.na(LDMC))
SLA_df <- subset(mean_traits2, !is.na(SLA))
count_df <- subset(mean_traits2, !is.na(stomata_count))
length_df <- subset(mean_traits2, !is.na(stomata_length))
cn_df <- subset(mean_traits2, !is.na(CtoN))

#thick
modGP1_thick <- lmer(thick ~  Garden*Genotype + (1|Block_Garden), data = thick_df)
R2_GPc_part1_thick <- partR2(modGP1_thick,  partvars = c("Garden:Genotype"), data = thick_df, nboot = 10)

modGP2_thick <- lmer(thick ~ Garden + Genotype + (1|Block_Garden), data = thick_df)
R2_GPc_part2_thick <- partR2(modGP2_thick, partvars = c("Garden", "Genotype"), data = thick_df,
                             max_level = 1, nboot = 10)

R2_GPc_thick <- mergeR2(R2_GPc_part1_thick, R2_GPc_part2_thick) 
R2_GPc_thick

#SLA
modGP1_SLA <- lmer(SLA ~  Garden*Genotype + (1|Block_Garden), data = SLA_df)
R2_GPc_part1_SLA <- partR2(modGP1_SLA,  partvars = c("Garden:Genotype"), data = SLA_df, nboot = 10)

modGP2_SLA <- lmer(SLA ~ Garden + Genotype + (1|Block_Garden), data = SLA_df)
R2_GPc_part2_SLA <- partR2(modGP2_SLA, partvars = c("Garden", "Genotype"), data = SLA_df,
                           max_level = 1, nboot = 10)

R2_GPc_SLA <- mergeR2(R2_GPc_part1_SLA, R2_GPc_part2_SLA) 
R2_GPc_SLA

#LDMC
modGP1_LDMC <- lmer(LDMC ~  Garden*Genotype + (1|Block_Garden), data = LDMC_df)
R2_GPc_part1_LDMC <- partR2(modGP1_LDMC,  partvars = c("Garden:Genotype"), data = LDMC_df, nboot = 10)

modGP2_LDMC <- lmer(LDMC ~ Garden + Genotype + (1|Block_Garden), data = LDMC_df)
R2_GPc_part2_LDMC <- partR2(modGP2_LDMC, partvars = c("Garden", "Genotype"), data = LDMC_df,
                            max_level = 1, nboot = 10)

R2_GPc_LDMC <- mergeR2(R2_GPc_part1_LDMC, R2_GPc_part2_LDMC) 
R2_GPc_LDMC

#stomata density
modGP1_count <- lmer(stomata_count ~  Garden*Genotype + (1|Block_Garden), data = count_df)
R2_GPc_part1_count <- partR2(modGP1_count,  partvars = c("Garden:Genotype"), data = count_df, nboot = 10)

modGP2_count <- lmer(stomata_count ~ Garden + Genotype + (1|Block_Garden), data = count_df)
R2_GPc_part2_count <- partR2(modGP2_count, partvars = c("Garden", "Genotype"), data = count_df,
                             max_level = 1, nboot = 10)

R2_GPc_count <- mergeR2(R2_GPc_part1_count, R2_GPc_part2_count) 
R2_GPc_count

#stomata length
modGP1_length <- lmer(stomata_length ~  Garden*Genotype + (1|Block_Garden), data = length_df)
R2_GPc_part1_length <- partR2(modGP1_length,  partvars = c("Garden:Genotype"), data = length_df, nboot = 10)

modGP2_length <- lmer(stomata_length ~ Garden + Genotype + (1|Block_Garden), data = length_df)
R2_GPc_part2_length <- partR2(modGP2_length, partvars = c("Garden", "Genotype"), data = length_df,
                              max_level = 1, nboot = 10)

R2_GPc_length <- mergeR2(R2_GPc_part1_length, R2_GPc_part2_length) 
R2_GPc_length

#CtoN - Can't include because samples were compiled by genotype


##### Genotype by environment interaction ANOVA FOR THESIS 6/16   ######

#' Make trait dataframe into list
traits_for_anova <- mean_traits2 %>% dplyr::select(thick:stomata_count, Garden, Genotype)

traits_for_anova %<>% pivot_longer(cols = thick:stomata_count,
                                   names_to = "trait",
                                   values_to = "value")

traits_for_anova %<>% remove_missing()

traits_for_anova <- split(traits_for_anova, f = traits_for_anova$trait)

anova <- function(df){
  fit1 <- aov(value ~ Genotype*Garden, data = df)
  summary(fit1) 
  fit1df <- tidy(fit1)
}

anova_list <- traits_for_anova %>% lapply(anova)

anova_df <- reshape2::melt(anova_list)

anova_df <- pivot_wider(anova_df, names_from = variable)

write_csv(anova_df, "output/csvs/leaf_trait_anovas.csv")

# Genotype effect for leaf traits within garden/region
wc_thick_df <- thick_df %>% filter(Region_Garden == "West_Corvallis")
ec_thick_df <- thick_df %>% filter(Region_Garden == "East_Corvallis")
wh_thick_df <- thick_df %>% filter(Region_Garden == "West_Hermiston")
eh_thick_df <- thick_df %>% filter(Region_Garden == "East_Hermiston")

wc_SLA_df <- SLA_df %>% filter(Region_Garden == "West_Corvallis")
ec_SLA_df <- SLA_df %>% filter(Region_Garden == "East_Corvallis")
wh_SLA_df <- SLA_df %>% filter(Region_Garden == "West_Hermiston")
eh_SLA_df <- SLA_df %>% filter(Region_Garden == "East_Hermiston")

wc_LDMC_df <- LDMC_df %>% filter(Region_Garden == "West_Corvallis")
ec_LDMC_df <- LDMC_df %>% filter(Region_Garden == "East_Corvallis")
wh_LDMC_df <- LDMC_df %>% filter(Region_Garden == "West_Hermiston")
eh_LDMC_df <- LDMC_df %>% filter(Region_Garden == "East_Hermiston")

wc_count_df <- count_df %>% filter(Region_Garden == "West_Corvallis")
ec_count_df <- count_df %>% filter(Region_Garden == "East_Corvallis")
wh_count_df <- count_df %>% filter(Region_Garden == "West_Hermiston")
eh_count_df <- count_df %>% filter(Region_Garden == "East_Hermiston")

wc_length_df <- length_df %>% filter(Region_Garden == "West_Corvallis")
ec_length_df <- length_df %>% filter(Region_Garden == "East_Corvallis")
wh_length_df <- length_df %>% filter(Region_Garden == "West_Hermiston")
eh_length_df <- length_df %>% filter(Region_Garden == "East_Hermiston")

#find p-values
anova_geno <- function(trait_name, df){
  trait_name.aov <- lm(trait_name ~ Genotype, data = df) 
  aov <- summary(aov(trait_name.aov))
  aov[[1]][["Pr(>F)"]][1]
}

anova_geno(wc_thick_df$thick, wc_thick_df)
anova_geno(ec_thick_df$thick, ec_thick_df)
anova_geno(wh_thick_df$thick, wh_thick_df) #ns
anova_geno(eh_thick_df$thick, eh_thick_df)

anova_geno(wc_SLA_df$SLA, wc_SLA_df)
anova_geno(ec_SLA_df$SLA, ec_SLA_df) #ns
anova_geno(wh_SLA_df$SLA, wh_SLA_df) 
anova_geno(eh_SLA_df$SLA, eh_SLA_df)

anova_geno(wc_LDMC_df$LDMC, wc_LDMC_df)
anova_geno(ec_LDMC_df$LDMC, ec_LDMC_df)
anova_geno(wh_LDMC_df$LDMC, wh_LDMC_df) #ns
anova_geno(eh_LDMC_df$LDMC, eh_LDMC_df)

anova_geno(wc_length_df$stomata_length, wc_length_df)
anova_geno(ec_length_df$stomata_length, ec_length_df)
anova_geno(wh_length_df$stomata_length, wh_length_df)
anova_geno(eh_length_df$stomata_length, eh_length_df)

anova_geno(wc_count_df$stomata_count, wc_count_df) #ns
anova_geno(ec_count_df$stomata_count, ec_count_df)
anova_geno(wh_count_df$stomata_count, wh_count_df) #ns
anova_geno(eh_count_df$stomata_count, eh_count_df)

#find R2
anova_geno <- function(trait_name, df){
  trait_name.aov <- lm(trait_name ~ Genotype, data = df) 
  anova(trait_name.aov)
  summary(trait_name.aov)$r.squared
}

wcthick <- anova_geno(wc_thick_df$thick, wc_thick_df)
ecthick <- anova_geno(ec_thick_df$thick, ec_thick_df)
whthick <- anova_geno(wh_thick_df$thick, wh_thick_df) #ns
ehthick <- anova_geno(eh_thick_df$thick, eh_thick_df)

wcSLA <- anova_geno(wc_SLA_df$SLA, wc_SLA_df)
ecSLA <- anova_geno(ec_SLA_df$SLA, ec_SLA_df) 
whSLA <- anova_geno(wh_SLA_df$SLA, wh_SLA_df) #ns
ehSLA <- anova_geno(eh_SLA_df$SLA, eh_SLA_df)

wcLDMC <- anova_geno(wc_LDMC_df$LDMC, wc_LDMC_df)
ecLDMC <- anova_geno(ec_LDMC_df$LDMC, ec_LDMC_df)
whLDMC <- anova_geno(wh_LDMC_df$LDMC, wh_LDMC_df) #ns
ehLDMC <- anova_geno(eh_LDMC_df$LDMC, eh_LDMC_df)

wclength <- anova_geno(wc_length_df$stomata_length, wc_length_df)
eclength <- anova_geno(ec_length_df$stomata_length, ec_length_df)
whlength <- anova_geno(wh_length_df$stomata_length, wh_length_df)
ehlength <- anova_geno(eh_length_df$stomata_length, eh_length_df)

wccount <- anova_geno(wc_count_df$stomata_count, wc_count_df) #ns
eccount <- anova_geno(ec_count_df$stomata_count, ec_count_df)
whcount <- anova_geno(wh_count_df$stomata_count, wh_count_df) #ns
ehcount <- anova_geno(eh_count_df$stomata_count, eh_count_df)

geno_effect <- c(wcthick, wcSLA, wcLDMC, wccount, wclength, whthick, whSLA, whLDMC, whcount, whlength, ecthick, ecSLA, ecLDMC, eccount, eclength, ehthick, ehSLA, ehLDMC, ehcount, ehlength)
garden <- c("Corvallis", "Corvallis", "Corvallis", "Corvallis", "Corvallis", "Hermiston", "Hermiston", "Hermiston", "Hermiston", "Hermiston", "Corvallis", "Corvallis", "Corvallis", "Corvallis", "Corvallis", "Hermiston", "Hermiston", "Hermiston", "Hermiston", "Hermiston")
measurement <- c("thick", "SLA", "LDMC", "count", "length", "thick", "SLA", "LDMC", "count", "length", "thick", "SLA", "LDMC", "count", "length", "thick", "SLA", "LDMC", "count", "length")
ecotype <- c("West", "West", "West", "West", "West", "West", "West", "West", "West", "West", "East", "East", "East", "East", "East", "East", "East", "East", "East", "East")
significance <- c(1, 1, 1, 0, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)

#this data will be used in microbiome_trait_geno_effect.R
ew_geno_df <- data.frame(geno_effect, garden, measurement, ecotype, significance)
write.csv(ew_geno_df, "output/csvs/ew_geno_effect.csv")
