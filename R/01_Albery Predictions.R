
# 01_Albery Predictions ####

rm(list = ls())

# Prediction ####

library(ggrepel); library(tidyverse); library(SpRanger); library(cowplot); library(patchwork)
library(ggregplot); library(data.table)

theme_set(theme_cowplot())

Panth1 <- read.delim("Data/PanTHERIA_1-0_WR05_Aug2008.txt") %>%
  dplyr::rename_all(~str_replace(.x, "MSW05_", "h")) %>%
  rename(Sp = hBinomial)

Panth1$Sp <- Panth1$Sp %>% str_replace(" ", "_")

# load("~/LargeFiles/MammalStackFullMercator.Rdata")

PredictedNetwork <- fread("Data/AlberyPredicted.csv") %>% mutate_all(as.numeric)

PredictedNetwork %>% as.matrix -> PredictedNetwork

PredictedNetwork <- 
  PredictedNetwork[,-1]

rownames(PredictedNetwork) <- colnames(PredictedNetwork)

read.csv("Data/Virionette/03_interaction_data/virionette.csv") ->
  
  BetaCovHosts

# Bat-Specific ####

BetaCovHosts %>% 
  filter(virus_genus == "Betacoronavirus") %>% 
  filter(host_order == "Chiroptera") %>% 
  pull(host_species) %>% as.character %>% 
  str_trim %>%
  str_replace_all(" ", "_") %>%
  intersect(rownames(PredictedNetwork)) %>% sort -> 
  BetaCovBats

NetworkPredict(BetaCovBats, (PredictedNetwork), IncludeObserved = T) %>%
  as.data.frame() %>% 
  left_join(Panth1, by = "Sp") %>%
  filter(hOrder == "Chiroptera") %>% 
  dplyr::select(1:3, Sp, Observed, hOrder, hFamily, hGenus) ->
  
  BetaCovPredictedBats

BetaCovPredictedBats %>% 
  select(Sp, Count) %>%
  write.csv("AlberyBats.csv")

# Non-Bats ####

BetaCovHosts %>% 
  filter(virus_genus == "Betacoronavirus") %>% 
  #filter(host_order == "Chiroptera") %>% 
  pull(host_species) %>% as.character %>% 
  str_trim %>%
  str_replace_all(" ", "_") %>%
  intersect(rownames(PredictedNetwork)) %>% sort -> 
  BetaCovMammals

NetworkPredict(BetaCovMammals, (PredictedNetwork), IncludeObserved = T) %>%
  as.data.frame() %>% 
  left_join(Panth1, by = "Sp") %>%
  #filter(!hOrder == "Chiroptera") %>% 
  dplyr::select(1:3, Sp, Observed, hOrder, hFamily, hGenus) ->
  
  BetaCovPredictedMammals

BetaCovPredictedMammals %>% 
  select(Sp, Count) %>%
  write.csv("AlberyNonBats.csv")


# R_affinis ####

NetworkPredict(c("Rhinolophus_affinis"), (PredictedNetwork)) %>%
  as.data.frame() %>% left_join(Panth1, by = "Sp") %>%
  #filter(!hOrder == "Chiroptera") %>% 
  dplyr::select(1:3, Sp, hOrder, hFamily, hGenus) ->
  R_affinisPredictedBats

R_affinisPredictedBats %>% nrow

NetworkPredict(c("Rhinolophus_affinis"), as.matrix(PredictedNetwork)) %>%
  as.data.frame() %>% left_join(Panth1, by = "Sp") %>%
  filter(!hOrder == "Chiroptera") %>%
  dplyr::select(1:3, Sp, hOrder, hFamily, hGenus) ->
  
  R_affinisPredictedNonBats

R_affinisPredictedNonBats %>% select(Sp, Count) %>%
  saveRDS(file = "Intermediate/AlberyPredictedNonBats_R_affinis.rds")

R_affinisPredictedBats %>% 
  saveRDS(file = "Intermediate/AlberyPredictedBats_R_affinis.rds")

# R_malayanus ####

NetworkPredict(c("Rhinolophus_malayanus"), (PredictedNetwork)) %>%
  as.data.frame() %>% left_join(Panth1, by = "Sp") %>%
  #filter(!hOrder == "Chiroptera") %>% 
  dplyr::select(1:3, Sp, hOrder, hFamily, hGenus) ->
  R_malayanusPredictedBats

NetworkPredict(c("Rhinolophus_malayanus"), as.matrix(PredictedNetwork)) %>%
  as.data.frame() %>% left_join(Panth1, by = "Sp") %>%
  filter(!hOrder == "Chiroptera") %>%
  dplyr::select(1:3, Sp, hOrder, hFamily, hGenus) ->
  
  R_malayanusPredictedNonBats

R_malayanusPredictedNonBats %>% select(Sp, Count) %>%
  saveRDS(file = "Intermediate/AlberyPredictedNonBats_R_malayanus.rds")

R_malayanusPredictedBats %>% 
  saveRDS(file = "Intermediate/AlberyPredictedBats_R_malayanus.rds")

R_affinisPredictedBats %>% 
  mutate_at("Sp", ~str_replace_all(.x, "_", " ")) %>%
  mutate(Sp = glue::glue("{1:n()}. {Sp} (P={Count})")) %>%
  slice(1:20) %>% select(Sp) %>%
  bind_cols(R_malayanusPredictedBats %>% 
              mutate_at("Sp", ~str_replace_all(.x, "_", " ")) %>%
              mutate(Sp = glue::glue("{1:n()}. {Sp} (P={Count})")) %>%
              slice(1:20) %>% select(Sp)) %>%
  rename(`R.affinis` = Sp, `R.malayanus` = Sp1) %>% 
  write.csv("Output Files/AlberyRhinolophusBatPredictions.csv", row.names = F)

R_affinisPredictedNonBats %>% 
  mutate_at("Sp", ~str_replace_all(.x, "_", " ")) %>%
  mutate(Sp = glue::glue("{1:n()}. {Sp} (P={Count})")) %>%
  slice(1:20) %>% select(Sp) %>%
  bind_cols(R_malayanusPredictedNonBats %>% 
              mutate_at("Sp", ~str_replace_all(.x, "_", " ")) %>%
              mutate(Sp = glue::glue("{1:n()}. {Sp} (P={Count})")) %>%
              slice(1:20) %>% select(Sp)) %>%
  rename(`R.affinis` = Sp, `R.malayanus` = Sp1) %>% 
  write.csv("Output Files/AlberyRhinolophusNonBatPredictions.csv", row.names = F)
