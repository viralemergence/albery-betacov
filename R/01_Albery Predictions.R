
# 01_Albery Predictions ####

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

read.delim("Data/BetaCov.txt", sep = ",") ->
  
  BetaCovHosts

BetaCovHosts %>% 
  filter(virus_genus == "Betacoronavirus") %>% 
  pull(clean_hostnames) %>% as.character %>% 
  intersect(rownames(PredictedNetwork)) %>% sort -> 
  BetaCovHosts2

NetworkPredict(BetaCovHosts2, (PredictedNetwork), IncludeObserved = T) %>%
  as.data.frame() %>% 
  left_join(Panth1, by = "Sp") %>%
  filter(hOrder == "Chiroptera") %>% 
  dplyr::select(1:3, Sp, Observed, hOrder, hFamily, hGenus) ->
  
  BetaCovPredictedBats

BetaCovPredictedBats %>% 
  select(Sp, Count) %>%
  write.csv("AlberyPredicted.csv")

# BetaCovPredictedBats %>% write.csv("Github/CSVs/AlberyPredicted.csv", row.names = F)

NetworkPredict(BetaCovHosts2, (PredictedNetwork), IncludeObserved = T) %>%
  as.data.frame() %>% left_join(Panth1, by = "Sp") %>%
  filter(!hOrder == "Chiroptera") %>% 
  dplyr::select(1:3, Sp, Observed, hOrder, hFamily, hGenus) ->
  
  BetaCovPredictedNonBats

NetworkValidate(BetaCovHosts2, (PredictedNetwork)) %>%
  as.data.frame() %>% left_join(Panth1, by = "Sp") %>%
  #filter(!hOrder == "Chiroptera") %>% 
  dplyr::select(1:3, Sp, hOrder, hFamily, hGenus) ->
  BetacovPredictedBats

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
