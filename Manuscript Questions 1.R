library(tidyverse)
library(tidyr)
library(compositions)
library(ggplot2)
library(dplyr)
library(vegan)
# load data

## A. Comp and div of settled microbes (and sub)
## 1. Tile/bar/point plots ##
#### Presence/absence heat maps ####
# 18S
mc_tmp_18 <- asv_wtax_18 |> 
  select(-SAMPLE) |> 
  filter(SAMPLETYPE == "Microcolonizer" | VENT == "Mt Edwards" | VENT == "Deep seawater") |> 
  filter(Domain == "Eukaryota") |> 
  mutate(SUBSTRATE = str_replace_all(Substrate, "z2", "z")) |> 
  mutate(SAMPLE = case_when(
    SAMPLETYPE == "Microcolonizer" ~ paste(MC, SUBSTRATE, sep = "-"),
    VENT == "Deep seawater" ~ "Background seawater",
    VENT == "Mt Edwards" ~ "Mt Edwards diffuse fluid"
  )) |>
  group_by(SAMPLE, FeatureID, Taxon, Domain, 
           Supergroup, Phylum, Class, Order, Family, Genus, Species) |> 
  summarize(SEQ_AVG = mean(SEQUENCE_COUNT)) |> 
  filter(SEQ_AVG > 0)

format_18s <- mc_tmp_18 %>%
  mutate(SAMPLE_TYPE = case_when(
    SAMPLE == "Background seawater" ~ "Background", 
    SAMPLE == "Mt Edwards diffuse fluid" ~ "Vent", 
    TRUE ~ "Substrate")) %>%
  select(FeatureID, SEQ_AVG, SAMPLE_TYPE) %>%
  pivot_wider(id_cols = FeatureID, names_from = SAMPLE_TYPE, values_fn = mean, values_from = SEQ_AVG, values_fill = NA)

view(format_18s)

#assign placement of ASVs
format_18s_assigned <- format_18s %>%
  mutate(DISTRIBUTION = case_when(
    (is.na(Background) & is.na(Substrate) & !(is.na(Vent))) ~ "Vent only", 
    (is.na(Background) & !(is.na(Substrate)) & !(is.na(Vent))) ~ "Vent & Substrate only", 
    (is.na(Background) & !(is.na(Substrate)) & is.na(Vent)) ~ "Substrate only", 
    (!(is.na(Background)) & is.na(Substrate) & !(is.na(Vent))) ~ "Vent & Background only",
    (!(is.na(Background)) & !(is.na(Substrate)) & is.na(Vent)) ~ "Substrate & Background only", 
    (!(is.na(Background)) & is.na(Substrate) & is.na(Vent)) ~ "Background only", 
    (!(is.na(Background)) & !(is.na(Substrate)) & !(is.na(Vent)) ~ "All")
  ))
view(format_18s_assigned)

key_dist_18 <- format_18s_assigned %>%
  select(FeatureID, DISTRIBUTION) %>%
  distinct()

df_18s_all <- mc_tmp_18 %>%
  left_join(key_dist_18)

mc_df_18s <- df_18s_all %>%
  filter(SAMPLE != "Background seawater" & SAMPLE != "Mt Edwards diffuse fluid") %>%
  filter(Phylum != "Metazoa" & Phylum != "Streptophyta") %>%
  filter(!is.na(Phylum)) %>%
  separate(SAMPLE, into = c("MC", "SUBSTRATE"), sep = "-", remove = FALSE) %>%
  add_column(PRESENCE = 1) 

ggplot(mc_df_18s, aes(x = MC, y = FeatureID, fill = PRESENCE)) +
  geom_tile(stat = "identity") + 
  theme_classic() +
  facet_grid(cols = vars(SUBSTRATE), rows = vars(Supergroup), scales = "free", space = "free") +
  theme(axis.text.y = element_blank(), axis.text.x = element_text(angle = 90), strip.text.y = element_text(angle = 0))

ggplot(mc_df_18s, aes(x = MC, y = FeatureID, fill = PRESENCE)) +
  geom_tile(stat = "identity") + 
  theme_classic() +
  facet_grid(cols = vars(SUBSTRATE), rows = vars(Phylum), scales = "free", space = "free") +
  theme(axis.text.y = element_blank(), axis.text.x = element_text(angle = 90), strip.text.y = element_text(angle = 0))

# 16s  
mc_tmp_16 <- asv16s_df |> 
  filter(Sampletype == "Microcolonizer" | LocationName == "Mt Edwards Vent" | LocationName == "Deep seawater" & STATUS == "keep") |> 
  mutate(SAMPLE = case_when(
    Sampletype == "Microcolonizer" ~ paste(MC, Substrate, sep = "-"),
    LocationName == "Deep seawater" ~ "Background seawater",
    LocationName == "Mt Edwards Vent" ~ "Mt Edwards diffuse fluid"
  )) |> 
  mutate(SAMPLE  = case_when(
    `SAMPLE` == "MC3-shell" ~ "MC3-Shell",
    .default = `SAMPLE`)) |>
  group_by(SAMPLE, FeatureID, Taxon, Domain, 
           Phylum, Class, Order, Family, Genus, Species) |> 
  summarize(SEQ_AVG = mean(SEQUENCE_COUNT)) |> 
  filter(!(Domain == "Unassigned") & !(Domain == "Eukaryota")) |> 
  filter(SEQ_AVG > 0)
format_16s <- mc_tmp_16 %>%
  filter(!(grepl("Olivine", SAMPLE)) & 
           !(grepl("Basalt", SAMPLE)) &
           !(grepl("Pyrite", SAMPLE))) %>%
  mutate(SAMPLE_TYPE = case_when(
    SAMPLE == "Background seawater" ~ "Background", 
    SAMPLE == "Mt Edwards diffuse fluid" ~ "Vent", 
    TRUE ~ "Substrate"
  )) %>%
  select(FeatureID, SEQ_AVG, SAMPLE_TYPE) %>%
  pivot_wider(id_cols = FeatureID, names_from = SAMPLE_TYPE, values_fn = mean, values_from = SEQ_AVG, values_fill = NA)

#assign placement of ASVs
format_16s_assigned <- format_16s %>%
  mutate(DISTRIBUTION = case_when(
    (is.na(Background) & is.na(Substrate) & !(is.na(Vent))) ~ "Vent only", 
    (is.na(Background) & !(is.na(Substrate)) & !(is.na(Vent))) ~ "Vent & Substrate only", 
    (is.na(Background) & !(is.na(Substrate)) & is.na(Vent)) ~ "Substrate only", 
    (!(is.na(Background)) & is.na(Substrate) & !(is.na(Vent))) ~ "Vent & Background only",
    (!(is.na(Background)) & !(is.na(Substrate)) & is.na(Vent)) ~ "Substrate & Background only", 
    (!(is.na(Background)) & is.na(Substrate) & is.na(Vent)) ~ "Background only", 
    (!(is.na(Background)) & !(is.na(Substrate)) & !(is.na(Vent)) ~ "All")
  ))

key_dist <- format_16s_assigned %>%
  select(FeatureID, DISTRIBUTION) %>%
  distinct()

df_16s_all <- mc_tmp_16 %>%
  filter(!(grepl("Olivine", SAMPLE)) & 
           !(grepl("Basalt", SAMPLE)) &
           !(grepl("Pyrite", SAMPLE))) %>%
  left_join(key_dist)

# on substrates
mc_df_16s <- df_16s_all %>%
  filter(SAMPLE != "Background seawater" & SAMPLE != "Mt Edwards diffuse fluid") %>%
  filter(!is.na(Phylum)) %>%
  separate(SAMPLE, into = c("MC", "SUBSTRATE"), sep = "-") %>%
  add_column(PRESENCE = 1) 
unique(mc_df_18s$Phylum)

ggplot(mc_df_16s, aes(x = MC, y = FeatureID, fill = PRESENCE)) +
  geom_tile(stat = "identity") + 
  theme_classic() +
  facet_grid(cols = vars(SUBSTRATE), rows = vars(Phylum), scales = "free", space = "free") +
  theme(axis.text.y = element_blank(), axis.text.x = element_text(angle = 45), strip.text.y = element_text(angle = 0))

ggplot(mc_df_16s, aes(x = SUBSTRATE, y = FeatureID, fill = PRESENCE)) +
  geom_tile(stat = "identity") + 
  theme_classic() +
  facet_grid(cols = vars(MC), rows = vars(Phylum), scales = "free", space = "free") +
  theme(axis.text.y = element_blank(), axis.text.x = element_text(angle = 45), strip.text.y = element_text(angle = 0))

####
#### Bar plots of Supergroup to show Phylum ####
subs_18s <- df_18s_all %>%
  filter(DISTRIBUTION == "Substrate only") %>%
  filter(!is.na(Phylum)) %>%
  filter(Phylum != "Metazoa" & Phylum != "Streptophyta") %>%
  add_column(COUNT = 1) %>%
  separate(SAMPLE, c("Microcolonizer", "Substrate"), sep = "-")
unique(subs_18s$Phylum)

rel_abund_18s <- subs_18s %>%
  ungroup() %>%
  group_by(Supergroup, Phylum) %>% #Microcolonizer, Substrate, 
  summarise(SUM_PHYLUM = sum(COUNT))

unique(rel_abund_18s$Supergroup)
# "Alveolata" "Amoebozoa" "Archaeplastida" "Hacrobia" "Opisthokonta" "Rhizaria" "Stramenopiles"  "Excavata"  "Apusozoa"

rel_abund_alv <- rel_abund_18s %>%
  filter(Supergroup == "Alveolata")

rel_abund_Amoeb <- rel_abund_18s %>%
  filter(Supergroup == "Amoebozoa")

rel_abund_Apus <- rel_abund_18s %>%
  filter(Supergroup == "Apusozoa")

rel_abund_Arch <- rel_abund_18s %>%
  filter(Supergroup == "Archaeplastida")

rel_abund_Ex <- rel_abund_18s %>%
  filter(Supergroup == "Excavata")

rel_abund_Hac <- rel_abund_18s %>%
  filter(Supergroup == "Hacrobia")

rel_abund_Opistho <- rel_abund_18s %>%
  filter(Supergroup == "Opisthokonta")

rel_abund_Rhiz <- rel_abund_18s %>%
  filter(Supergroup == "Rhizaria")

rel_abund_Stram <- rel_abund_18s %>%
  filter(Supergroup == "Stramenopiles")

ggplot(rel_abund_alv, aes(x = Supergroup, SUM_PHYLUM, fill = Phylum)) +
  geom_bar(stat = "identity", position = "fill", color = "black") +
  labs(x = "", y = "") +
  theme_minimal() +
  coord_flip()

ggplot(rel_abund_Amoeb, aes(x = Supergroup, SUM_PHYLUM, fill = Phylum)) +
  geom_bar(stat = "identity", position = "fill", color = "black") +
  labs(x = "", y = "") +
  theme_minimal() +
  coord_flip()

ggplot(rel_abund_Apus, aes(x = Supergroup, SUM_PHYLUM, fill = Phylum)) +
  geom_bar(stat = "identity", position = "fill", color = "black") +
  labs(x = "", y = "") +
  theme_minimal() +
  coord_flip()

ggplot(rel_abund_Arch, aes(x = Supergroup, SUM_PHYLUM, fill = Phylum)) +
  geom_bar(stat = "identity", position = "fill", color = "black") +
  labs(x = "", y = "") +
  theme_minimal() +
  coord_flip()

ggplot(rel_abund_Ex, aes(x = Supergroup, SUM_PHYLUM, fill = Phylum)) +
  geom_bar(stat = "identity", position = "fill", color = "black") +
  labs(x = "", y = "") +
  theme_minimal() +
  coord_flip()

ggplot(rel_abund_Hac, aes(x = Supergroup, SUM_PHYLUM, fill = Phylum)) +
  geom_bar(stat = "identity", position = "fill", color = "black") +
  labs(x = "", y = "") +
  theme_minimal() +
  coord_flip()

ggplot(rel_abund_Opistho, aes(x = Supergroup, SUM_PHYLUM, fill = Phylum)) +
  geom_bar(stat = "identity", position = "fill", color = "black") +
  labs(x = "", y = "") +
  theme_minimal() +
  coord_flip()

ggplot(rel_abund_Rhiz, aes(x = Supergroup, SUM_PHYLUM, fill = Phylum)) +
  geom_bar(stat = "identity", position = "fill", color = "black") +
  labs(x = "", y = "") +
  theme_minimal() +
  coord_flip()

ggplot(rel_abund_Stram, aes(x = Supergroup, SUM_PHYLUM, fill = Phylum)) +
  geom_bar(stat = "identity", position = "fill", color = "black") +
  labs(x = "", y = "") +
  theme_minimal() +
  coord_flip()
#### Relative ASV abundance across sampled sites - stacked bar (w and w/out Streptophyta and Metazoa) ####
#16S 
stacked_asvs_16S <- df_16s_all %>%
  filter(SEQ_AVG > 10) %>%
  group_by(SAMPLE, Phylum) %>%
  summarise(Count = n(), .groups = 'drop') %>%
  select(SAMPLE, Phylum, Count) %>%
  separate(SAMPLE, into = c("MC", "SUBSTRATE"), sep = "-", remove = FALSE, fill = "right")

#  filter(!(is.na(Phylum)))
stacked_asvs_16S$ordered_SAMPLES_16S <- factor(stacked_asvs_16S$SAMPLE, levels = c("MC2-Quartz", "MC2-Riftia", "MC3-Quartz","MC3-Riftia", "MC3-Shell", "MC1-Quartz", "MC1-Riftia", "MC1-Shell", "MC4-Quartz", "MC4-Riftia", "MC5-Quartz", "MC5-Riftia", "MC5-Shell", "MC6-Quartz", "MC6-Riftia", "MC6-Shell", "Mt Edwards diffuse fluid", "Background seawater"))

ggplot(stacked_asvs_16S, aes(x = ordered_SAMPLES_16S, y = Count, fill = Phylum)) +
  geom_bar(stat = "identity", position = "fill", color = "black") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = -45, hjust = 0)) + 
  labs(x = "", y = "Relative ASV Abundance")
ggplot(stacked_asvs_16S, aes(x = ordered_SAMPLES_16S, y = Count, fill = Phylum)) +
  geom_bar(stat = "identity", position = "fill", color = "black") +
  theme_classic() +
  facet_grid(cols = vars(MC), scales = "free", space = "free") +
  theme(axis.text.x = element_text(angle = -45, hjust = 0)) + 
  labs(x = "", y = "Relative ASV Abundance")
ggplot(stacked_asvs_16S, aes(x = ordered_SAMPLES_16S, y = Count, fill = Phylum)) +
  geom_bar(stat = "identity", position = "fill", color = "black") +
  theme_classic() +
  facet_grid(cols = vars(SUBSTRATE), scales = "free", space = "free") +
  theme(axis.text.x = element_text(angle = -45, hjust = 0)) + 
  labs(x = "", y = "Relative ASV Abundance")

#18S
stacked_asvs <- mc_18S_df %>%
  group_by(SAMPLE, Phylum, Supergroup)%>%
  summarise(Count = n(), .groups = 'drop') %>%
  select(SAMPLE, Supergroup, Phylum, Count) %>%
  unite(TAXA_LEVEL, Supergroup, Phylum, sep = "-", remove = FALSE) %>%
  filter(!(is.na(Phylum)))
stacked_asvs$ordered_SAMPLES_18S <- factor(stacked_asvs$SAMPLE, levels = c("MC2-Quartz", "MC2-Riftia", "MC3-Riftia", "MC3-Shell", "MC1-Quartz", "MC1-Riftia", "MC1-Shell", "MC4-Quartz", "MC4-Shell", "MC5-Quartz", "MC5-Riftia", "MC5-Shell", "MC6-Quartz", "MC6-Riftia", "MC6-Shell", "Mt Edwards diffuse fluid", "Background seawater"))

ggplot(stacked_asvs, aes(x = ordered_SAMPLES_18S, y = Count, fill = TAXA_LEVEL)) +
  geom_bar(stat = "identity", position = "fill", color = "black") +
  theme_classic() +
  scale_fill_manual(values = c("#f1eef6", "#d7b5d8", "#df65b0", "#dd1c77", 
                               "#fde0dd", "#fa9fb5", "#c51b8a", "#edf8fb", "#bfd3e6", "#9ebcda", 
                               "#8c96c6", "#8856a7", "#000000", "#810f7c","#ffffcc", "#f0f9e8", "#c7e9b4", "#7fcdbb", 
                               "#41b6c4", "#2c7fb8", "#253494",  "#fee391",  "#fec44f", "#fe9929", "#000000","#000000","#d95f0e", "#993404", "#fc9272", "#fb6a4a", "#de2d26", "#a50f15",  "#000000")) +
  theme(axis.text.x = element_text(angle = -45, hjust = 0)) + 
  labs(x = "", y = "Relative ASV Abundance")

filt_stacked_asvs$order_MC <- factor(filt_stacked_asvs$MC, levels = c("MC2", "MC3", "MC1", "MC4", "MC5", "MC6", "Mt Edwards diffuse fluid", "Background seawater"))
filt_stacked_asvs <- stacked_asvs %>%
  filter(!(Phylum == "Metazoa" | Phylum == "Streptophyta" | Phylum == "Opisthokonta_X" | Phylum == "Stramenopiles_X")) %>%
  unite(TAXA_LEVEL, Supergroup, Phylum, sep = "-", remove = FALSE) %>%
  separate(SAMPLE, into = c("MC", "SUBSTRATE"), sep = "-", remove = FALSE, fill = "left") %>%
  mutate(SITE = case_when(
    SUBSTRATE == "Background seawater" ~ "Background",
    SUBSTRATE == "Mt Edwards diffuse fluid" ~ "Vent",
    .default = SUBSTRATE
  ))
ggplot(filt_stacked_asvs, aes(x = ordered_SAMPLES_18S, y = Count, fill = TAXA_LEVEL)) +
  geom_bar(stat = "identity", position = "fill", color = "black") +
  theme_classic() +
  scale_fill_manual(values = c("#f1eef6", "#d7b5d8", "#df65b0", "#dd1c77", 
                               "#fde0dd", "#fa9fb5", "#c51b8a", "#edf8fb", "#bfd3e6", "#9ebcda", 
                               "#8c96c6", "#8856a7", "#810f7c","#ffffcc", "#f0f9e8", "#c7e9b4", "#7fcdbb", 
                               "#41b6c4", "#2c7fb8", "#253494", "#fee391",  "#fec44f",  "#fe9929", 
                               "#d95f0e", "#993404", "#fc9272", "#fb6a4a", "#de2d26", "#a50f15")) +
  facet_grid(cols = vars(SITE), scales = "free", space = "free") +
  theme(axis.text.x = element_text(angle = -45, hjust = 0)) + 
  labs(x = "", y = "Relative ASV Abundance")

ggplot(filt_stacked_asvs, aes(x = ordered_SAMPLES_18S, y = Count, fill = TAXA_LEVEL)) +
  geom_bar(stat = "identity", position = "fill", color = "black") +
  theme_classic() +
  scale_fill_manual(values = c("#f1eef6", "#d7b5d8", "#df65b0", "#dd1c77", 
                               "#fde0dd", "#fa9fb5", "#c51b8a", "#edf8fb", "#bfd3e6", "#9ebcda", 
                               "#8c96c6", "#8856a7", "#810f7c","#ffffcc", "#f0f9e8", "#c7e9b4", "#7fcdbb", 
                               "#41b6c4", "#2c7fb8", "#253494", "#fee391",  "#fec44f",  "#fe9929", 
                               "#d95f0e", "#993404", "#fc9272", "#fb6a4a", "#de2d26", "#a50f15")) +
  facet_grid(cols = vars(order_MC), scales = "free", space = "free") +
  theme(axis.text.x = element_text(angle = -45, hjust = 0)) + 
  labs(x = "", y = "Relative ASV Abundance")

#### Relative SEQ abundance across sampled sites - stacked bar (w and w/out Streptophyta and Metazoa) ####
# 16S 
rel_seq_16S <- mc_16S_df %>%
  filter(SEQ_AVG > 10) %>%
  separate(SAMPLE, into = c("MC", "SUBSTRATE"), sep = "-", remove = FALSE, fill = "right") %>%
  filter(SUBSTRATE == "Quartz" | SUBSTRATE == "Riftia" | SUBSTRATE == "Shell") %>%
  group_by(SAMPLE, Phylum) %>%
  summarise(SUM_SEQ = sum(SEQ_AVG)) %>%
  separate(SAMPLE, into = c("MC", "SUBSTRATE"), sep = "-", remove = FALSE, fill = "right")

df_taxonomy_16S <- mc_16S_df %>%
  ungroup() %>%
  select(Phylum, Class)%>%
  distinct()
rel_seq_16S_2 <- rel_seq_16S %>%
  left_join(df_taxonomy_16S, by = "Phylum", relationship = "many-to-many") %>%
  unite(TAXA_LEVEL, Phylum, Class, sep = "-", remove = FALSE)
ggplot(rel_seq_16S, aes(x = SAMPLE, y = SUM_SEQ, fill = Phylum)) +
  geom_bar(stat = "identity", position = "fill", color = "black") +
  facet_grid(cols = vars(ordered_MC), scales = "free", space = "free") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = -45, hjust = 0)) + 
  labs(x = "", y = "16S Relative Seqeunce Abundance")
rel_seq_16S$ordered_MC <- factor(rel_seq_16S$MC, levels = c("MC2", "MC3", "MC1", "MC4", "MC5", "MC6"))

#18S
rel_seq_18S_2 <- mc_18S_df %>%
  group_by(SAMPLE, Phylum) %>%
  summarise(SUM_SEQ = sum(SEQ_AVG)) 
main_rel_seq_18S_2 <- rel_seq_18S_2 %>%
  filter(SUM_SEQ > 300) %>%
  separate(SAMPLE, into = c("MC", "SUBSTRATE"), sep = "-", remove = FALSE, fill = "left")

rel_seq_18S_2$ordered_SAMPLES_18S <- factor(rel_seq_18S_2$SAMPLE, levels = c("MC2-Quartz", "MC2-Riftia", "MC3-Riftia", "MC3-Shell", "MC1-Quartz", "MC1-Riftia", "MC1-Shell", "MC4-Quartz", "MC4-Shell", "MC5-Quartz", "MC5-Riftia", "MC5-Shell", "MC6-Quartz", "MC6-Riftia", "MC6-Shell", "Mt Edwards diffuse fluid", "Background seawater"))

ggplot(rel_seq_18S_2, aes(x = ordered_SAMPLES_18S, y = SUM_SEQ, fill = Phylum)) +
  geom_bar(stat = "identity", position = "fill", color = "black") +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = -45, hjust = 0))

ggplot(main_rel_seq_18S_2, aes(x = ordered_SAMPLES_18S, y = SUM_SEQ, fill = Phylum)) +
  geom_bar(stat = "identity", position = "fill", color = "black") +
  facet_grid(cols = vars(SUBSTRATE), scales = "free", space = "free") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = -45, hjust = 0))

df_taxonomy <- mc_18S_df %>%
  ungroup() %>%
  select(Supergroup, Phylum)%>%
  distinct()

filt_rel_seq_18S <- rel_seq_18S_2 %>%
  filter(!(Phylum == "Metazoa" | Phylum == "Streptophyta" | Phylum == "Opisthokonta_X" | Phylum == "Stramenopiles_X")) %>%
  separate(SAMPLE, into = c("MC", "SUBSTRATE"), sep = "-", remove = FALSE, fill = "left") %>%
  mutate(SITE = case_when(
    SUBSTRATE == "Background seawater" ~ "Background",
    SUBSTRATE == "Mt Edwards diffuse fluid" ~ "Vent",
    .default = SUBSTRATE
  )) %>%
  left_join(df_taxonomy, by = "Phylum", relationship = "many-to-many") %>%
  unite(TAXA_LEVEL, Supergroup, Phylum, sep = "-", remove = FALSE)

ggplot(filt_rel_seq_18S, aes(x = ordered_SAMPLES_18S, y = SUM_SEQ, fill = TAXA_LEVEL)) +
  geom_bar(stat = "identity", position = "fill", color = "black") +
  facet_grid(cols = vars(SITE), scales = "free", space = "free") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = -45, hjust = 0)) + 
  labs(x = "", y = "Relative Seqeunce Abundance")
ggplot(filt_rel_seq_18S, aes(x = ordered_SAMPLES_18S, y = SUM_SEQ, fill = TAXA_LEVEL)) +
  geom_bar(stat = "identity", position = "fill", color = "black") +
  facet_grid(cols = vars(MC), scales = "free", space = "free") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = -45, hjust = 0)) + 
  labs(x = "", y = "Relative Seqeunce Abundance")
#### ASV count of each phyla across substrates - point/dot plot ####
ASV_onsubs <- stacked_asvs %>%
  separate(SAMPLE, into = c("MC", "SUBSTRATE"), sep = "-", remove = FALSE) %>%
  filter(SUBSTRATE == "Quartz" | SUBSTRATE == "Riftia" | SUBSTRATE == "Shell") %>%
  filter(!(Phylum == "Opisthokonta_X" | Phylum == "Stramenopiles_X"))
ggplot(ASV_onsubs, aes(x = SUBSTRATE, y = Phylum)) +
  geom_point(aes(size = Count), color = "blue", alpha= 0.7) +
  facet_grid(rows = vars(Supergroup), scales = "free", space = "free") + 
  scale_size_continuous(name = "ASV Count", range = c(2, 10)) +
  theme_classic()
# w out Meta/Strep
ASV_onsubs <- stacked_asvs %>%
  separate(SAMPLE, into = c("MC", "SUBSTRATE"), sep = "-", remove = FALSE) %>%
  filter(SUBSTRATE == "Quartz" | SUBSTRATE == "Riftia" | SUBSTRATE == "Shell") %>%
  filter(!(Phylum == "Opisthokonta_X" | Phylum == "Stramenopiles_X" | Phylum == "Metazoa" | Phylum == "Streptophyta"))
ggplot(ASV_onsubs, aes(x = SUBSTRATE, y = Phylum)) +
  geom_point(aes(size = Count), color = "blue", alpha= 0.7) +
  facet_grid(rows = vars(Supergroup), scales = "free", space = "free") + 
  scale_size_continuous(name = "ASV Count", range = c(2, 10)) +
  theme_classic()
#### Relative SEQ abundance of each phylum across MCs - tile plot (grouped by substrate and supergroup, no Metazoa or Streptophyta) ####
seq_rel_ab_18S<- mc_df_18s %>%
  mutate(Relative_Abundance = SEQ_AVG / sum(SEQ_AVG))  

seq_rel_ab_18S$MC_ordered_18S <- factor(seq_rel_ab_18S$MC, levels = c("MC2", "MC3", "MC1", "MC4", "MC5", "MC6"))

ggplot(seq_rel_ab_18S, aes(x = MC_ordered_18S, y = Phylum, fill = Relative_Abundance)) +
  geom_tile(stat = "identity") +
  theme_classic() +
  facet_grid(cols = vars(SUBSTRATE), rows = vars(Supergroup), scales = "free", space = "free") +
  scale_fill_gradient(low = "#c6dbef",
                      high = "#08519c",
                      guide = "colorbar")


#### 2. Alpha and beta diversity ####
### ALL SAMPLED SITES ###
#18S-Shannon
wide_18S <- mc_18S_df %>% ungroup() %>%
  select(SAMPLE, FeatureID, SEQ_AVG) %>%
  pivot_wider(names_from = SAMPLE, values_from = SEQ_AVG, values_fill = 0) %>%
  column_to_rownames(var = "FeatureID")
wide_18S_mat <- as.matrix(wide_18S)
shannon_18S <- diversity(wide_18S_mat, index = "shannon", MARGIN = 2)
plot_shannon_18S <- as.data.frame(shannon_18S) %>%
  rownames_to_column(var = "SAMPLE")
sep_plot_shannon_18S <- plot_shannon_18S %>%
  separate(SAMPLE, into = c("MC", "SUBSTRATE"), sep = "-", remove = FALSE)
sep_plot_shannon_18S
unique(sep_plot_shannon_18S$SAMPLE)
sep_plot_shannon_18S$ordered_SAMPLES_18S <- factor(sep_plot_shannon_18S$SAMPLE, levels = c("MC2-Quartz", "MC2-Riftia", "MC3-Riftia", "MC3-Shell", "MC1-Quartz", "MC1-Riftia", "MC1-Shell", "MC4-Quartz", "MC4-Shell", "MC5-Quartz", "MC5-Riftia", "MC5-Shell", "MC6-Quartz", "MC6-Riftia", "MC6-Shell", "Mt Edwards diffuse fluid", "Background seawater"))

ggplot(sep_plot_shannon_18S, aes(x = ordered_SAMPLES_18S, y = shannon_18S), fill = MC, shape = SUBSTRATE) +
  geom_point(fill = c( "#2a2b47", "#8ac926","#8ac926", "#8ac926", "#ff595e", "#ff595e", "#ff924c", "#ff924c", "#8ac926", "#8ac926", "#1982c4", "#1982c4", "#1982c4", "#6a4c93", "#6a4c93",  "#6a4c93", "#2a2b47"), shape = c(4, 21, 24, 23, 21, 24, 24, 23, 21, 23, 21, 24, 23, 21, 24, 23, 4)) +
  theme(axis.text.x = element_text(angle = -45, hjust = 0)) +
  labs(x= "", y = "Shannon diversity")

onlyMCs_shann <- sep_plot_shannon_18S %>%
  filter(!(SAMPLE == "Background seawater"| SAMPLE == "Mt Edwards diffuse fluid"))
ggplot(onlyMCs_shann, aes(x = ordered_SAMPLES_18S, y = shannon_18S), fill = MC, shape = SUBSTRATE) +
  geom_point(fill = c( "#8ac926","#8ac926", "#8ac926", "#ff595e", "#ff595e", "#ff924c", "#ff924c", "#8ac926", "#8ac926", "#1982c4", "#1982c4", "#1982c4", "#6a4c93", "#6a4c93",  "#6a4c93"), shape = c(21, 24, 23, 21, 24, 24, 23, 21, 23, 21, 24, 23, 21, 24, 23)) +
  theme(axis.text.x = element_text(angle = -45, hjust = 0)) +
  labs(x= "", y = "Shannon diversity")

# 18S-Inv Simp
stand_wide_18S_mat <- decostand(wide_18S_mat, method = "hellinger", MARGIN = 2)

inv_simp_18S <- diversity(stand_wide_18S_mat, index = "invsimpson", MARGIN = 2, equalize.groups = TRUE)
inv_simp_18S
plot_simpson_18S <- as.data.frame(inv_simp_18S) %>%
  rownames_to_column(var = "SAMPLE")

sep_plot_simpson_18S <- plot_simpson_18S %>%
  separate(SAMPLE, into = c("MC", "SUBSTRATE"), sep = "-", remove = FALSE)

filt_sep_plot_simpson_18S <- sep_plot_simpson_18S %>%
  filter(!(SAMPLE == "Background seawater" | SAMPLE == "Mt Edwards diffuse fluid"))
sep_plot_simpson_18S$ordered_SAMPLES_18S <- factor(sep_plot_simpson_18S$SAMPLE, levels = c("MC2-Quartz", "MC2-Riftia", "MC3-Riftia", "MC3-Shell", "MC1-Quartz", "MC1-Riftia", "MC1-Shell", "MC4-Quartz", "MC4-Shell", "MC5-Quartz", "MC5-Riftia", "MC5-Shell", "MC6-Quartz", "MC6-Riftia", "MC6-Shell", "Mt Edwards diffuse fluid", "Background seawater"))
filt_sep_plot_simpson_18S$ordered_SAMPLES_18S <- factor(filt_sep_plot_simpson_18S$SAMPLE, levels = c("MC2-Quartz", "MC2-Riftia", "MC3-Riftia", "MC3-Shell", "MC1-Quartz", "MC1-Riftia", "MC1-Shell", "MC4-Quartz", "MC4-Shell", "MC5-Quartz", "MC5-Riftia", "MC5-Shell", "MC6-Quartz", "MC6-Riftia", "MC6-Shell"))

ggplot(sep_plot_simpson_18S, aes(x = ordered_SAMPLES_18S, y = inv_simp_18S), fill = MC, shape = SUBSTRATE) +
  geom_point(fill = c( "#2a2b47", "#8ac926","#8ac926", "#8ac926", "#ff595e", "#ff595e", "#ff924c", "#ff924c", "#8ac926", "#8ac926", "#1982c4", "#1982c4", "#1982c4", "#6a4c93", "#6a4c93",  "#6a4c93", "#2a2b47"), shape = c(4, 21, 24, 23, 21, 24, 24, 23, 21, 23, 21, 24, 23, 21, 24, 23, 4)) +
  theme(axis.text.x = element_text(angle = -45, hjust = 0)) +
  labs(x= "", y = "Inverse Simpson diversity")

ggplot(filt_sep_plot_simpson_18S, aes(x = ordered_SAMPLES_18S, y = inv_simp_18S), fill = MC, shape = SUBSTRATE) +
  geom_point(fill = c( "#8ac926","#8ac926", "#8ac926", "#ff595e", "#ff595e", "#ff924c", "#ff924c", "#8ac926", "#8ac926", "#1982c4", "#1982c4", "#1982c4", "#6a4c93", "#6a4c93",  "#6a4c93"), shape = c(21, 24, 23, 21, 24, 24, 23, 21, 23, 21, 24, 23, 21, 24, 23)) +
  theme(axis.text.x = element_text(angle = -45, hjust = 0)) +
  labs(x= "", y = "Shannon diversity")
#
ggplot(plot_shannon_microcol_18S, aes(x= SAMPLE, y= shannon_microcol_18s)) +
  geom_point(color = "purple", size = 4, shape = 17) +
  labs(x = "", y= "Shannon diversity") +
  theme(axis.text.x = element_text(angle = -45, hjust = 0))

ggplot(plot_inv_simp_microcol_18S, aes(x= MC, y= inv_simp_microcol_18S)) +
  geom_point(color = "#2d00f7", size = 4, shape = 17) +
  labs(x= "", y= "Inverse Simpson diversity") +
  theme(axis.text.x = element_text(angle = -45, hjust = 0))
#
### 16S
wide_microcol_only_16S <- microcol_only_16S %>% ungroup() %>%
  select(SAMPLE, FeatureID, SUM) %>%
  pivot_wider(names_from = SAMPLE, values_from = SUM, values_fill = 0) %>%
  column_to_rownames(var = "FeatureID")
wide_microcol_only_16S_mat <- as.matrix(wide_microcol_only_16S)  
# 16S-Shannon 
shannon_microcol_16s <- diversity(wide_microcol_only_16S_mat, index = "shannon", MARGIN = 2)
shannon_microcol_16s
plot_shannon_microcol_16S <- as.data.frame(shannon_microcol_16s) %>%
  rownames_to_column(var = "SAMPLE") %>%
  separate(col = SAMPLE, into = c("MC", "SUBSTRATE"), sep = "-", remove = FALSE)

# 16S-Inv Simp
stand_wide_microcol_16S_mat<- decostand(wide_microcol_only_16S_mat, method = "hellinger", MARGIN = 2)
inv_simp_microcol_16S <- diversity(stand_wide_microcol_16S_mat, index = "invsimpson", MARGIN = 2, equalize.groups = TRUE)
plot_inv_simp_microcol_16S
plot_inv_simp_microcol_16S <- as.data.frame(inv_simp_microcol_16S) %>%
  rownames_to_column(var = "SAMPLE") %>%
  separate(col = SAMPLE, into = c("MC", "SUBSTRATE"), sep = "-", remove = FALSE) %>%
  filter(!(SUBSTRATE == "Basalt"| SUBSTRATE == "Olivine" | SUBSTRATE == "Pyrite"))
unique(plot_inv_simp_microcol_16S$SAMPLE)
plot_inv_simp_microcol_16S$ordered_SAMPLES_16S <- factor(plot_inv_simp_microcol_16S$SAMPLE, levels = c("MC2-Quartz", "MC2-Riftia", "MC3-Quartz", "MC3-Riftia", "MC3-Shell", "MC1-Quartz", "MC1-Riftia", "MC1-Shell", "MC4-Quartz", "MC4-Riftia", "MC5-Quartz", "MC5-Riftia", "MC5-Shell", "MC6-Quartz", "MC6-Riftia", "MC6-Shell"))

ggplot(plot_inv_simp_microcol_16S, aes(x = ordered_SAMPLES_16S, y = inv_simp_microcol_16S), fill = MC, shape = SUBSTRATE) +
  geom_point(fill = c( "#8ac926","#8ac926", "#8ac926", "#ff595e", "#ff595e", "#ff924c", "#ff924c", "#ff924c", "#8ac926", "#8ac926", "#1982c4", "#1982c4", "#1982c4", "#6a4c93", "#6a4c93",  "#6a4c93"), shape = c(21, 24, 23, 21, 24, 21, 24, 23, 21, 23, 21, 24, 23, 21, 24, 23)) +
  theme(axis.text.x = element_text(angle = -45, hjust = 0)) +
  labs(x= "", y = "Shannon diversity")
#
ggplot(plot_shannon_microcol_16S, aes(x= SAMPLE, y= shannon_microcol_16s)) +
  geom_point(color = "purple", size = 4, shape = 17) +
  labs(x = "", y= "16S Shannon diversity") +
  theme(axis.text.x = element_text(angle = -45, hjust = 0))
ggplot(plot_inv_simp_microcol_16S, aes(x= SAMPLE, y= inv_simp_microcol_16S)) +
  geom_point(color = "#2d00f7", size = 4, shape = 17) +
  labs(x= "", y= "16S Inverse Simpson diversity") +
  theme(axis.text.x = element_text(angle = -45, hjust = 0))
#
### Total ASVs / Sp richness (B/MC/V, MCs, SUBS)###
# 18S-B/MC/V
samp_sites_asvs_18S <- mc_18S_df %>%
  mutate(SAMPLE_SITE = case_when(
    SAMPLE == "Background seawater" ~ "Background",
    SAMPLE == "Mt Edwards diffuse fluid" ~ "Vent",
    grepl("MC", SAMPLE) ~ "MC",
    TRUE ~ "Other"))
unique(samp_sites_asvs_18S$SAMPLE_SITE)  
plot_samp_asvs_18s <- samp_sites_asvs_18S %>%
  group_by(SAMPLE_SITE) %>%
  summarise(Total_ASVs = n_distinct(FeatureID))

ggplot(plot_samp_asvs_18s, aes(x= SAMPLE_SITE, y= Total_ASVs)) +
  geom_point(color = "#3eabf4", size = 4, shape = 18) +
  labs(x= "Sample Site", y= "Total ASVs")

# 18S-MCs
microcol_only <- mc_18S_df %>%
  separate(SAMPLE, into = c("SAMPLE_ID", "SUBSTRATE"), sep = "-", remove = FALSE) %>%
  filter(!is.na(SUBSTRATE)) %>%
  group_by(SAMPLE, SUBSTRATE, FeatureID) %>%
  summarise(SUM = sum(SEQ_AVG))
sep_microcol_only <- microcol_only %>%
  separate(SAMPLE, into = c("MC", "SUBSTRATE"), sep = "-", remove = FALSE)
plot_sub_asvs_18S <- sep_microcol_only %>%
  group_by(SUBSTRATE) %>%
  summarise(Total_ASVs = n_distinct(FeatureID))
ggplot(plot_sub_asvs_18S, aes(x= SUBSTRATE, y= Total_ASVs)) +
  geom_point(color = "#3eabf4", size = 4, shape = 15) +
  labs(x= "", y= "Total ASVs")

# 18S-Q/R/S
plot_sub_asvs_18S <- sep_microcol_only %>%
  group_by(SUBSTRATE) %>%
  summarise(Total_ASVs = n_distinct(FeatureID))
ggplot(plot_sub_asvs_18S, aes(x= SUBSTRATE, y= Total_ASVs)) +
  geom_point(color = "#3eabf4", size = 4, shape = 15) +
  labs(x= "", y= "Total ASVs")

## 16S-B/MC/V
MC_only_16S <- mc_16S_df %>%
  filter(SAMPLE != "Background seawater" & SAMPLE != "Mt Edwards diffuse fluid") %>%
  separate(SAMPLE, into = c("MC", "SUBSTRATE"), sep = "-") %>%
  filter(SUBSTRATE == "Quartz" | SUBSTRATE == "Riftia" | SUBSTRATE == "Shell")
MC_only_16S$MC_ordered_16S <- factor(MC_only_16S$MC, levels = c("MC2", "MC3", "MC1", "MC4", "MC5", "MC6"))

plot_mc_asvs_16S <- MC_only_16S %>%
  group_by(MC_ordered_16S) %>%
  summarise(Total_ASVs = n_distinct(FeatureID))
ggplot(plot_mc_asvs_16S, aes(x= MC_ordered_16S, y= Total_ASVs)) +
  geom_point(color = "#3eabf4", size = 4, shape = 17) +
  labs(x= "", y= "Total ASVs")
# 16S-MCs
MC_only_16S <- mc_16S_df %>%
  filter(SAMPLE != "Background seawater" & SAMPLE != "Mt Edwards diffuse fluid") %>%
  separate(SAMPLE, into = c("MC", "SUBSTRATE"), sep = "-") %>%
  filter(SUBSTRATE == "Quartz" | SUBSTRATE == "Riftia" | SUBSTRATE == "Shell")

MC_only_16S$MC_ordered_16S <- factor(MC_only_16S$MC, levels = c("MC2", "MC3", "MC1", "MC4", "MC5", "MC6"))

plot_mc_asvs_16S <- MC_only_16S %>%
  group_by(MC_ordered_16S) %>%
  summarise(Total_ASVs = n_distinct(FeatureID))
ggplot(plot_mc_asvs_16S, aes(x= MC_ordered_16S, y= Total_ASVs)) +
  geom_point(color = "#3eabf4", size = 4, shape = 17) +
  labs(x= "", y= "Total ASVs")

# 16S-Q/R/S
plot_sub_asvs_18S <- sep_microcol_only %>%
  group_by(SUBSTRATE) %>%
  summarise(Total_ASVs = n_distinct(FeatureID))
ggplot(plot_sub_asvs_18S, aes(x= SUBSTRATE, y= Total_ASVs)) +
  geom_point(color = "#3eabf4", size = 4, shape = 15) +
  labs(x= "", y= "Total ASVs")

### Shannon div (B/MC/V, MCs1-6, Q/R/S) ###
# 18S-B/MC/V
sample_sites_18S <- mc_18S_df %>% 
  mutate(SAMPLE_SITE = case_when(
    SAMPLE == "Background seawater" ~ "Background",
    SAMPLE == "Mt Edwards diffuse fluid" ~ "Vent",
    grepl("MC", SAMPLE) ~ "MC",
    TRUE ~ "Other"))
summ_samp_sit_18S <- sample_sites_18S %>%
  group_by(SAMPLE_SITE, FeatureID) %>%
  summarize(SUM = sum(SEQ_AVG))

wide_samp_sites_18s <- summ_samp_sit_18S %>% ungroup() %>%
  select(SAMPLE_SITE, FeatureID, SUM) %>%
  pivot_wider(names_from = SAMPLE_SITE, values_from = SUM, values_fill = 0) %>%
  column_to_rownames(var = "FeatureID")

wide_samp_sites_18s_mat <- as.matrix(wide_samp_sites_18s)
stand_wide_samp_site_18S_mat <- decostand(wide_samp_sites_18s_mat, method = "hellinger", MARGIN = 2)
shannon_samp_18S <- diversity(wide_samp_sites_18s_mat, index = "shannon", MARGIN = 2)
shannon_samp_18S
plot_shannon_samp_18S <- as.data.frame(shannon_samp_18S) %>%
  rownames_to_column(var = "SAMPLE_SITE")
plot_shannon_samp_18S  

ggplot(plot_shannon_samp_18S, aes(x= SAMPLE_SITE, y= shannon_samp_18S)) +
  geom_point(color = "purple", size = 4, shape = 18) +
  labs(x = "Sample Site", y= "Shannon diversity")
# 18S-MCs
MC_shann <- sep_microcol_only %>%
  group_by(MC, FeatureID) %>%
  summarize(SUM = sum(SUM))
wide_MCs <- MC_shann %>% ungroup() %>%
  select(MC, FeatureID, SUM) %>%
  pivot_wider(names_from = MC, values_from = SUM, values_fill = 0) %>%
  column_to_rownames(var = "FeatureID")
wide_MCs_mat <- as.matrix(wide_MCs)  
shannon_MCs_18s <- diversity(wide_MCs_mat, index = "shannon", MARGIN = 2)
shannon_MCs_18s
plot_shannon_MCs_18S <- as.data.frame(shannon_MCs_18s) %>%
  rownames_to_column(var = "MC") 
ggplot(plot_shannon_MCs_18S, aes(x= MC, y= shannon_MCs_18s)) +
  geom_point(color = "purple", size = 4, shape = 17) +
  labs(x = "", y= "Shannon diversity")

# 18s-Q/R/S
grp_subs_18s <- sep_microcol_only %>%
  group_by(SUBSTRATE, FeatureID) %>%
  summarize(SUM = sum(SUM))
grp_subs_18s
wide_subs_18s <- grp_subs_18s %>% ungroup() %>%
  select(SUBSTRATE, FeatureID, SUM) %>%
  pivot_wider(names_from = SUBSTRATE, values_from = SUM, values_fill = 0) %>%
  column_to_rownames(var = "FeatureID")
wide_subs_18s_mat <- as.matrix(wide_subs_18s)  
shannon_subs_18s <- diversity(wide_subs_18s_mat, index = "shannon", MARGIN = 2)
shannon_subs_18s
plot_shannon_subs_18s <- as.data.frame(shannon_subs_18s) %>%
  rownames_to_column(var = "SUBSTRATE")
ggplot(plot_shannon_subs_18s, aes(x= SUBSTRATE, y= shannon_subs_18s)) +
  geom_point(color = "purple", size = 4, shape = 15) +
  labs(x = "", y= "Shannon diversity")

# 16S-B/MC/V
sample_sites_16S <- mc_16S_df %>% 
  mutate(SAMPLE_SITE = case_when(
    SAMPLE == "Background seawater" ~ "Background",
    SAMPLE == "Mt Edwards diffuse fluid" ~ "Vent",
    grepl("MC", SAMPLE) ~ "MC",
    TRUE ~ "Other"))
summ_samp_sit_16S <- sample_sites_16S %>%
  group_by(SAMPLE_SITE, FeatureID) %>%
  summarize(SUM = sum(SEQ_AVG))

wide_samp_sites_16s <- summ_samp_sit_16S %>% ungroup() %>%
  select(SAMPLE_SITE, FeatureID, SUM) %>%
  pivot_wider(names_from = SAMPLE_SITE, values_from = SUM, values_fill = 0) %>%
  column_to_rownames(var = "FeatureID")

wide_samp_sites_16s_mat <- as.matrix(wide_samp_sites_16s)
stand_wide_samp_site_16S_mat <- decostand(wide_samp_sites_16s_mat, method = "hellinger", MARGIN = 2)
shannon_samp_16S <- diversity(wide_samp_sites_16s_mat, index = "shannon", MARGIN = 2)
shannon_samp_16S
plot_shannon_samp_16S <- as.data.frame(shannon_samp_16S) %>%
  rownames_to_column(var = "SAMPLE_SITE")
plot_shannon_samp_16S
ggplot(plot_shannon_samp_16S, aes(x= SAMPLE_SITE, y= shannon_samp_16S)) +
  geom_point(color = "purple", size = 4, shape = 18) +
  labs(x = "Sample Site", y= "16S Shannon diversity")

# 16S-MCs
MC_only_16S <- mc_16S_df %>%
  separate(SAMPLE, into = c("MC", "SUBSTRATE"), sep = "-", remove = FALSE) %>%
  filter(!is.na(SUBSTRATE)) %>%
  group_by(MC, FeatureID) %>%
  summarise(SUM = sum(SEQ_AVG), .groups = 'drop')
wide_MC_only_16S <- MC_only_16S %>% ungroup() %>%
  select(MC, FeatureID, SUM) %>%
  pivot_wider(names_from = MC, values_from = SUM, values_fill = 0) %>%
  column_to_rownames(var = "FeatureID")
wide_MC_only_16S_mat <- as.matrix(wide_MC_only_16S)  
shannon_MC_16s <- diversity(wide_MC_only_16S_mat, index = "shannon", MARGIN = 2)
shannon_MC_16s
plot_shannon_MC_16S <- as.data.frame(shannon_MC_16s) %>%
  rownames_to_column(var = "MC")
plot_shannon_MC_16S$MC_ordered_16S <- factor(plot_shannon_MC_16S$MC, levels = c("MC2", "MC3", "MC1", "MC4", "MC5", "MC6"))

ggplot(plot_shannon_MC_16S, aes(x= MC_ordered_16S, y= shannon_MC_16s)) +
  geom_point(color = "purple", size = 4, shape = 17) +
  labs(x = "", y= "16S Shannon diversity")

# 16S-Q/R/S
grp_subs_16S <- mc_16S_df %>%
  separate(SAMPLE, into = c("MC", "SUBSTRATE"), sep = "-", remove = FALSE) %>%
  filter(!is.na(SUBSTRATE)) %>%
  group_by(SUBSTRATE, FeatureID) %>%
  summarise(SUM = sum(SEQ_AVG))

wide_subs_16S <- grp_subs_16S %>% ungroup() %>%
  select(SUBSTRATE, FeatureID, SUM) %>%
  pivot_wider(names_from = SUBSTRATE, values_from = SUM, values_fill = 0) %>%
  column_to_rownames(var = "FeatureID")
wide_subs_16S_mat <- as.matrix(wide_subs_16S)  

ggplot(plot_shannon_subs_16S, aes(x= SUBSTRATE, y= shannon_subs_16S)) +
  geom_point(color = "purple", size = 4, shape = 15) +
  labs(x = "", y= "16S Shannon diversity")

### Inv Simp div (B/MC/V, MCs, Q/R/S)
# 18S-B/MC/V
inv_simp_samp_18S <- diversity(stand_wide_samp_site_18S_mat, index = "invsimpson", MARGIN = 2, equalize.groups = TRUE)
plot_inv_simp_samp_18S
plot_inv_simp_samp_18S <- as.data.frame(inv_simp_samp_18S) %>%
  rownames_to_column(var = "SAMPLE_SITE")

ggplot(plot_inv_simp_samp_18S, aes(x= SAMPLE_SITE, y= inv_simp_samp_18S)) +
  geom_point(color = "#2d00f7", size = 4, shape = 18) +
  labs(x= "Sample Site", y= "Inverse Simpson diversity")

# 18S-MCs
stand_wide_MCs_18S_mat<- decostand(wide_MCs_mat, method = "hellinger", MARGIN = 2)
inv_simp_MCs_18S <- diversity(stand_wide_MCs_18S_mat, index = "invsimpson", MARGIN = 2, equalize.groups = TRUE)
plot_inv_simp_MCs_18S
plot_inv_simp_MCs_18S <- as.data.frame(inv_simp_MCs_18S) %>%
  rownames_to_column(var = "MC")
ggplot(plot_inv_simp_MCs_18S, aes(x= MC, y= inv_simp_MCs_18S)) +
  geom_point(color = "#2d00f7", size = 4, shape = 17) +
  labs(x= "", y= "Inverse Simpson diversity")

# 18S-Q/R/S
stand_wide_subs_18s_mat<- decostand(wide_subs_18s_mat, method = "hellinger", MARGIN = 2)
inv_simp_subs_18s <- diversity(stand_wide_subs_18s_mat, index = "invsimpson", MARGIN = 2, equalize.groups = TRUE)
inv_simp_subs_18s
plot_inv_simp_subs_18s <- as.data.frame(inv_simp_subs_18s) %>%
  rownames_to_column(var = "SUBSTRATE")
ggplot(plot_inv_simp_subs_18s, aes(x= SUBSTRATE, y= inv_simp_subs_18s)) +
  geom_point(color = "#2d00f7", size = 4, shape = 15) +
  labs(x= "", y= "Inverse Simpson diversity")

# 16S-B/MC/V
inv_simp_samp_16S <- diversity(stand_wide_samp_site_16S_mat, index = "invsimpson", MARGIN = 2, equalize.groups = TRUE)
inv_simp_samp_16S
plot_inv_simp_samp_16S <- as.data.frame(inv_simp_samp_16S) %>%
  rownames_to_column(var = "SAMPLE_SITE")
ggplot(plot_inv_simp_samp_16S, aes(x= SAMPLE_SITE, y= inv_simp_samp_16S)) +
  geom_point(color = "#2d00f7", size = 4, shape = 18) +
  labs(x= "Sample Site", y= "Inverse Simpson diversity")

# 16S-MCs
stand_wide_MC_16S_mat<- decostand(wide_MC_only_16S_mat, method = "hellinger", MARGIN = 2)
inv_simp_MC_16S <- diversity(stand_wide_MC_16S_mat, index = "invsimpson", MARGIN = 2, equalize.groups = TRUE)
inv_simp_MC_16S
plot_inv_simp_MC_16S <- as.data.frame(inv_simp_MC_16S) %>%
  rownames_to_column(var = "MC")
ggplot(plot_inv_simp_MC_16S, aes(x= MC, y= inv_simp_MC_16S)) +
  geom_point(color = "#2d00f7", size = 4, shape = 17) +
  labs(x= "", y= "16S Inverse Simpson diversity")

# 16S-Q/R/S
stand_wide_subs_16S_mat<- decostand(wide_subs_16S_mat, method = "hellinger", MARGIN = 2)
inv_simp_subs_16S <- diversity(stand_wide_subs_16S_mat, index = "invsimpson", MARGIN = 2, equalize.groups = TRUE)
inv_simp_subs_16S
plot_inv_simp_subs_16S <- as.data.frame(inv_simp_subs_16S) %>%
  rownames_to_column(var = "SUBSTRATE")




#### Microbes that might have been recruited from the vent and settled on substrate ####
## Find FeatureIDs that are present in both locations
#18S
df_18s_vs <- mc_tmp_18 %>%
  left_join(key_dist_18) %>%
  filter(DISTRIBUTION == "Vent & Substrate only")
sep_18s_vs <- df_18s_vs %>%
  separate(col = SAMPLE, into = c("MC", "SUBSTRATE"), sep = "-", remove = FALSE, fill = "right") %>%
  filter(Phylum != "Metazoa")
shared_MC2 <- sep_18s_vs %>% 
  filter(MC %in% c("MC2", "Mt Edwards diffuse fluid")) %>%
  group_by(FeatureID) %>%
  filter(n_distinct(MC) == 2) %>%  # Keeps only FeatureIDs appearing in both MC2 and vent
  ungroup()
shared_MC3 <- sep_18s_vs %>%
  filter(MC %in% c("MC3", "Mt Edwards diffuse fluid")) %>%
  group_by(FeatureID) %>%
  filter(n_distinct(MC) == 2) %>%  # Keeps only FeatureIDs appearing in both MC2 and vent
  ungroup()
shared_MC1 <- sep_18s_vs %>%
  filter(MC %in% c("MC1", "Mt Edwards diffuse fluid")) %>%
  group_by(FeatureID) %>%
  filter(n_distinct(MC) == 2) %>%  # Keeps only FeatureIDs appearing in both MC2 and vent
  ungroup()
shared_MC4 <- sep_18s_vs %>%
  filter(MC %in% c("MC4", "Mt Edwards diffuse fluid")) %>%
  group_by(FeatureID) %>%
  filter(n_distinct(MC) == 2) %>%  # Keeps only FeatureIDs appearing in both MC2 and vent
  ungroup()
shared_MC5 <- sep_18s_vs %>%
  filter(MC %in% c("MC5", "Mt Edwards diffuse fluid")) %>%
  group_by(FeatureID) %>%
  filter(n_distinct(MC) == 2) %>%  # Keeps only FeatureIDs appearing in both MC2 and vent
  ungroup()
shared_MC6 <- sep_18s_vs %>%
  filter(MC %in% c("MC6", "Mt Edwards diffuse fluid")) %>%
  group_by(FeatureID) %>%
  filter(n_distinct(MC) == 2) %>%  # Keeps only FeatureIDs appearing in both MC2 and vent
  ungroup()

# 16S
sep_16s_vs <- df_16s_all %>%
  separate(col = SAMPLE, into = c("MC", "SUBSTRATE"), sep = "-", remove = FALSE, fill = "right")

shared_MC2_16s <- sep_16s_vs %>%
  filter(MC %in% c("MC2", "Mt Edwards diffuse fluid")) %>%
  group_by(FeatureID) %>%
  filter(n_distinct(MC) == 2) %>%  # Keeps only FeatureIDs appearing in both MC2 and vent
  ungroup()

shared_MC3_16s <- sep_16s_vs %>%
  filter(MC %in% c("MC3", "Mt Edwards diffuse fluid")) %>%
  group_by(FeatureID) %>%
  filter(n_distinct(MC) == 2) %>%  # Keeps only FeatureIDs appearing in both MC2 and vent
  ungroup()
shared_MC1_16s <- sep_16s_vs %>%
  filter(MC %in% c("MC1", "Mt Edwards diffuse fluid")) %>% 
  group_by(FeatureID) %>%
  filter(n_distinct(MC) == 2) %>%  # Keeps only FeatureIDs appearing in both MC2 and vent
  ungroup()
shared_MC4_16s <- sep_16s_vs %>%
  filter(MC %in% c("MC4", "Mt Edwards diffuse fluid")) %>%
  group_by(FeatureID) %>%
  filter(n_distinct(MC) == 2) %>%  # Keeps only FeatureIDs appearing in both MC2 and vent
  ungroup()
shared_MC5_16s <- sep_16s_vs %>%
  filter(MC %in% c("MC5", "Mt Edwards diffuse fluid")) %>% 
  group_by(FeatureID) %>%
  filter(n_distinct(MC) == 2) %>%  # Keeps only FeatureIDs appearing in both MC2 and vent
  ungroup()
shared_MC6_16s <- sep_16s_vs %>%
  filter(MC %in% c("MC6", "Mt Edwards diffuse fluid")) %>%
  group_by(FeatureID) %>%
  filter(n_distinct(MC) == 2) %>%  # Keeps only FeatureIDs appearing in both MC2 and vent
  ungroup()



#### RESULTS: Project stats ####
## 1. Total sequences
#16S - total, filtered, unassigned
total_16s_sequences <- sum(asv16s_df$SEQUENCE_COUNT)
total_16s_sequences

filt_asv16s_df <- asv16s_df %>%
  filter(Domain == "Bacteria" | Domain == "Archaea")
sum_16s_seq <- sum(filt_asv16s_df$SEQUENCE_COUNT)
sum_16s_seq

unassigned_16s <- asv16s_df %>%
  filter(Domain == "Unassigned")
sum(unassigned_16s$SEQUENCE_COUNT) / sum(filt_asv16s_df$SEQUENCE_COUNT) * 100

#18S - total, filtered, unassigned 
total_18s_sequences <- sum(asv_wtax_18$SEQUENCE_COUNT)
total_18s_sequences

filt_asv18s <- asv_wtax_18 %>%
  filter(Domain == "Eukaryota")
sum(filt_asv18s$SEQUENCE_COUNT)

unassigned_18S <- asv_wtax_18 %>%
  filter(Domain == "Unassigned")
sum(unassigned_18S$SEQUENCE_COUNT) / sum(asv_wtax_18$SEQUENCE_COUNT) * 100
sum(unassigned_18S$SEQUENCE_COUNT) / sum(filt_asv18s$SEQUENCE_COUNT) * 100

## 2. Total ASVs
# 18S - total, filtered, unassigned
length(unique(asv_wtax_18$FeatureID))
length(unique(filt_asv18s$FeatureID))
length(unique(unassigned_18S$FeatureID)) / length(unique(filt_asv18s$FeatureID)) * 100

# 16S - total, filtered, unassigned 
total_16s_ASVs <- length(unique(asv16s_df$FeatureID))
length(unique(filt_asv16s_df$FeatureID))
length(unique(unassigned_16s$FeatureID)) / length(unique(asv16s_df$FeatureID)) * 100
length(unique(unassigned_16s$FeatureID)) / length(unique(filt_asv16s_df$FeatureID)) * 100

## 3. ASV Distribution on MC substrates
### 18S
# Total ASVs - real total(contains repeats) & unique ASVs
length(df_18s_all$FeatureID)
length(unique(df_18s_all$FeatureID))

# ASVs on substrates (MCs)
ASVs_MC_18s <- df_18s_all %>%
  separate(SAMPLE, into = c("MC", "SUBSTRATE"), sep = "-", remove = FALSE) %>%
  filter(!(is.na(SUBSTRATE)))
length(ASVs_MC_18s$FeatureID)
length(unique(ASVs_MC_18s$FeatureID))

# ASVs on substrate only (MC only - no B/V)
df_18s_subs <- mc_tmp_18 %>%
  left_join(key_dist_18) %>%
  filter(DISTRIBUTION == "Substrate only")
length(df_18s_subs$FeatureID)
length(unique(df_18s_subs$FeatureID))

# MC only ASVs on more than one substrate
asvs_subs_only_18S <- df_18s_subs %>%
  separate(SAMPLE, into = c("MC", "SUBSTRATE"), sep = "-", remove = FALSE)

format_18s_subonly <- asvs_subs_only_18S %>%
  select(FeatureID, SEQ_AVG, SUBSTRATE) %>%
  pivot_wider(id_cols = FeatureID, names_from = SUBSTRATE, values_fn = mean, values_from = SEQ_AVG, values_fill = NA)

format_18s_subonly_assigned <- format_18s_subonly %>%
  mutate(Distribution = case_when(
    (is.na(Riftia) & is.na(Shell) & !(is.na(Quartz))) ~ "Quartz only", 
    (is.na(Riftia) & !(is.na(Shell)) & !(is.na(Quartz))) ~ "Quartz & Shell only", 
    (is.na(Riftia) & !(is.na(Shell)) & is.na(Quartz)) ~ "Shell only", 
    (!(is.na(Riftia)) & is.na(Shell) & !(is.na(Quartz))) ~ "Quartz & Riftia only",
    (!(is.na(Riftia)) & !(is.na(Shell)) & is.na(Quartz)) ~ "Shell & Riftia only", 
    (!(is.na(Riftia)) & is.na(Shell) & is.na(Quartz)) ~ "Riftia only", 
    (!(is.na(Riftia)) & !(is.na(Shell)) & !(is.na(Quartz)) ~ "All")
  ))

key_sub_dist_18 <- format_18s_subonly_assigned %>%
  select(FeatureID, Distribution) %>%
  distinct()
unique(key_sub_dist_18$Distribution)

more_than_one <- key_sub_dist_18 %>%
  filter(!(Distribution == "Quartz only" | Distribution == "Riftia only" | Distribution == "Shell only"))
length(more_than_one$FeatureID)

# MC only ASVs on shell only
shell_only_18s <- key_sub_dist_18 %>%
  filter(Distribution == "Shell only")
length(shell_only_18s$FeatureID)

# MC only ASVs on riftia only
riftia_only_18s <- key_sub_dist_18 %>%
  filter(Distribution == "Riftia only")
length(riftia_only_18s$FeatureID)

# MC only ASVs on quartz only
quartz_only_18s <- key_sub_dist_18 %>%
  filter(Distribution == "Quartz only")
length(quartz_only_18s$FeatureID)

# MC only ASVs on all substrates 
onall_sub_18s <- key_sub_dist_18 %>%
  filter(Distribution == "All")
length(unique(onall_sub_18s$FeatureID))

### 16S
# Total ASVs - real total(contains repeats) & unique ASVs
length(df_16s_all$FeatureID)
length(unique(df_16s_all$FeatureID))

# ASVs on substrates (MCs)
ASVs_MC_16s <- df_16s_all %>%
  separate(SAMPLE, into = c("MC", "SUBSTRATE"), sep = "-", remove = FALSE) %>%
  filter(!(is.na(SUBSTRATE)))
length(ASVs_MC_16s$FeatureID)
length(unique(ASVs_MC_16s$FeatureID))

# ASVs on substrate only (MC only - no B/V)
df_16s_subs <- mc_tmp_16 %>%
  left_join(key_dist) %>%
  filter(DISTRIBUTION == "Substrate only")
length(df_16s_subs$FeatureID)
length(unique(df_16s_subs$FeatureID))

# MC only ASVs on more than one substrate
asvs_subs_only_16S <- df_16s_subs %>%
  separate(SAMPLE, into = c("MC", "SUBSTRATE"), sep = "-", remove = FALSE)

format_16s_subonly <- asvs_subs_only_16S %>%
  select(FeatureID, SEQ_AVG, SUBSTRATE) %>%
  pivot_wider(id_cols = FeatureID, names_from = SUBSTRATE, values_fn = mean, values_from = SEQ_AVG, values_fill = NA)

format_16s_subonly_assigned <- format_16s_subonly %>%
  mutate(Distribution = case_when(
    (is.na(Riftia) & is.na(Shell) & !(is.na(Quartz))) ~ "Quartz only", 
    (is.na(Riftia) & !(is.na(Shell)) & !(is.na(Quartz))) ~ "Quartz & Shell only", 
    (is.na(Riftia) & !(is.na(Shell)) & is.na(Quartz)) ~ "Shell only", 
    (!(is.na(Riftia)) & is.na(Shell) & !(is.na(Quartz))) ~ "Quartz & Riftia only",
    (!(is.na(Riftia)) & !(is.na(Shell)) & is.na(Quartz)) ~ "Shell & Riftia only", 
    (!(is.na(Riftia)) & is.na(Shell) & is.na(Quartz)) ~ "Riftia only", 
    (!(is.na(Riftia)) & !(is.na(Shell)) & !(is.na(Quartz)) ~ "All")
  ))

key_sub_dist_16s <- format_16s_subonly_assigned %>%
  select(FeatureID, Distribution) %>%
  distinct()
unique(key_sub_dist_16s$Distribution)

more_than_one_16s <- key_sub_dist_16s %>%
  filter(!(Distribution == "Quartz only" | Distribution == "Riftia only" | Distribution == "Shell only"))
length(more_than_one_16s$FeatureID)

# MC only ASVs on shell only
shell_only_16s <- key_sub_dist_16s %>%
  filter(Distribution == "Shell only")
length(shell_only_16s$FeatureID)

# MC only ASVs on riftia only
riftia_only_16s <- key_sub_dist_16s %>%
  filter(Distribution == "Riftia only")
length(riftia_only_16s$FeatureID)

# MC only ASVs on quartz only
quartz_only_16s <- key_sub_dist_16s %>%
  filter(Distribution == "Quartz only")
length(quartz_only_16s$FeatureID)

# MC only ASVs on all substrates 
onall_sub_16s <- key_sub_dist_16s %>%
  filter(Distribution == "All")
length(unique(onall_sub_16s$FeatureID))


#### Results: Microbe diversity-comparing Vent & Sub vs. Sub only ####
# 18S-stacked bar plot of community compostion (Relative ASVs)
svs_taxa_18S <- df_18s_svs %>%
  group_by(DISTRIBUTION, Phylum) %>%
  summarise(Count = n(), .groups = 'drop') %>%
  filter(!(is.na(Phylum)))

filt_df_18s_svs <- df_18s_svs %>%
  filter(SEQ_AVG > 10) %>%
  filter(!(is.na(Phylum)))

filt_svs_taxa_18S <- filt_df_18s_svs %>%
  group_by(DISTRIBUTION, Phylum) %>%
  summarise(Count = n(), .groups = 'drop')

ggplot(svs_taxa_18S, aes(x= DISTRIBUTION, y = Count, fill = Phylum)) +
  geom_bar(stat = "identity", position = "fill", color = "black") +
  theme_classic()
ggplot(filt_svs_taxa_18S, aes(x= DISTRIBUTION, y = Count, fill = Phylum)) +
  geom_bar(stat = "identity", position = "fill", color = "black") +
  theme_classic()

# 18S- Shannon diversity
wide_svs_18s <- df_18s_svs %>%
  group_by(DISTRIBUTION, FeatureID) %>%
  summarize(SUM = sum(SEQ_AVG)) %>%
  ungroup() %>%
  select(DISTRIBUTION, FeatureID, SUM) %>%
  pivot_wider(names_from = DISTRIBUTION, values_from = SUM, values_fill = 0) %>%
  column_to_rownames(var = "FeatureID")

wide_svs_18s_mat <- as.matrix(wide_svs_18s)
stand_wide_svs_18S_mat <- decostand(wide_svs_18s_mat, method = "hellinger", MARGIN = 2)
shannon_svs_18S <- diversity(wide_svs_18s_mat, index = "shannon", MARGIN = 2)
shannon_svs_18S
plot_shannon_svs_18S <- as.data.frame(shannon_svs_18S) %>%
  rownames_to_column(var = "DISTRIBUTION")
plot_shannon_svs_18S
ggplot(plot_shannon_svs_18S, aes(x= DISTRIBUTION, y= shannon_svs_16S)) +
  geom_point(color = "purple", size = 4, shape = 18) +
  labs(x = "", y= "18S Shannon diversity")

# 16S- stacked bar plot of community composition (Relative ASVs) 
df_16s_svs <- mc_tmp_16 %>%
  filter(!(grepl("Olivine", SAMPLE)) & 
           !(grepl("Basalt", SAMPLE)) &
           !(grepl("Pyrite", SAMPLE))) %>%
  left_join(key_dist) %>%
  filter(DISTRIBUTION == "Substrate only" | DISTRIBUTION == "Vent & Substrate only")

svs_taxa <- df_16s_svs %>%
  group_by(DISTRIBUTION, Phylum) %>%
  summarise(Count = n(), .groups = 'drop')

filt_df_16s_svs <- df_16s_svs %>%
  filter(SEQ_AVG > 10)

filt_svs_taxa <- filt_df_16s_svs %>%
  group_by(DISTRIBUTION, Phylum) %>%
  summarise(Count = n(), .groups = 'drop')

ggplot(svs_taxa, aes(x= DISTRIBUTION, y = Count, fill = Phylum)) +
  geom_bar(stat = "identity", position = "fill", color = "black") +
  theme_classic()
ggplot(filt_svs_taxa, aes(x= DISTRIBUTION, y = Count, fill = Phylum)) +
  geom_bar(stat = "identity", position = "fill", color = "black") +
  theme_classic()

# 16S- Shannon diversity sub only vs vent&sub
wide_svs_16s <- df_16s_svs %>%
  group_by(DISTRIBUTION, FeatureID) %>%
  summarize(SUM = sum(SEQ_AVG)) %>%
  ungroup() %>%
  select(DISTRIBUTION, FeatureID, SUM) %>%
  pivot_wider(names_from = DISTRIBUTION, values_from = SUM, values_fill = 0) %>%
  column_to_rownames(var = "FeatureID")

wide_svs_16s_mat <- as.matrix(wide_svs_16s)
stand_wide_svs_16S_mat <- decostand(wide_svs_16s_mat, method = "hellinger", MARGIN = 2)
shannon_svs_16S <- diversity(wide_svs_16s_mat, index = "shannon", MARGIN = 2)
shannon_svs_16S
plot_shannon_svs_16S <- as.data.frame(shannon_svs_16S) %>%
  rownames_to_column(var = "DISTRIBUTION")
plot_shannon_svs_16S
ggplot(plot_shannon_svs_16S, aes(x= DISTRIBUTION, y= shannon_svs_16S)) +
  geom_point(color = "purple", size = 4, shape = 18) +
  labs(x = "", y= "16S Shannon diversity")




#### IDing primary colonizers ####
#ASVs ONLY on MCs
MC_only_18S <- mc_df_18s %>%
  filter(DISTRIBUTION == "Substrate only")

MC_only_16S <- mc_df_16s %>%
  filter(DISTRIBUTION == "Substrate only")
MC_only_18S
MC_only_16S

plot_MC_only_18S <- MC_only_18S %>%
  group_by(SAMPLE, Phylum) %>%
  summarise(Count = n(), .groups = 'drop')
ggplot(plot_MC_only_18S, aes(x= SAMPLE, y = Count, fill = Phylum)) +
  geom_bar(stat = "identity", position = "fill", color = "black") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = -45))

#Of these, ASVs at more than one site
ASVs_18S_mult_sites <- MC_only_18S %>%
  group_by(FeatureID) %>%
  summarize(site_count = n_distinct(SAMPLE)) %>%
  filter(site_count > 1)

ASVs_18S_df <- MC_only_18S %>%
  left_join(ASVs_18S_mult_sites) %>%
  filter(!(is.na(site_count)))

plot_MC_ASVs_18S <- ASVs_18S_df %>%
  group_by(SAMPLE, Phylum) %>%
  summarise(Count = n(), .groups = 'drop') 
ggplot(plot_MC_ASVs_18S, aes(x= SAMPLE, y = Count, fill = Phylum)) +
  geom_bar(stat = "identity", position = "fill", color = "black") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = -45))

# 16S
ASVs_16S_mult_sites <- MC_only_16S %>%
  group_by(FeatureID) %>%
  summarize(site_count = n_distinct(SAMPLE)) %>%
  filter(site_count > 1)

ASVs_16S_df <- MC_only_16S %>%
  left_join(ASVs_16S_mult_sites) %>%
  filter(!(is.na(site_count)))

plot_MC_ASVs_16S <- ASVs_16S_df %>%
  group_by(SAMPLE, Phylum) %>%
  summarise(Count = n(), .groups = 'drop') 
ggplot(plot_MC_ASVs_16S, aes(x= SAMPLE, y = Count, fill = Phylum)) +
  geom_bar(stat = "identity", position = "fill", color = "black") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = -45))


### point plot of 18S ASVs on substrates ####
ASV_onsubs <- stacked_asvs %>%
  separate(SAMPLE, into = c("MC", "SUBSTRATE"), sep = "-", remove = FALSE) %>%
  filter(SUBSTRATE == "Quartz" | SUBSTRATE == "Riftia" | SUBSTRATE == "Shell") %>%
  filter(!(Phylum == "Opisthokonta_X" | Phylum == "Stramenopiles_X"))
ggplot(ASV_onsubs, aes(x = SUBSTRATE, y = Phylum)) +
  geom_point(aes(size = Count), color = "blue", alpha= 0.7) +
  facet_grid(rows = vars(Supergroup), scales = "free", space = "free") + 
  scale_size_continuous(name = "ASV Count", range = c(2, 10)) +
  theme_classic()


#### PCA Analysis ####
## 18S
subset_substrate_df <- asv_wtax_18 %>%
  filter(SAMPLETYPE == "Microcolonizer") #%>%

matrix_substrate_wide <- subset_substrate_df %>%
  select(SAMPLE, FeatureID, SEQUENCE_COUNT) %>%
  pivot_wider(id_cols = FeatureID, names_from = SAMPLE, values_from = SEQUENCE_COUNT, values_fill = 0) %>%
  column_to_rownames(var = "FeatureID") %>%
  as.matrix()
covariance_substrate <- as.matrix(matrix_substrate_wide) %*% t(matrix_substrate_wide)
det(covariance_substrate)

# PCA by SAMPLE MC-SUB
clr_substrate_sample <- compositions::clr(t(matrix_substrate_wide))
class(clr_substrate_sample)

# Transform data back to a df
clr_df_sample <- data.frame(clr_substrate_sample)
colnames(clr_df_sample)
head(clr_df)

covar_clr <- t(clr_df_sample) %*% as.matrix(clr_df_sample)
det(covar_clr)

clr_pca_sample <- prcomp(clr_df_sample)

# Make scree plot
clr_variances <- as.data.frame(clr_pca_sample$sdev^2/sum(clr_pca_sample$sdev^2)) %>% #Extract axes
  # Format to plot
  select(PercVar = 'clr_pca_sample$sdev^2/sum(clr_pca_sample$sdev^2)') %>% 
  rownames_to_column(var = "PCaxis") %>% # Pull out PC axis
  data.frame

ggplot(clr_variances, aes(x = as.numeric(PCaxis), y = PercVar)) + 
  geom_bar(stat = "identity") +
  theme_classic() +
  labs(x = "PC axis", y = "% Variance", title = "Log-Ratio PCA Screeplot")

# plot
pca_df_toplot <- data.frame(clr_pca_sample$x) %>%
  rownames_to_column(var = "SAMPLE") %>%
  separate(SAMPLE, c("number", "VentSite", "tmp", "Microcolonizer_id", "Substrate", "excess"), sep = "_") 

ggplot(pca_df_toplot, aes(x = PC1, y = PC2, fill = Microcolonizer_id, shape = Substrate)) +
  geom_point(color = "black", size = 3) +
  scale_shape_manual(values = c(21, 22, 24)) +
  theme_classic() +
  ggtitle('CLR PCA Ordination for 18s') +
  guides(fill = guide_legend(override.aes = list(shape = 22)))

# PCA for each SUBSTRATE + VENT
df_18s_prePCA <- asv_wtax_18 |> 
  select(-SAMPLE) |> 
  filter(SAMPLETYPE == "Microcolonizer" | VENT == "Mt Edwards" | VENT == "Deep seawater") |> 
  filter(Domain == "Eukaryota") |> 
  filter(!(Phylum == "NA")) %>%
  mutate(SUBSTRATE = str_replace_all(Substrate, "z2", "z")) |> 
  mutate(SAMPLE = case_when(
    SAMPLETYPE == "Microcolonizer" ~ paste(MC, SUBSTRATE, sep = "-"),
    VENT == "Deep seawater" ~ "Background seawater",
    VENT == "Mt Edwards" ~ "Mt Edwards diffuse fluid"
  )) |>
  # group_by(SAMPLE, FeatureID, Taxon, Domain, 
  #          Supergroup, Phylum, Class, Order, Family, Genus, Species, MC, SEQUENCE_COUNT) %>%
  mutate(TAXA_REVISED = case_when(
    Supergroup == "Alveolata" ~ paste(Supergroup, Phylum, sep = "-"),
    TRUE ~ Supergroup
  )) %>% 
  unite(TAX_LEVEL, Supergroup, Phylum, sep = "-", remove = FALSE) %>%  
  unite(TAX_MC_SUBSTRATE, TAXA_REVISED, MC, SUBSTRATE, sep = "_", remove = FALSE) 

df_18s_pre_PCA <- asv_wtax_18 |> 
  select(-SAMPLE) |> 
  filter(SAMPLETYPE == "Microcolonizer" | VENT == "Mt Edwards" | VENT == "Deep seawater") |> 
  filter(Domain == "Eukaryota") |> 
  filter(!(Phylum == "NA")) %>%
  mutate(SUBSTRATE = str_replace_all(Substrate, "z2", "z")) |> 
  mutate(SAMPLE = case_when(
    SAMPLETYPE == "Microcolonizer" ~ paste(MC, SUBSTRATE, sep = "-"),
    VENT == "Deep seawater" ~ "Background seawater",
    VENT == "Mt Edwards" ~ "Mt Edwards diffuse fluid"
  )) |>
  mutate(TAXA_REVISED = case_when(
    Supergroup == "Alveolata" ~ paste(Supergroup, Phylum, sep = "-"),
    TRUE ~ Supergroup
  )) %>% 
  unite(TAX_LEVEL, Supergroup, Phylum, sep = "-", remove = FALSE)

unique(df_18s_prePCA$TAX_LEVEL)
unique(df_18s_prePCA$SAMPLE)

# Make 3 separate matrices
matrix_shell_wide <- df_18s_prePCA %>%
  filter(SUBSTRATE == "Shell") %>% 
  ungroup() %>%
  select(FeatureID, TAX_MC_SUBSTRATE, SEQUENCE_COUNT) %>%
  mutate(SEQUENCE_COUNT = as.numeric(SEQUENCE_COUNT)) %>% 
  pivot_wider(id_cols = FeatureID, names_from = TAX_MC_SUBSTRATE, values_from = SEQUENCE_COUNT, values_fill = 0, values_fn = mean) %>%
  column_to_rownames(var = "FeatureID") %>%
  as.matrix()

## Repeat this with other substrates
matrix_riftia_wide <- df_18s_prePCA %>%
  filter(SUBSTRATE == "Riftia") %>% 
  ungroup() %>%
  select(FeatureID, TAX_MC_SUBSTRATE, SEQUENCE_COUNT) %>%
  mutate(SEQUENCE_COUNT = as.numeric(SEQUENCE_COUNT)) %>% 
  pivot_wider(id_cols = FeatureID, names_from = TAX_MC_SUBSTRATE, values_from = SEQUENCE_COUNT, values_fill = 0, values_fn = mean) %>%
  column_to_rownames(var = "FeatureID") %>%
  as.matrix()

matrix_quartz_wide <- df_18s_prePCA %>%
  filter(SUBSTRATE == "Quartz") %>% 
  ungroup() %>%
  select(FeatureID, TAX_MC_SUBSTRATE, SEQUENCE_COUNT) %>%
  mutate(SEQUENCE_COUNT = as.numeric(SEQUENCE_COUNT)) %>% 
  pivot_wider(id_cols = FeatureID, names_from = TAX_MC_SUBSTRATE, values_from = SEQUENCE_COUNT, values_fill = 0, values_fn = mean) %>%
  column_to_rownames(var = "FeatureID") %>%
  as.matrix()

# vent
matrix_vent_wide <- df_18s_pre_PCA %>%
  filter(SAMPLE == "Mt Edwards diffuse fluid") %>% 
  ungroup() %>%
  select(FeatureID, TAXA_REVISED, SEQUENCE_COUNT) %>%
  mutate(SEQUENCE_COUNT = as.numeric(SEQUENCE_COUNT)) %>% 
  pivot_wider(id_cols = FeatureID, names_from = TAXA_REVISED, values_from = SEQUENCE_COUNT, values_fill = 0, values_fn = mean) %>%
  column_to_rownames(var = "FeatureID") %>%
  as.matrix()

view(matrix_shell_wide)
class(matrix_shell_wide)


covariance_tax <- as.matrix(matrix_shell_wide) %*% t(matrix_shell_wide)
det(covariance_tax)

covariance_riftia <- as.matrix(matrix_riftia_wide) %*% t(matrix_riftia_wide)
det(covariance_riftia)

covariance_quartz <- as.matrix(matrix_quartz_wide) %*% t(matrix_quartz_wide)
det(covariance_quartz)

covariance_vent <- as.matrix(matrix_vent_wide) %*% t(matrix_vent_wide)
det(covariance_vent)

# plot by tax level
clr_tax <- compositions::clr(t(matrix_shell_wide))
class(clr_tax)

clr_riftia <- compositions::clr(t(matrix_riftia_wide))
class(clr_riftia)

clr_quartz <- compositions::clr(t(matrix_quartz_wide))
class(clr_quartz)

clr_vent <- compositions::clr(t(matrix_vent_wide))
class(clr_vent)

# Transform data class back to a Data Frame
clr_df_tax <- data.frame(clr_tax)
colnames(clr_df_tax)
covar_clr_tax <- t(clr_df_tax) %*% as.matrix(clr_df_tax)
det(covar_clr_tax)

clr_df_riftia <- data.frame(clr_riftia)
colnames(clr_df_riftia)
covar_clr_riftia <- t(clr_df_riftia) %*% as.matrix(clr_df_riftia)
det(covar_clr_riftia)

clr_df_quartz <- data.frame(clr_quartz)
colnames(clr_df_quartz)
covar_clr_quartz <- t(clr_df_quartz) %*% as.matrix(clr_df_quartz)
det(covar_clr_tax)

clr_df_vent <- data.frame(clr_vent)
colnames(clr_df_vent)
covar_clr_vent <- t(clr_df_vent) %*% as.matrix(clr_df_vent)
det(covar_clr_vent)

# Histogram check
# hist(clr_df_tax$X00b72d1a5fefb03bc39edbfba4566d03)

clr_pca_tax <- prcomp(clr_df_tax)
class(clr_pca_tax)
summary(clr_pca_tax)

clr_pca_riftia <- prcomp(clr_df_riftia)
class(clr_pca_riftia)
summary(clr_pca_riftia)

clr_pca_quartz <- prcomp(clr_df_quartz)
class(clr_pca_quartz)
summary(clr_pca_quartz)

clr_pca_vent <- prcomp(clr_df_vent)
class(clr_pca_vent)
summary(clr_pca_vent)

# Make scree plot
clr_var_tax <- as.data.frame(clr_pca_tax$sdev^2/sum(clr_pca_tax$sdev^2)) %>% 
  select(PercVar = 'clr_pca_tax$sdev^2/sum(clr_pca_tax$sdev^2)') %>% 
  rownames_to_column(var = "PCaxis") %>% 
  data.frame

clr_var_riftia <- as.data.frame(clr_pca_riftia$sdev^2/sum(clr_pca_riftia$sdev^2)) %>% 
  select(PercVar = 'clr_pca_riftia$sdev^2/sum(clr_pca_riftia$sdev^2)') %>% 
  rownames_to_column(var = "PCaxis") %>% 
  data.frame

clr_var_quartz <- as.data.frame(clr_pca_quartz$sdev^2/sum(clr_pca_quartz$sdev^2)) %>% 
  select(PercVar = 'clr_pca_quartz$sdev^2/sum(clr_pca_quartz$sdev^2)') %>% 
  rownames_to_column(var = "PCaxis") %>% 
  data.frame

clr_var_vent <- as.data.frame(clr_pca_vent$sdev^2/sum(clr_pca_vent$sdev^2)) %>% 
  select(PercVar = 'clr_pca_vent$sdev^2/sum(clr_pca_vent$sdev^2)') %>% 
  rownames_to_column(var = "PCaxis") %>% 
  data.frame

#ggplot Scree plot
ggplot(clr_var_tax, aes(x = as.numeric(PCaxis), y = PercVar)) + 
  geom_bar(stat = "identity") +
  theme_classic() +
  labs(x = "PC axis", y = "% Variance", title = "Log-Ratio PCA Screeplot")

ggplot(clr_var_riftia, aes(x = as.numeric(PCaxis), y = PercVar)) + 
  geom_bar(stat = "identity") +
  theme_classic() +
  labs(x = "PC axis", y = "% Variance", title = "Log-Ratio PCA Screeplot")

ggplot(clr_var_quartz, aes(x = as.numeric(PCaxis), y = PercVar)) + 
  geom_bar(stat = "identity") +
  theme_classic() +
  labs(x = "PC axis", y = "% Variance", title = "Log-Ratio PCA Screeplot")

ggplot(clr_var_vent, aes(x = as.numeric(PCaxis), y = PercVar)) + 
  geom_bar(stat = "identity") +
  theme_classic() +
  labs(x = "PC axis", y = "% Variance", title = "Log-Ratio PCA Screeplot")

head(clr_pca_tax$x)

## Extract df for PCA plotting, and separate our "sample"
pca_tax_toplot <- data.frame(clr_pca_tax$x) %>%
  rownames_to_column(var = "TMP") %>%
  separate(TMP, c("TAX", "MC", "Substrate"), sep = "_")

pca_riftia_toplot <- data.frame(clr_pca_riftia$x) %>%
  rownames_to_column(var = "TMP") %>%
  separate(TMP, c("TAX", "MC", "Substrate"), sep = "_")

pca_quartz_toplot <- data.frame(clr_pca_quartz$x) %>%
  rownames_to_column(var = "TMP") %>%
  separate(TMP, c("TAX", "MC", "Substrate"), sep = "_")

pca_vent_toplot <- data.frame(clr_pca_vent$x) %>%
  rownames_to_column(var = "Supergroup")
unique(pca_vent_toplot$Supergroup)
view(pca_tax_toplot)

ggplot(pca_tax_toplot, aes(x = PC1, y = PC2, fill = TAX, shape = MC)) +
  geom_point(size = 3, color = "black") +
  scale_shape_manual(values = c(21, 22, 23, 24, 25)) +
  theme_classic() +
  ggtitle('CLR PCA Ordination for 18s on Shell') +
  guides(fill = guide_legend(override.aes = list(shape = 22))) +
  theme(legend.position = "right")

ggplot(pca_riftia_toplot, aes(x = PC1, y = PC2, fill = TAX, shape = MC)) +
  geom_point(size = 3, color = "black") +
  scale_shape_manual(values = c(21, 22, 23, 24, 25)) +
  theme_classic() +
  ggtitle('CLR PCA Ordination for 18s on Riftia') +
  guides(fill = guide_legend(override.aes = list(shape = 22))) +
  theme(legend.position = "right")

ggplot(pca_quartz_toplot, aes(x = PC1, y = PC2, fill = TAX, shape = MC)) +
  geom_point(size = 3, color = "black") +
  scale_shape_manual(values = c(21, 22, 23, 24, 25)) +
  theme_classic() +
  ggtitle('CLR PCA Ordination for 18s on Quartz') +
  guides(fill = guide_legend(override.aes = list(shape = 22))) +
  theme(legend.position = "right")

ggplot(pca_vent_toplot, aes(x = PC1, y = PC2, fill = Supergroup, shape = Supergroup)) +
  geom_point(size = 3, color = "black", shape = 21) +
  theme_classic() +
  scale_fill_manual(values = c("#ff595e", "#ff924c", "#ffca3a", "#c5ca30",  "#8ac926",  "#52a675",  "#1982c4",  "#4267ac",  "#6a4c93",  "#b5a6c9")) +
  ggtitle('CLR PCA Ordination for 18s Vent fluid') +
  guides(fill = guide_legend(override.aes = list(shape = 22))) +
  theme(legend.position = "right")

## 16S
mc_tmp_16s <- asv16s_df |> 
  filter(Sampletype == "Microcolonizer" | LocationName == "Mt Edwards Vent" | LocationName == "Deep seawater" & STATUS == "keep") |> 
  mutate(SAMPLE = case_when(
    Sampletype == "Microcolonizer" ~ paste(MC, Substrate, sep = "-"),
    LocationName == "Deep seawater" ~ "Background seawater",
    LocationName == "Mt Edwards Vent" ~ "Mt Edwards diffuse fluid"
  )) |> 
  mutate(SAMPLE  = case_when(
    `SAMPLE` == "MC3-shell" ~ "MC3-Shell",
    .default = `SAMPLE`)) |>
  group_by(SAMPLE, FeatureID, Taxon, Domain, 
           Phylum, Class, Order, Family, Genus, Species) |> 
  filter(!(Domain == "Unassigned") & !(Domain == "Eukaryota")) |> 
  filter(SEQUENCE_COUNT > 0)

df_16s_to_use <- mc_tmp_16s %>%
  filter(!(grepl("Olivine", SAMPLE)) & 
           !(grepl("Basalt", SAMPLE)) &
           !(grepl("Pyrite", SAMPLE))) %>%
  mutate(SAMPLE_TYPE = case_when(
    SAMPLE == "Background seawater" ~ "Background", 
    SAMPLE == "Mt Edwards diffuse fluid" ~ "Vent", 
    TRUE ~ "Substrate"
  ))

subset_sub_16s <- df_16s_to_use %>%
  filter(Sampletype == "Microcolonizer") 

matrix_sub_wide <- subset_sub_16s %>%
  ungroup() %>%  
  select(SAMPLE, FeatureID, SEQUENCE_COUNT) %>% 
  mutate(SEQUENCE_COUNT = as.numeric(SEQUENCE_COUNT)) %>% 
  pivot_wider(id_cols = FeatureID, names_from = SAMPLE, values_from = SEQUENCE_COUNT, values_fill = 0, values_fn = mean) %>% 
  column_to_rownames(var = "FeatureID") %>% 
  as.matrix()

covariance_sub_16s <- as.matrix(matrix_sub_wide) %*%
  t(matrix_sub_wide)
det(covariance_sub_16s)

# data transformation # 
#plot by samples rows=sample cols=FeatureID
clr_subs_16s_sample <- compositions::clr(t(matrix_sub_wide))

clr_df_16s_samp <- data.frame(clr_subs_16s_sample)
colnames(clr_df_16s_samp)

covar_clr_16s_samp <- t(clr_df_16s_samp) %*% as.matrix(clr_df_16s_samp)
det(covar_clr_16s_samp)

# hist check
hist(clr_df_16s_samp$X00253cc12b812cda98d681c0a72a4b20)

#PCA
clr_16s_pca_2 <- prcomp(clr_df_16s_samp)
class(clr_16s_pca_2)
summary(clr_16s_pca_2)

# Scree
clr_variances2 <- as.data.frame(clr_16s_pca_2$sdev^2/sum(clr_16s_pca_2$sdev^2)) %>%
  select(PercVar = 'clr_16s_pca_2$sdev^2/sum(clr_16s_pca_2$sdev^2)') %>%
  rownames_to_column(var = "PCaxis") %>% 
  data.frame

ggplot(clr_variances2, aes(x = as.numeric(PCaxis), y = PercVar)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  labs(x = "PC axis", y = "Percent Variance")

df_16s_pca2_p <- data.frame(clr_16s_pca_2$x) %>%
  rownames_to_column(var = "SAMPLE") %>% 
  separate(SAMPLE, c("Microcolonizer", "Substrate"), sep = "-")

ggplot(df_16s_pca2_p, aes(x = PC1, y = PC2, fill = Microcolonizer, shape = Substrate)) +
  geom_point(color = "black", size = 3) +
  scale_shape_manual(values = c(21, 22, 24)) +
  theme_classic() +
  ggtitle('CLR PCA Ordination for 16s') +
  guides(fill = guide_legend(override.aes = list(shape = 22)))


#### ####