# Species occurance check list

# Read in libraries
library(tidyverse)

# Read in data
CoCos_false  <- readRDS('data/phyloseq_objects/CoCosV10I_16S_phyloseq_nt_FALSE_decontaminated.rds')
CoCos_true  <- readRDS('data/phyloseq_objects/CoCosV10I_16S_phyloseq_nt_TRUE_decontaminated.rds')
CoCos_pseudo <- readRDS('data/phyloseq_objects/CoCosV10I_16S_phyloseq_nt_pseudo_decontaminated.rds')

tax_false_cocos <- as.data.frame(CoCos_false@tax_table)
tax_pseudo_cocos <- as.data.frame(CoCos_pseudo@tax_table)
tax_true_cocos <- as.data.frame(CoCos_true@tax_table)

RS_false  <- readRDS('data/phyloseq_objects/RS21AUG_16S_phyloseq_nt_FALSE_decontam.rds')
RS_true  <- readRDS('data/phyloseq_objects/RS21AUG_16S_phyloseq_nt_TRUE_decontam.rds')
RS_pseudo <- readRDS('data/phyloseq_objects/RS21AUG_16S_phyloseq_nt_pseudo_decontam.rds')

tax_false_RS <- as.data.frame(RS_false@tax_table)
tax_pseudo_RS <- as.data.frame(RS_pseudo@tax_table)
tax_true_RS <- as.data.frame(RS_true@tax_table)

NW_false  <- readRDS('data/phyloseq_objects/Pool_FALSE_KWEST_16S_phyloseq_nt_decontaminated.rds')
NW_true  <- readRDS('data/phyloseq_objects/Pool_TRUE_KWEST_16S_phyloseq_nt_decontaminated.rds')
NW_pseudo <- readRDS('data/phyloseq_objects/Pool_pseudo_KWEST_16S_phyloseq_nt_decontaminated.rds')

tax_false_NW <- as.data.frame(NW_false@tax_table)
tax_pseudo_NW <- as.data.frame(NW_pseudo@tax_table)
tax_true_NW <- as.data.frame(NW_true@tax_table)

checklist <- read_csv(file = "data/Aust_fish_species_list.csv")
head(checklist)


## create list of species

### Cocos dataset
tax_false_cocos %>%
  filter(!(species %in% tax_false_cocos$species[grep(NA, (tax_false_cocos$species))])) %>%
  filter(!(species %in% tax_false_cocos$species[grep("dropped", (tax_false_cocos$species))])) %>%
  distinct(species) %>%
  mutate(Dataset = 'Cocos Islands',
         Mode = 'Independent') -> cocos_false_species

tax_true_cocos %>%
  filter(!(species %in% tax_true_cocos$species[grep(NA, (tax_true_cocos$species))])) %>%
  filter(!(species %in% tax_true_cocos$species[grep("dropped", (tax_true_cocos$species))])) %>%
  distinct(species) %>%
  mutate(Dataset = 'Cocos Islands',
         Mode = 'Pooled') -> cocos_true_species

tax_pseudo_cocos %>%
  filter(!(species %in% tax_pseudo_cocos$species[grep(NA, (tax_pseudo_cocos$species))])) %>%
  filter(!(species %in% tax_pseudo_cocos$species[grep("dropped", (tax_pseudo_cocos$species))])) %>%
  distinct(species) %>%
  mutate(Dataset = 'Cocos Islands',
         Mode = 'Pseudo') -> cocos_pseudo_species


### Rowley Shoals dataset
tax_false_RS %>%
  # select(species) %>% 
  filter(!(species %in% tax_false_RS$species[grep(NA, (tax_false_RS$species))])) %>%
  filter(!(species %in% tax_false_RS$species[grep("dropped", (tax_false_RS$species))])) %>%
  distinct(species) %>%
  mutate(Dataset = 'Rowley Shoals',
         Mode = 'Independent') -> RS_false_species

tax_true_RS %>%
  # select(species) %>% 
  distinct(species) %>%
  filter(!(species %in% tax_true_RS$species[grep(NA, (tax_true_RS$species))])) %>%
  filter(!(species %in% tax_true_RS$species[grep("dropped", (tax_true_RS$species))])) %>%
  mutate(Dataset = 'Rowley Shoals',
         Mode = 'Pooled') -> RS_true_species

tax_pseudo_RS %>%
  filter(!(species %in% tax_pseudo_RS$species[grep(NA, (tax_pseudo_RS$species))])) %>%
  filter(!(species %in% tax_pseudo_RS$species[grep("dropped", (tax_pseudo_RS$species))])) %>%
  # select(species) %>% 
  distinct(species) %>%
  mutate(Dataset = 'Rowley Shoals',
         Mode = 'Pseudo') -> RS_pseudo_species



### North-west dataset
tax_false_NW %>%
  # select(species) %>% 
  distinct(species) %>%
  filter(!(species %in% tax_false_NW$species[grep(NA, (tax_false_NW$species))])) %>%
  filter(!(species %in% tax_false_NW$species[grep("dropped", (tax_false_NW$species))])) %>%
  mutate(Dataset = 'North-west',
         Mode = 'Independent') -> NW_false_species

tax_true_NW %>%
  # select(species) %>% 
  distinct(species) %>%
  filter(!(species %in% tax_true_NW$species[grep(NA, (tax_true_NW$species))])) %>%
  filter(!(species %in% tax_true_NW$species[grep("dropped", (tax_true_NW$species))])) %>%
  mutate(Dataset = 'North-west',
         Mode = 'Pooled') -> NW_true_species

tax_pseudo_NW %>%
  filter(!(species %in% tax_pseudo_NW$species[grep(NA, (tax_pseudo_NW$species))])) %>%
  filter(!(species %in% tax_pseudo_NW$species[grep("dropped", (tax_pseudo_NW$species))])) %>%
  # select(species) %>% 
  distinct(species) %>%
  mutate(Dataset = 'North-west',
         Mode = 'Pseudo') -> NW_pseudo_species


observed_species_list <- as_tibble(rbind(cocos_false_species,
                               cocos_true_species,
                               cocos_pseudo_species,
                               RS_false_species,
                               RS_true_species,
                               RS_pseudo_species,
                               NW_false_species,
                               NW_true_species,
                               NW_pseudo_species))

## intersection of checklist and observed species
intersect(observed_species_list$species, checklist$Species)
observed_species_list$checklist <- observed_species_list$species %in% checklist$Species

write.csv(observed_species_list, file = "Observed_species_list.csv")

length(observed_species_list$species)
length(unique(observed_species_list$species))

observed_species_list %>% 
  distinct(species, checklist) %>%
  filter(checklist == FALSE)  -> species_false
length(species_false$species)

observed_species_list %>% 
  distinct(species, checklist, Mode) %>%
  filter(checklist == FALSE) %>%
  filter(Mode == "Pooled")  -> species_false_pooled
length(species_false_pooled$species)

observed_species_list %>% 
  distinct(species, checklist) %>%
  filter(checklist == TRUE)  -> species_true
length(species_true$species)


### Investigate which species were found different between modes and if they are in the checklist
## Cocos
diff <- setdiff(tax_true_cocos$species, tax_false_cocos$species)
# [1] "Bolinichthys distofax"      "Notoscopelus caudispinosus" "Kali indica"               
# [4] "Lethrinus genivittatus"     "Vinciguerria nimbaria"      "Desmodema polystictum"     
# [7] "Trachinocephalus myops"     "Diretmus argenteus"         "Magnisudis atlantica"      
# [10] "Upeneichthys stotti"        "Kyphosus bigibbus" 
cocos_diff_false <- subset(observed_species_list, species %in% diff) %>%
  filter(checklist == FALSE) %>%
  filter(Dataset == "Cocos Islands") %>%
  distinct(species, Mode, checklist)
# 1 Bolinichthys distofax  Pooled FALSE
# 2 Kali indica            Pooled FALSE
# 3 Trachinocephalus myops Pooled FALSE    
# 4 Bolinichthys distofax  Pseudo FALSE    
# 5 Trachinocephalus myops Pseudo FALSE 

cocos_diff_true <- subset(observed_species_list, species %in% diff) %>%
  filter(checklist == TRUE) %>%
  filter(Dataset == "Cocos Islands") %>%
  distinct(species, Mode, checklist) 
length(unique(cocos_diff_true$species))
length(unique(cocos_diff_false$species))

diff <- setdiff(tax_false_cocos$species, tax_true_cocos$species)
# 0
cocos_diff <- subset(observed_species_list, species %in% diff) %>% 
  filter(checklist == FALSE) %>%
  distinct(species, Mode, checklist) ## 0

## Rowley Shoals
diff <- setdiff(tax_true_RS$species, tax_false_RS$species)

rs_diff_false <- subset(observed_species_list, species %in% diff) %>% 
  filter(checklist == FALSE) %>%
  filter(Dataset == "Rowley Shoals") %>%
  distinct(species, Mode, checklist)
# 1 Exocoetus environmental sample Pooled FALSE    
# 2 Pseudojuloides cerasinus       Pooled FALSE 
# 3 Amblyeleotris sp. SOK-2013     Pooled FALSE    
# 4 cf. Coris sp. AWFS-F15-1952    Pooled FALSE    
# 7 cf. Coris sp. AWFS-F15-1952    Pseudo FALSE

rs_diff_true <- subset(observed_species_list, species %in% diff) %>% 
  filter(checklist == TRUE) %>%
  filter(Dataset == "Rowley Shoals") %>%
  distinct(species, Mode, checklist)
length(unique(rs_diff_false$species))
length(unique(rs_diff_true$species))

diff <- setdiff(tax_false_RS$species, tax_true_RS$species)
# 0
rs_diff <- subset(observed_species_list, species %in% diff) %>% 
  filter(checklist == FALSE) %>%
  distinct(species, Mode, checklist) ## 0

## North-west
diff <- setdiff(tax_true_NW$species, tax_false_NW$species)

nw_diff_false <- subset(observed_species_list, species %in% diff) %>% 
  filter(checklist == FALSE) %>%
  filter(Dataset == "North-west") %>%
  distinct(species, Mode, checklist)
# 1 Ophicthus sp. Fj-1                       Pooled FALSE    
# 2 Chelonodontops patoca                    Pooled FALSE    
# 3 Epinephelus quoyans                      Pooled FALSE    
# 4 cf. Yongeichthys sp. AWFS-F16-1517       Pooled FALSE    
# 6 Pomadasys argyreus                       Pooled FALSE    
# 7 Plicofollis aff. argyropleuron 2 RB-2009 Pooled FALSE    
# 8 Neotrygon kuhlii                         Pooled FALSE    
# 9 Arripis georgianus                       Pooled FALSE    
# 10 Chelon planiceps                         Pooled FALSE    
# 11 Gazza sp. Fiji                           Pooled FALSE    
# 12 Aetobatus narinari                       Pooled FALSE    
# 13 cf. Polydactylus sp. AWFS-F16-1491       Pooled FALSE    
# 14 Aetobatus narinari                       Pseudo FALSE    
# 15 cf. Polydactylus sp. AWFS-F16-1491       Pseudo FALSE    

nw_diff_true <- subset(observed_species_list, species %in% diff) %>% 
  filter(checklist == TRUE) %>%
  filter(Dataset == "North-west") %>%
  distinct(species, Mode, checklist)
length(unique(nw_diff_false$species))
length(unique(nw_diff_true$species))

diff <- setdiff(tax_false_NW$species, tax_true_NW$species)
# [1] "Trachinocephalus myops"         "Decapterus akaadsi"             "Platax batavianus"              "Paramonacanthus choirocephalus"
# [5] "Torquigener whitleyi"           "Eleutheronema tetradactylum"    "Pomacentrus milleri"            "Ulua mentalis"                 
# [9] "Pseudorhombus jenynsii"
nw_diff <- subset(observed_species_list, species %in% diff) %>% 
  # filter(checklist == FALSE) %>%
  filter(Dataset == "North-west") %>%
  distinct(species, Mode, checklist)
# 1 Trachinocephalus myops         Independent FALSE    
# 2 Decapterus akaadsi             Independent FALSE    
# 3 Platax batavianus              Independent TRUE     
# 4 Paramonacanthus choirocephalus Independent TRUE     
# 5 Torquigener whitleyi           Independent TRUE     
# 6 Eleutheronema tetradactylum    Independent TRUE     
# 7 Pomacentrus milleri            Independent TRUE     
# 8 Ulua mentalis                  Independent TRUE     
# 9 Pseudorhombus jenynsii         Independent TRUE     
# 10 Platax batavianus              Pseudo      TRUE     
# 11 Paramonacanthus choirocephalus Pseudo      TRUE     
# 12 Torquigener whitleyi           Pseudo      TRUE     
# 13 Ulua mentalis                  Pseudo      TRUE     
# 14 Pseudorhombus jenynsii         Pseudo      TRUE

true_hits <- rbind(cocos_diff_true, rs_diff_true, nw_diff_true) %>%
  distinct(species)
length(unique(true_hits$species))

false_hits <- rbind(cocos_diff_false, rs_diff_false, nw_diff_false) %>%
  distinct(species)
length(unique(false_hits$species))
