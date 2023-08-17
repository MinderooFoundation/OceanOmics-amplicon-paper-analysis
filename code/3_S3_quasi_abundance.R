# Quasi abundance figure for paper

library(tidyverse)
library(tidyr)
library(phyloseq)
library(cowplot)

#==============================================================================
# Functions

add_zero_one <- function(x) {
  ifelse(x > 0, 1, 0)
}

# define a function to split a string by the last underscore and return the first element
split_string <- function(x) {
  split_x <- strsplit(x, "_(?!.*_)", perl=TRUE)
  return(split_x[[1]][1])
}

#==============================================================================
# CococsI
#==============================================================================
Cocos_independent  <- readRDS('data/phyloseq_objects/CoCosV10I_16S_phyloseq_nt_FALSE_decontaminated.rds')
Cocos_pooled  <- readRDS('data/phyloseq_objects/CoCosV10I_16S_phyloseq_nt_TRUE_decontaminated.rds')
Cocos_pseudo <- readRDS('data/phyloseq_objects/CoCosV10I_16S_phyloseq_nt_pseudo_decontaminated.rds')


format_quasi_abundance <- function(phyloseq){
  # Extract data tables 
  otu <- rownames_to_column(as.data.frame(otu_table(phyloseq)@.Data, var = "ASV")) %>% 
    rename(ASV = rowname) %>%
    select(-(contains(c('WC', 'BC', 'EB', 'DI'))))
  tax <- rownames_to_column(as.data.frame(tax_table(phyloseq)@.Data, var = "ASV")) %>% 
    rename(ASV = rowname)
  
  
  # apply the function to each string in the vector using lapply
  sample_names    <- names(otu)
  sample_names <- sample_names[-grep('ASV', sample_names)]
  sample_names <- unique(unlist(lapply(sample_names, split_string)))
  
  for(i in sample_names){
    
    otu_sample <- otu[,grep(paste0(i,'_'), names(otu))]
    
    # loop over the columns and apply the function
    df_new <- data.frame(lapply(otu_sample, add_zero_one)) 
    otu[,i] <- rowSums(df_new)
  }
  
  otu %>%
    dplyr::select(c(ASV, sample_names)) -> otu_sums
  
  merge(otu_sums, tax, by="ASV")
  
  
  # Reformat to long format
  otu_sums %>% 
    pivot_longer(sample_names) %>%
    rename(site_id = name,
           Abundance = value) -> otu_long
  
  merge(otu_long, tax, by="ASV")
  
}


pooled     <- format_quasi_abundance(Cocos_pooled) %>% mutate(Option = 'Pooled')
independent <- format_quasi_abundance(Cocos_independent) %>% mutate(Option = 'Independent')
pseudo     <- format_quasi_abundance(Cocos_pseudo) %>% mutate(Option = 'Pseudo')

cocos_species <- rbind(pooled, independent, pseudo)
head(cocos_species)
length(cocos_species$species)
length(unique(cocos_species$species))

cocos_species_ind <- cocos_species %>% 
  distinct(species, Abundance, Option) %>%
  filter(Option == "Independent")
length(cocos_species_ind$Abundance)

cocos_species_ind_1 <- cocos_species_ind %>% filter(Abundance <= 1)
length(cocos_species_ind_1$Abundance) # 130 = 48% of ASVs only found in 1 or 0 reps
cocos_species_ind_2 <- cocos_species_ind %>% filter(Abundance >= 2)
length(cocos_species_ind_2$Abundance) # 144
cocos_species_ind_3 <- cocos_species_ind %>% filter(Abundance >= 3)
length(cocos_species_ind_3$Abundance) # 95 = 35% of asvs found in 3 or more reps
cocos_species_ind_5 <- cocos_species_ind %>% filter(Abundance == 5)
length(cocos_species_ind_5$Abundance) # 19 = 7% of asvs found in 5 reps

cocos_species_pooled <- cocos_species %>% 
  distinct(species, Abundance, Option) %>%
  filter(Option == "Pooled")
length(cocos_species_pooled$Abundance)

cocos_species_pooled_1 <- cocos_species_pooled %>% filter(Abundance <= 1)
length(cocos_species_pooled_1$Abundance) # 154 = 41% of ASVs only found in 1 or 0 reps
cocos_species_pooled_2 <- cocos_species_pooled %>% filter(Abundance >= 2)
length(cocos_species_pooled_2$Abundance) # 219
cocos_species_pooled_3 <- cocos_species_pooled %>% filter(Abundance >= 3)
length(cocos_species_pooled_3$Abundance) # 153 = 41% of asvs found in 3 or more reps
cocos_species_pooled_5 <- cocos_species_pooled %>% filter(Abundance == 5)
length(cocos_species_pooled_5$Abundance) # 38 = 10% of asvs found in 5 reps

cocos_species <- rbind(pooled, independent, pseudo) %>%
  ggplot(aes(x=site_id, y=species, fill = Abundance)) +
  geom_tile() +
  theme_classic(base_size = 12) +
  facet_wrap(~Option) +
  theme(axis.text.x = element_blank()) +
  theme(axis.text.y = element_blank()) +
  scale_fill_gradient(low = "#8FBCBB", high = "#B48EAD")

Cocos_family <- 
  rbind(pooled, independent) %>%
  group_by(site_id,family, Option) %>% 
  summarise(Abundance_family = mean(Abundance)) %>%
  ggplot(aes(x=site_id, y=family, fill = Abundance_family)) +
  xlab(NULL) + 
  ylab(NULL) +
  geom_tile() +
  theme_classic(base_size = 12) +
  facet_wrap(~Option) +
  theme(axis.text.x = element_blank()) +
  theme(axis.text.y = element_blank()) +
  scale_fill_gradient(low = "#2E3440", high = "#B48EAD")


#==============================================================================
# North west WA data
#==============================================================================

NW_pooled     <- readRDS('data/phyloseq_objects/Pool_TRUE_KWEST_16S_phyloseq_nt_decontaminated.rds')
NW_independent <- readRDS('data/phyloseq_objects/Pool_FALSE_KWEST_16S_phyloseq_nt_decontaminated.rds')
NW_pseudo     <- readRDS('data/phyloseq_objects/Pool_pseudo_KWEST_16S_phyloseq_nt_decontaminated.rds')

format_quasi_abundance <- function(phyloseq){
  # Extract data tables 
  otu <- rownames_to_column(as.data.frame(otu_table(phyloseq)@.Data, var = "ASV")) %>% 
    rename(ASV = rowname) %>%
    select(-(contains(c('WC', 'BC', 'Ext', 'FC'))))
  tax <- rownames_to_column(as.data.frame(tax_table(phyloseq)@.Data, var = "ASV")) %>% rename(ASV = rowname)
  
  
  # apply the function to each string in the vector using lapply
  sample_names    <- names(otu)
  sample_names <- sample_names[-grep('ASV', sample_names)]
  sample_names <- unique(gsub("[^0-9]","", sample_names)) #unique(unlist(lapply(sample_names, split_string)))
  
  for(i in sample_names){
    
    otu_sample <- otu[,grep(paste0('^[',i,']','[a-z]'), names(otu))]
    
    # loop over the columns and apply the function
    df_new <- data.frame(lapply(otu_sample, add_zero_one)) 
    otu[,i] <- rowSums(df_new)
  }
  
  otu %>%
    dplyr::select(c(ASV, sample_names)) -> otu_sums
  
  merge(otu_sums, tax, by="ASV")
  
  
  # Reformat to long format
  otu_sums %>% 
    pivot_longer(sample_names) %>%
    rename(Site = name,
           Abundance = value) -> otu_long
  
  merge(otu_long, tax, by="ASV")
  
}


pooled     <- format_quasi_abundance(NW_pooled) %>% mutate(Option = 'Pooled')
independent <- format_quasi_abundance(NW_independent) %>% mutate(Option = 'Independent')
pseudo     <- format_quasi_abundance(NW_pseudo) %>% mutate(Option = 'Pseudo')

nw_species <- rbind(pooled, independent, pseudo)
nw_species_ind <- nw_species %>% 
  distinct(species, Abundance, Option) %>%
  filter(Option == "Independent")
length(nw_species_ind$Abundance)

nw_species_ind_1 <- nw_species_ind %>% filter(Abundance <= 1)
length(nw_species_ind_1$Abundance) # 115 = 72% of ASVs only found in 1 or 0 reps
nw_species_ind_2 <- nw_species_ind %>% filter(Abundance >= 2)
length(nw_species_ind_2$Abundance) # 45
nw_species_ind_3 <- nw_species_ind %>% filter(Abundance >= 3)
length(nw_species_ind_3$Abundance) # 26 = 16% of asvs found in 3 or more reps
nw_species_ind_5 <- nw_species_ind %>% filter(Abundance == 5)
length(nw_species_ind_5$Abundance) # 4 = 3% of asvs found in 5 reps

nw_species_pooled <- nw_species %>% 
  distinct(species, Abundance, Option) %>%
  filter(Option == "Pooled")
length(nw_species_pooled$Abundance)

nw_species_pooled_1 <- nw_species_pooled %>% filter(Abundance <= 1)
length(nw_species_pooled_1$Abundance) # 236 = 64% of ASVs only found in 1 or 0 reps
nw_species_pooled_2 <- nw_species_pooled %>% filter(Abundance >= 2)
length(nw_species_pooled_2$Abundance) # 131
nw_species_pooled_3 <- nw_species_pooled %>% filter(Abundance >= 3)
length(nw_species_pooled_3$Abundance) # 76 = 21% of asvs found in 3 or more reps
nw_species_pooled_5 <- nw_species_pooled %>% filter(Abundance == 5)
length(nw_species_pooled_5$Abundance) # 12 = 10% of asvs found in 5 reps

NW_species <- rbind(pooled, independent, pseudo) %>%
  ggplot(aes(x=Site, y=species, fill = Abundance)) +
  geom_tile() +
  theme_classic(base_size = 12) +
  facet_wrap(~Option) +
  theme(axis.text.x = element_blank()) +
  theme(axis.text.y = element_blank()) +
  scale_fill_gradient(low = "#8FBCBB", high = "#B48EAD")

NW_family <- 
  rbind(pooled, independent) %>%
  group_by(Site,family, Option) %>% 
  summarise(Abundance_family = mean(Abundance)) %>%
  ggplot(aes(x=Site, y=family, fill = Abundance_family)) +
  geom_tile() +
  xlab(NULL) + 
  ylab(NULL) +
  theme_classic(base_size = 12) +
  facet_wrap(~Option) +  
  theme(axis.text.x = element_blank()) +
  theme(axis.text.y = element_blank()) +
  scale_fill_gradient(low = "#2E3440", high = "#B48EAD")



#==============================================================================
# Rowley Shoals data
#==============================================================================

RS_pooled     <- readRDS('data/phyloseq_objects/RS21AUG_16S_phyloseq_nt_TRUE_decontam.rds')
RS_independent <- readRDS('data/phyloseq_objects/RS21AUG_16S_phyloseq_nt_FALSE_decontam.rds')
RS_pseudo     <- readRDS('data/phyloseq_objects/RS21AUG_16S_phyloseq_nt_pseudo_decontam.rds')

RS_pooled = subset_samples(RS_pooled, sample_names(RS_pooled) != "RS1_ME_S4_1_2_")
RS_independent = subset_samples(RS_independent, sample_names(RS_independent) != "RS1_ME_S4_1_2_")
RS_pseudo = subset_samples(RS_pseudo, sample_names(RS_pseudo) != "RS1_ME_S4_1_2_")

format_quasi_abundance <- function(phyloseq){
  # Extract data tables 
  otu <- rownames_to_column(as.data.frame(otu_table(phyloseq)@.Data, var = "ASV")) %>% 
    rename(ASV = rowname) %>%
    select(-(contains(c('WC', 'BC'))))
  tax <- rownames_to_column(as.data.frame(tax_table(phyloseq)@.Data, var = "ASV")) %>% rename(ASV = rowname)
  
  
  # apply the function to each string in the vector using lapply
  sample_names    <- names(otu)
  sample_names <- sample_names[-grep('ASV', sample_names)]
  sample_names <- unique(unlist(lapply(sample_names, split_string)))
  
  for(i in sample_names){
    
    otu_sample <- otu[,grep(i, names(otu))]
    
    # loop over the columns and apply the function
    df_new <- data.frame(lapply(otu_sample, add_zero_one)) 
    otu[,i] <- rowSums(df_new)
  }
  
  otu %>%
    dplyr::select(c(ASV, sample_names)) -> otu_sums
  
  merge(otu_sums, tax, by="ASV")
  
  
  # Reformat to long format
  otu_sums %>% 
    pivot_longer(sample_names) %>%
    rename(Site = name,
           Abundance = value) -> otu_long
  
  merge(otu_long, tax, by="ASV")
  
}


pooled     <- format_quasi_abundance(RS_pooled) %>% mutate(Option = 'Pooled')
independent <- format_quasi_abundance(RS_independent) %>% mutate(Option = 'Independent')
pseudo     <- format_quasi_abundance(RS_pseudo) %>% mutate(Option = 'Pseudo')

rs_species <- rbind(pooled, independent, pseudo)
rs_species_ind <- rs_species %>% 
  distinct(species, Abundance, Option) %>%
  filter(Option == "Independent")
length(rs_species_ind$Abundance)

rs_species_ind_1 <- rs_species_ind %>% filter(Abundance <= 1)
length(rs_species_ind_1$Abundance) # 329 = 49% of ASVs only found in 1 or 0 reps
rs_species_ind_2 <- rs_species_ind %>% filter(Abundance >= 2)
length(rs_species_ind_2$Abundance) # 344
rs_species_ind_3 <- rs_species_ind %>% filter(Abundance >= 3)
length(rs_species_ind_3$Abundance) # 209 = 31% of asvs found in 3 or more reps
rs_species_ind_5 <- rs_species_ind %>% filter(Abundance == 5)
length(rs_species_ind_5$Abundance) # 46 = 7% of asvs found in 5 reps

rs_species_pooled <- rs_species %>% 
  distinct(species, Abundance, Option) %>%
  filter(Option == "Pooled")
length(rs_species_pooled$Abundance)

rs_species_pooled_1 <- rs_species_pooled %>% filter(Abundance <= 1)
length(rs_species_pooled_1$Abundance) # 536 = 44% of ASVs only found in 1 or 0 reps
rs_species_pooled_2 <- rs_species_pooled %>% filter(Abundance >= 2)
length(rs_species_pooled_2$Abundance) # 690
rs_species_pooled_3 <- rs_species_pooled %>% filter(Abundance >= 3)
length(rs_species_pooled_3$Abundance) # 452 = 37% of asvs found in 3 or more reps
rs_species_pooled_5 <- rs_species_pooled %>% filter(Abundance == 5)
length(rs_species_pooled_5$Abundance) # 113 = 9% of asvs found in 5 reps

RS_species <- rbind(pooled, independent, pseudo) %>%
  ggplot(aes(x=Site, y=species, fill = Abundance)) +
  geom_tile() +
  theme_classic(base_size = 12) +
  facet_wrap(~Option) +
  theme(axis.text.x = element_blank()) +
  theme(axis.text.y = element_blank()) +
  scale_fill_gradient(low = "#8FBCBB", high = "#B48EAD")

RS_family <- 
  rbind(pooled, independent) %>%
  group_by(Site,family, Option) %>% 
  summarise(Abundance_family = mean(Abundance)) %>%
  ggplot(aes(x=Site, y=family, fill = Abundance_family)) +
  geom_tile() +
  xlab(NULL) + 
  ylab(NULL) +
  theme_classic(base_size = 12) +
  facet_wrap(~Option) +
  theme(axis.text.x = element_blank()) +
  theme(axis.text.y = element_blank()) +
  scale_fill_gradient(low = "#2E3440", high = "#B48EAD")

RS_family + labs(fill = "Replicates")


# Create plot grid
plot_grid(Cocos_family + theme(legend.position="bottom") + labs(fill = "Replicates") + ylab("Families") + xlab("Sites") + theme(axis.ticks = element_blank()),
                  RS_family + theme(legend.position="bottom") + labs(fill = "Replicates") + xlab("Sites") + theme(axis.ticks = element_blank()),
                  NW_family + theme(legend.position="bottom") + labs(fill = "Replicates") + xlab("Sites") + theme(axis.ticks = element_blank()),
                  labels = 'AUTO', nrow = 1, label_size = 12, hjust = -1.5, vjust = 2)


plot_grid(cocos_species + theme(legend.position="right") + labs(fill = "Replicates") + xlab("Sites"),
          RS_species + theme(legend.position="right") + labs(fill = "Replicates") + xlab("Sites"),
          NW_species + theme(legend.position="right") + labs(fill = "Replicates") + xlab("Sites"),
          labels = 'AUTO', nrow = 3, ncol = 1, label_size = 12, hjust = -1.5, vjust = 2)
