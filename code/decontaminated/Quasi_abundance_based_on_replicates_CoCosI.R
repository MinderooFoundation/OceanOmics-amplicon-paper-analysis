#==============================================================================
# Script to represent "quasi" abundance based on the replicates:
# 0/5, 1/5, ..., 5/5 for presence
#
# Seb Rauschert
#==============================================================================

library(tidyverse)
library(phyloseq)

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

RS_pooled     <- readRDS('data/phyloseq_objects_decontaminated/CoCosV10I_16S_phyloseq_nt_decontaminated.rds')
RS_individual <- readRDS('data/phyloseq_objects_decontaminated/CoCosV10I_16S_phyloseq_nt_FALSE_decontaminated.rds')
RS_pseudo     <- readRDS('data/phyloseq_objects_decontaminated/CoCosV10I_16S_phyloseq_nt_pseudo_decontaminated.rds')

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


pooled     <- format_quasi_abundance(RS_pooled) %>% mutate(Option = 'pooled')
individual <- format_quasi_abundance(RS_individual) %>% mutate(Option = 'individual')
pseudo     <- format_quasi_abundance(RS_pseudo) %>% mutate(Option = 'pseudo')

# rbind(pooled, individual, pseudo) %>%
#   #mutate(Presence = factor(ifelse(Abundance==0, "absent", "present"))) %>%
#   #filter(!is.na(Presence)) 
#   ggplot(aes(x=site_id, y=species, fill = Abundance)) +
#   geom_tile() +
#   theme_bw(base_size = 8) +
#   facet_wrap(~Option) +
#   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
#   scale_fill_gradient(low = "#8FBCBB", high = "#B48EAD")
#   



family <- 
rbind(pooled, individual, pseudo) %>%
  group_by(site_id,family, Option) %>% 
  summarise(Abundance_family = mean(Abundance)) %>%
  #mutate(Presence = factor(ifelse(Abundance==0, "absent", "present"))) %>%
  #filter(!is.na(Presence)) 
  ggplot(aes(x=site_id, y=family, fill = Abundance_family)) +
  geom_tile() +
  theme_bw(base_size = 8) +
  facet_wrap(~Option) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  scale_fill_gradient(low = "#2E3440", high = "#B48EAD")

