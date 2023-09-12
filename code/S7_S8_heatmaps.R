#==============================================================================
# ASV count heatmaps for paper: Supplement 7 & 8

#==============================================================================
# Setup
#==============================================================================

library(tidyverse)
library(tidyr)
library(phyloseq)
library(cowplot)
library(stringr)
WIDTH      <- 16
HEIGHT     <- 9
CCSPLITPAT <- "_(?!.*_)"
CCGREPPAT  <- "_"
RSSPLITPAT <- "_(?!.*_)"
RSGREPPAT  <- "_"
NWSPLITPAT <- "(?<=\\d)(?=[[:alpha:]])"
NWGREPPAT  <- "[A-Za-z]$"
N          <- 20

#==============================================================================
# Functions
#==============================================================================

# Split a string using a pattern
split_string <- function(x, pattern) {
  split_x <- strsplit(x, pattern, perl=TRUE)
  return(split_x[[1]][1])
}

# This function is used to get presence/absence
add_zero_one <- function(x) {
  ifelse(x > 0, 1, 0)
}

# This function converts replicates to sites, it lengthens the df, 
# it can replace the ASV column with values from a different column
# and it can convert counts to presence absence
asv_counts_longer <- function(phyloseq, split_pattern, grep_pattern, col, pres_absence){
  # Extract data tables 
  otu <- rownames_to_column(as.data.frame(otu_table(phyloseq)@.Data, var = "ASV")) %>% 
    dplyr::rename(ASV = rowname) %>%
    select(-(contains(c('WC', 'BC', 'EB', 'DI', 'Ext', 'FC'))))
  
  # CODE to use LCAs instead of ASVs
  tax <- rownames_to_column(as.data.frame(tax_table(phyloseq)@.Data, var = "ASV")) %>% 
    dplyr::rename(ASV = rowname)
  merged_df <- merge(otu, tax[, c("ASV", col)], by = "ASV", all.x = TRUE)
  otu$ASV   <- merged_df[[col]]
  otu       <- otu[!(otu$ASV %in% c("dropped", NA)), ]
  
  # Get the site names from the sample names
  sample_names <- names(otu)
  sample_names <- sample_names[-grep('ASV', sample_names)]
  sites        <- unique(unlist(lapply(sample_names, split_string, pattern=split_pattern)))
  #sites <- sample_names
  
  # Create a new df with the sites instead of samples
  new_otu <- data.frame(matrix(NA, nrow = nrow(otu), ncol = length(sites)))
  colnames(new_otu) <- sites
  if (pres_absence) {
    for(i in sites){
      otu_sample <- otu[,grep(paste0('^', i, grep_pattern), names(otu))]
      otu_sample <- data.frame(lapply(otu_sample, add_zero_one)) 
      new_otu[,i] <- rowSums(otu_sample)
      #new_otu[,i] <- otu[,i]
    }
  } else {
    for(i in sites){
      otu_sample <- otu[,grep(paste0('^', i, grep_pattern), names(otu))]
      #otu_sample <- data.frame(lapply(otu_sample, add_zero_one)) 
      new_otu[,i] <- rowSums(otu_sample)
      #new_otu[,i] <- otu[,i]
    }
  }
  
  # Add the ASV column to the left of the new df
  asv_column <- otu[["ASV"]]
  new_otu <- cbind(ASV = asv_column, new_otu)
  
  # Transpose the df
  trans_otu           <- as.data.frame(t(new_otu))
  colnames(trans_otu) <- trans_otu[1, ]
  trans_otu           <- trans_otu[-1, ]
  asvs                <- names(trans_otu)
  
  # Make sure the columns are integer columns
  for (i in 1:(ncol(trans_otu))) {
    trans_otu[[i]]    <- as.integer(trans_otu[[i]])
  }
  
  trans_otu$sample    <- rownames(trans_otu)
  
  # Reformat to long format
  trans_otu %>% 
    pivot_longer(c(asvs)) %>%
    dplyr::rename(asv = 'name',
                  count = 'value') -> otu_long
  return(otu_long)
}

# This function takes three dfs and creates a heatmap with an optional fill column (count or log_count)
# This function can also, optionally, hide the y label
create_heatmap <- function(pooled, independent, pseudo, col, ylabels){
  if (ylabels) {
    rbind(pooled, independent, pseudo) %>%
      ggplot(aes(x=sample, y=asv, fill = .data[[col]])) +
      geom_tile() +
      theme_classic(base_size = 8) +
      facet_wrap(~Option) +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
      scale_fill_gradient(low = "#8FBCBB", high = "#B48EAD")
  } else {
    rbind(pooled, independent, pseudo) %>%
      ggplot(aes(x=sample, y=asv, fill = .data[[col]])) +
      geom_tile() +
      theme_classic(base_size = 8) +
      facet_wrap(~Option) +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), axis.text.y = element_blank()) +
      scale_fill_gradient(low = "#8FBCBB", high = "#B48EAD")
  }
}

# This function normalizes the asv counts
normalize_asvs <- function(df) {
  df_log <- df %>%
    mutate(log_count = log(count + 1))
  
  return(df_log)
}

# This function can create a df with the top/bottom n asvs
# This function can optionally normalize the counts
subset_asvs <- function(df, section, normalize = FALSE, n = 20) {
  grouped_df <- df %>%
    group_by(asv) %>%
    summarize(total_count = sum(count))
  
  if (section == "top") {
    subset_asvs <- grouped_df %>%
      arrange(desc(total_count)) %>%
      head(n) %>%
      select(asv)
  } else if (section == "bottom") {
    subset_asvs <- grouped_df %>%
      arrange(desc(total_count)) %>%
      tail(n) %>%
      select(asv)
  } else {
    subset_asvs <- grouped_df %>%
      arrange(desc(total_count)) %>%
      select(asv)
  }
  
  subset_df <- df %>%
    filter(asv %in% subset_asvs$asv)
  
  if (normalize) {
    subset_df_log <- normalize_asvs(subset_df)
    
    return(subset_df_log)
    
  } else {
    return(subset_df)
  }
}

#==============================================================================
# Import data
#==============================================================================

Cocos_indepe <- readRDS('data/phyloseq_objects/CoCosV10I_16S_phyloseq_nt_FALSE_decontaminated.rds')
Cocos_pooled <- readRDS('data/phyloseq_objects/CoCosV10I_16S_phyloseq_nt_TRUE_decontaminated.rds')
Cocos_pseudo <- readRDS('data/phyloseq_objects/CoCosV10I_16S_phyloseq_nt_pseudo_decontaminated.rds')

NW_pooled    <- readRDS('data/phyloseq_objects/NWWA_16S_phyloseq_nt_TRUE_decontam.rds')
NW_indepe    <- readRDS('data/phyloseq_objects/NWWA_16S_phyloseq_nt_FALSE_decontam.rds')
NW_pseudo    <- readRDS('data/phyloseq_objects/NWWA_16S_phyloseq_nt_pseudo_decontam.rds')
NW_pseudo    <- subset_samples(NW_pseudo, sample_names(NW_pseudo) != "71a")

RS_pooled    <- readRDS('data/phyloseq_objects/RS21AUG_16S_phyloseq_nt_TRUE_decontam.rds')
RS_indepe    <- readRDS('data/phyloseq_objects/RS21AUG_16S_phyloseq_nt_FALSE_decontam.rds')
RS_pseudo    <- readRDS('data/phyloseq_objects/RS21AUG_16S_phyloseq_nt_pseudo_decontam.rds')
RS_pooled    <- subset_samples(RS_pooled, sample_names(RS_pooled) != "RS1_ME_S4_1_2_")
RS_indepe    <- subset_samples(RS_indepe, sample_names(RS_indepe) != "RS1_ME_S4_1_2_")
RS_pseudo    <- subset_samples(RS_pseudo, sample_names(RS_pseudo) != "RS1_ME_S4_1_2_")

#==============================================================================
# CococsI
#==============================================================================

# Make all long dfs
pooled_presnc_species <- asv_counts_longer(Cocos_pooled, CCSPLITPAT, CCGREPPAT, "species", TRUE) %>% mutate(Option = 'Pooled')
indepe_presnc_species <- asv_counts_longer(Cocos_indepe, CCSPLITPAT, CCGREPPAT, "species", TRUE) %>% mutate(Option = 'Independent')
pseudo_presnc_species <- asv_counts_longer(Cocos_pseudo, CCSPLITPAT, CCGREPPAT, "species", TRUE) %>% mutate(Option = 'Pseudo')
pooled_presnc_asvs    <- asv_counts_longer(Cocos_pooled, CCSPLITPAT, CCGREPPAT, "ASV_sequence", TRUE) %>% mutate(Option = 'Pooled')
indepe_presnc_asvs    <- asv_counts_longer(Cocos_indepe, CCSPLITPAT, CCGREPPAT, "ASV_sequence", TRUE) %>% mutate(Option = 'Independent')
pseudo_presnc_asvs    <- asv_counts_longer(Cocos_pseudo, CCSPLITPAT, CCGREPPAT, "ASV_sequence", TRUE) %>% mutate(Option = 'Pseudo')
pooled_counts_species <- asv_counts_longer(Cocos_pooled, CCSPLITPAT, CCGREPPAT, "species", FALSE) %>% mutate(Option = 'Pooled')
indepe_counts_species <- asv_counts_longer(Cocos_indepe, CCSPLITPAT, CCGREPPAT, "species", FALSE) %>% mutate(Option = 'Independent')
pseudo_counts_species <- asv_counts_longer(Cocos_pseudo, CCSPLITPAT, CCGREPPAT, "species", FALSE) %>% mutate(Option = 'Pseudo')
pooled_counts_asvs    <- asv_counts_longer(Cocos_pooled, CCSPLITPAT, CCGREPPAT, "ASV_sequence", FALSE) %>% mutate(Option = 'Pooled')
indepe_counts_asvs    <- asv_counts_longer(Cocos_indepe, CCSPLITPAT, CCGREPPAT, "ASV_sequence", FALSE) %>% mutate(Option = 'Independent')
pseudo_counts_asvs    <- asv_counts_longer(Cocos_pseudo, CCSPLITPAT, CCGREPPAT, "ASV_sequence", FALSE) %>% mutate(Option = 'Pseudo')

# We need to start with the top/bottom independent dfs; some of these will be normalized
top_indepe_presnc_species <- subset_asvs(indepe_presnc_species, "top", FALSE, N)
bot_indepe_presnc_species <- subset_asvs(indepe_presnc_species, "bottom", FALSE, N)
all_indepe_presnc_species <- subset_asvs(indepe_presnc_species, "all", FALSE, N)
top_indepe_presnc_asvs    <- subset_asvs(indepe_presnc_asvs, "top", FALSE, N)
bot_indepe_presnc_asvs    <- subset_asvs(indepe_presnc_asvs, "bottom", FALSE, N)
all_indepe_presnc_asvs    <- subset_asvs(indepe_presnc_asvs, "all", FALSE, N)
top_indepe_counts_species <- subset_asvs(indepe_counts_species, "top", TRUE, N)
bot_indepe_counts_species <- subset_asvs(indepe_counts_species, "bottom", TRUE, N)
all_indepe_counts_species <- subset_asvs(indepe_counts_species, "all", TRUE, N)
top_indepe_counts_asvs    <- subset_asvs(indepe_counts_asvs, "top", TRUE, N)
bot_indepe_counts_asvs    <- subset_asvs(indepe_counts_asvs, "bottom", TRUE, N)
all_indepe_counts_asvs    <- subset_asvs(indepe_counts_asvs, "all", TRUE, N)

# We need to filter pooled and pseudo to match independent
top_pooled_presnc_species <- filter(pooled_presnc_species, asv %in% top_indepe_presnc_species$asv)
top_pseudo_presnc_species <- filter(pseudo_presnc_species, asv %in% top_indepe_presnc_species$asv)
bot_pooled_presnc_species <- filter(pooled_presnc_species, asv %in% bot_indepe_presnc_species$asv)
bot_pseudo_presnc_species <- filter(pseudo_presnc_species, asv %in% bot_indepe_presnc_species$asv)
all_pooled_presnc_species <- filter(pooled_presnc_species, asv %in% all_indepe_presnc_species$asv)
all_pseudo_presnc_species <- filter(pseudo_presnc_species, asv %in% all_indepe_presnc_species$asv)
top_pooled_presnc_asvs    <- filter(pooled_presnc_asvs, asv %in% top_indepe_presnc_asvs$asv)
top_pseudo_presnc_asvs    <- filter(pseudo_presnc_asvs, asv %in% top_indepe_presnc_asvs$asv)
bot_pooled_presnc_asvs    <- filter(pooled_presnc_asvs, asv %in% bot_indepe_presnc_asvs$asv)
bot_pseudo_presnc_asvs    <- filter(pseudo_presnc_asvs, asv %in% bot_indepe_presnc_asvs$asv)
all_pooled_presnc_asvs    <- filter(pooled_presnc_asvs, asv %in% all_indepe_presnc_asvs$asv)
all_pseudo_presnc_asvs    <- filter(pseudo_presnc_asvs, asv %in% all_indepe_presnc_asvs$asv)
top_pooled_counts_species <- filter(pooled_counts_species, asv %in% top_indepe_counts_species$asv)
top_pseudo_counts_species <- filter(pseudo_counts_species, asv %in% top_indepe_counts_species$asv)
bot_pooled_counts_species <- filter(pooled_counts_species, asv %in% bot_indepe_counts_species$asv)
bot_pseudo_counts_species <- filter(pseudo_counts_species, asv %in% bot_indepe_counts_species$asv)
all_pooled_counts_species <- filter(pooled_counts_species, asv %in% all_indepe_counts_species$asv)
all_pseudo_counts_species <- filter(pseudo_counts_species, asv %in% all_indepe_counts_species$asv)
top_pooled_counts_asvs    <- filter(pooled_counts_asvs, asv %in% top_indepe_counts_asvs$asv)
top_pseudo_counts_asvs    <- filter(pseudo_counts_asvs, asv %in% top_indepe_counts_asvs$asv)
bot_pooled_counts_asvs    <- filter(pooled_counts_asvs, asv %in% bot_indepe_counts_asvs$asv)
bot_pseudo_counts_asvs    <- filter(pseudo_counts_asvs, asv %in% bot_indepe_counts_asvs$asv)
all_pooled_counts_asvs    <- filter(pooled_counts_asvs, asv %in% all_indepe_counts_asvs$asv)
all_pseudo_counts_asvs    <- filter(pseudo_counts_asvs, asv %in% all_indepe_counts_asvs$asv)

# We need to normalize the pooled and pseudo counts dfs; independent already got normalized
top_pooled_counts_species <- normalize_asvs(top_pooled_counts_species)
top_pseudo_counts_species <- normalize_asvs(top_pseudo_counts_species)
bot_pooled_counts_species <- normalize_asvs(bot_pooled_counts_species)
bot_pseudo_counts_species <- normalize_asvs(bot_pseudo_counts_species)
all_pooled_counts_species <- normalize_asvs(all_pooled_counts_species)
all_pseudo_counts_species <- normalize_asvs(all_pseudo_counts_species)
top_pooled_counts_asvs    <- normalize_asvs(top_pooled_counts_asvs)
top_pseudo_counts_asvs    <- normalize_asvs(top_pseudo_counts_asvs)
bot_pooled_counts_asvs    <- normalize_asvs(bot_pooled_counts_asvs)
bot_pseudo_counts_asvs    <- normalize_asvs(bot_pseudo_counts_asvs)
all_pooled_counts_asvs    <- normalize_asvs(all_pooled_counts_asvs)
all_pseudo_counts_asvs    <- normalize_asvs(all_pseudo_counts_asvs)

# Now create the heatmaps
pdf(file="pres_absence_species_top_independant_cocos.pdf", width=WIDTH, height=HEIGHT)
create_heatmap(top_pooled_presnc_species, top_indepe_presnc_species, top_pseudo_presnc_species, "count", TRUE)
dev.off()
pdf(file="pres_absence_species_bottom_independant_cocos.pdf", width=WIDTH, height=HEIGHT)
create_heatmap(bot_pooled_presnc_species, bot_indepe_presnc_species, bot_pseudo_presnc_species, "count", TRUE)
dev.off()
pdf(file="pres_absence_species_all_independant_cocos.pdf", width=WIDTH, height=HEIGHT)
create_heatmap(all_pooled_presnc_species, all_indepe_presnc_species, all_pseudo_presnc_species, "count", TRUE)
dev.off()
pdf(file="pres_absence_asv_top_independant_cocos.pdf", width=WIDTH, height=HEIGHT)
cocos_top <- create_heatmap(top_pooled_presnc_asvs, top_indepe_presnc_asvs, top_pseudo_presnc_asvs, "count", FALSE)
dev.off()
pdf(file="pres_absence_asv_bottom_independant_cocos.pdf", width=WIDTH, height=HEIGHT)
cocos_bot <- create_heatmap(bot_pooled_presnc_asvs, bot_indepe_presnc_asvs, bot_pseudo_presnc_asvs, "count", FALSE)
dev.off()
pdf(file="pres_absence_asv_all_independant_cocos.pdf", width=WIDTH, height=HEIGHT)
create_heatmap(all_pooled_presnc_asvs, all_indepe_presnc_asvs, all_pseudo_presnc_asvs, "count", FALSE)
dev.off()
pdf(file="log_counts_species_top_independant_cocos.pdf", width=WIDTH, height=HEIGHT)
create_heatmap(top_pooled_counts_species, top_indepe_counts_species, top_pseudo_counts_species, "log_count", TRUE)
dev.off()
pdf(file="log_counts_species_bottom_independant_cocos.pdf", width=WIDTH, height=HEIGHT)
create_heatmap(bot_pooled_counts_species, bot_indepe_counts_species, bot_pseudo_counts_species, "log_count", TRUE)
dev.off()
pdf(file="log_counts_species_all_independant_cocos.pdf", width=WIDTH, height=HEIGHT)
create_heatmap(all_pooled_counts_species, all_indepe_counts_species, all_pseudo_counts_species, "log_count", TRUE)
dev.off()
pdf(file="log_counts_asv_top_independant_cocos.pdf", width=WIDTH, height=HEIGHT)
create_heatmap(top_pooled_counts_asvs, top_indepe_counts_asvs, top_pseudo_counts_asvs, "log_count", FALSE)
dev.off()
pdf(file="log_counts_asv_bottom_independant_cocos.pdf", width=WIDTH, height=HEIGHT)
create_heatmap(bot_pooled_counts_asvs, bot_indepe_counts_asvs, bot_pseudo_counts_asvs, "log_count", FALSE)
dev.off()
pdf(file="log_counts_asv_all_independant_cocos.pdf", width=WIDTH, height=HEIGHT)
create_heatmap(all_pooled_counts_asvs, all_indepe_counts_asvs, all_pseudo_counts_asvs, "log_count", FALSE)
dev.off()

#==============================================================================
# North west WA data
#==============================================================================

# Make all long dfs
pooled_presnc_species <- asv_counts_longer(NW_pooled, NWSPLITPAT, NWGREPPAT, "species", TRUE) %>% mutate(Option = 'Pooled')
indepe_presnc_species <- asv_counts_longer(NW_indepe, NWSPLITPAT, NWGREPPAT, "species", TRUE) %>% mutate(Option = 'Independent')
pseudo_presnc_species <- asv_counts_longer(NW_pseudo, NWSPLITPAT, NWGREPPAT, "species", TRUE) %>% mutate(Option = 'Pseudo')
pooled_presnc_asvs    <- asv_counts_longer(NW_pooled, NWSPLITPAT, NWGREPPAT, "ASV_sequence", TRUE) %>% mutate(Option = 'Pooled')
indepe_presnc_asvs    <- asv_counts_longer(NW_indepe, NWSPLITPAT, NWGREPPAT, "ASV_sequence", TRUE) %>% mutate(Option = 'Independent')
pseudo_presnc_asvs    <- asv_counts_longer(NW_pseudo, NWSPLITPAT, NWGREPPAT, "ASV_sequence", TRUE) %>% mutate(Option = 'Pseudo')
pooled_counts_species <- asv_counts_longer(NW_pooled, NWSPLITPAT, NWGREPPAT, "species", FALSE) %>% mutate(Option = 'Pooled')
indepe_counts_species <- asv_counts_longer(NW_indepe, NWSPLITPAT, NWGREPPAT, "species", FALSE) %>% mutate(Option = 'Independent')
pseudo_counts_species <- asv_counts_longer(NW_pseudo, NWSPLITPAT, NWGREPPAT, "species", FALSE) %>% mutate(Option = 'Pseudo')
pooled_counts_asvs    <- asv_counts_longer(NW_pooled, NWSPLITPAT, NWGREPPAT, "ASV_sequence", FALSE) %>% mutate(Option = 'Pooled')
indepe_counts_asvs    <- asv_counts_longer(NW_indepe, NWSPLITPAT, NWGREPPAT, "ASV_sequence", FALSE) %>% mutate(Option = 'Independent')
pseudo_counts_asvs    <- asv_counts_longer(NW_pseudo, NWSPLITPAT, NWGREPPAT, "ASV_sequence", FALSE) %>% mutate(Option = 'Pseudo')

# We need to start with the top/bottom independent dfs; some of these will be normalized
top_indepe_presnc_species <- subset_asvs(indepe_presnc_species, "top", FALSE, N)
bot_indepe_presnc_species <- subset_asvs(indepe_presnc_species, "bottom", FALSE, N)
all_indepe_presnc_species <- subset_asvs(indepe_presnc_species, "all", FALSE, N)
top_indepe_presnc_asvs    <- subset_asvs(indepe_presnc_asvs, "top", FALSE, N)
bot_indepe_presnc_asvs    <- subset_asvs(indepe_presnc_asvs, "bottom", FALSE, N)
all_indepe_presnc_asvs    <- subset_asvs(indepe_presnc_asvs, "all", FALSE, N)
top_indepe_counts_species <- subset_asvs(indepe_counts_species, "top", TRUE, N)
bot_indepe_counts_species <- subset_asvs(indepe_counts_species, "bottom", TRUE, N)
all_indepe_counts_species <- subset_asvs(indepe_counts_species, "all", TRUE, N)
top_indepe_counts_asvs    <- subset_asvs(indepe_counts_asvs, "top", TRUE, N)
bot_indepe_counts_asvs    <- subset_asvs(indepe_counts_asvs, "bottom", TRUE, N)
all_indepe_counts_asvs    <- subset_asvs(indepe_counts_asvs, "all", TRUE, N)

# We need to filter pooled and pseudo to match independent
top_pooled_presnc_species <- filter(pooled_presnc_species, asv %in% top_indepe_presnc_species$asv)
top_pseudo_presnc_species <- filter(pseudo_presnc_species, asv %in% top_indepe_presnc_species$asv)
bot_pooled_presnc_species <- filter(pooled_presnc_species, asv %in% bot_indepe_presnc_species$asv)
bot_pseudo_presnc_species <- filter(pseudo_presnc_species, asv %in% bot_indepe_presnc_species$asv)
all_pooled_presnc_species <- filter(pooled_presnc_species, asv %in% all_indepe_presnc_species$asv)
all_pseudo_presnc_species <- filter(pseudo_presnc_species, asv %in% all_indepe_presnc_species$asv)
top_pooled_presnc_asvs    <- filter(pooled_presnc_asvs, asv %in% top_indepe_presnc_asvs$asv)
top_pseudo_presnc_asvs    <- filter(pseudo_presnc_asvs, asv %in% top_indepe_presnc_asvs$asv)
bot_pooled_presnc_asvs    <- filter(pooled_presnc_asvs, asv %in% bot_indepe_presnc_asvs$asv)
bot_pseudo_presnc_asvs    <- filter(pseudo_presnc_asvs, asv %in% bot_indepe_presnc_asvs$asv)
all_pooled_presnc_asvs    <- filter(pooled_presnc_asvs, asv %in% all_indepe_presnc_asvs$asv)
all_pseudo_presnc_asvs    <- filter(pseudo_presnc_asvs, asv %in% all_indepe_presnc_asvs$asv)
top_pooled_counts_species <- filter(pooled_counts_species, asv %in% top_indepe_counts_species$asv)
top_pseudo_counts_species <- filter(pseudo_counts_species, asv %in% top_indepe_counts_species$asv)
bot_pooled_counts_species <- filter(pooled_counts_species, asv %in% bot_indepe_counts_species$asv)
bot_pseudo_counts_species <- filter(pseudo_counts_species, asv %in% bot_indepe_counts_species$asv)
all_pooled_counts_species <- filter(pooled_counts_species, asv %in% all_indepe_counts_species$asv)
all_pseudo_counts_species <- filter(pseudo_counts_species, asv %in% all_indepe_counts_species$asv)
top_pooled_counts_asvs    <- filter(pooled_counts_asvs, asv %in% top_indepe_counts_asvs$asv)
top_pseudo_counts_asvs    <- filter(pseudo_counts_asvs, asv %in% top_indepe_counts_asvs$asv)
bot_pooled_counts_asvs    <- filter(pooled_counts_asvs, asv %in% bot_indepe_counts_asvs$asv)
bot_pseudo_counts_asvs    <- filter(pseudo_counts_asvs, asv %in% bot_indepe_counts_asvs$asv)
all_pooled_counts_asvs    <- filter(pooled_counts_asvs, asv %in% all_indepe_counts_asvs$asv)
all_pseudo_counts_asvs    <- filter(pseudo_counts_asvs, asv %in% all_indepe_counts_asvs$asv)

# We need to normalize the pooled and pseudo counts dfs; independent already got normalized
top_pooled_counts_species <- normalize_asvs(top_pooled_counts_species)
top_pseudo_counts_species <- normalize_asvs(top_pseudo_counts_species)
bot_pooled_counts_species <- normalize_asvs(bot_pooled_counts_species)
bot_pseudo_counts_species <- normalize_asvs(bot_pseudo_counts_species)
all_pooled_counts_species <- normalize_asvs(all_pooled_counts_species)
all_pseudo_counts_species <- normalize_asvs(all_pseudo_counts_species)
top_pooled_counts_asvs    <- normalize_asvs(top_pooled_counts_asvs)
top_pseudo_counts_asvs    <- normalize_asvs(top_pseudo_counts_asvs)
bot_pooled_counts_asvs    <- normalize_asvs(bot_pooled_counts_asvs)
bot_pseudo_counts_asvs    <- normalize_asvs(bot_pseudo_counts_asvs)
all_pooled_counts_asvs    <- normalize_asvs(all_pooled_counts_asvs)
all_pseudo_counts_asvs    <- normalize_asvs(all_pseudo_counts_asvs)

# Now create the heatmaps
pdf(file="pres_absence_species_top_independant_nw.pdf", width=WIDTH, height=HEIGHT)
create_heatmap(top_pooled_presnc_species, top_indepe_presnc_species, top_pseudo_presnc_species, "count", TRUE)
dev.off()
pdf(file="pres_absence_species_bottom_independant_nw.pdf", width=WIDTH, height=HEIGHT)
create_heatmap(bot_pooled_presnc_species, bot_indepe_presnc_species, bot_pseudo_presnc_species, "count", TRUE)
dev.off()
pdf(file="pres_absence_species_all_independant_nw.pdf", width=WIDTH, height=HEIGHT)
create_heatmap(all_pooled_presnc_species, all_indepe_presnc_species, all_pseudo_presnc_species, "count", TRUE)
dev.off()
pdf(file="pres_absence_asv_top_independant_nw.pdf", width=WIDTH, height=HEIGHT)
nw_top <- create_heatmap(top_pooled_presnc_asvs, top_indepe_presnc_asvs, top_pseudo_presnc_asvs, "count", FALSE)
dev.off()
pdf(file="pres_absence_asv_bottom_independant_nw.pdf", width=WIDTH, height=HEIGHT)
nw_bot <- create_heatmap(bot_pooled_presnc_asvs, bot_indepe_presnc_asvs, bot_pseudo_presnc_asvs, "count", FALSE)
dev.off()
pdf(file="pres_absence_asv_all_independant_nw.pdf", width=WIDTH, height=HEIGHT)
create_heatmap(all_pooled_presnc_asvs, all_indepe_presnc_asvs, all_pseudo_presnc_asvs, "count", FALSE)
dev.off()
pdf(file="log_counts_species_top_independant_nw.pdf", width=WIDTH, height=HEIGHT)
create_heatmap(top_pooled_counts_species, top_indepe_counts_species, top_pseudo_counts_species, "log_count", TRUE)
dev.off()
pdf(file="log_counts_species_bottom_independant_nw.pdf", width=WIDTH, height=HEIGHT)
create_heatmap(bot_pooled_counts_species, bot_indepe_counts_species, bot_pseudo_counts_species, "log_count", TRUE)
dev.off()
pdf(file="log_counts_species_all_independant_nw.pdf", width=WIDTH, height=HEIGHT)
create_heatmap(all_pooled_counts_species, all_indepe_counts_species, all_pseudo_counts_species, "log_count", TRUE)
dev.off()
pdf(file="log_counts_asv_top_independant_nw.pdf", width=WIDTH, height=HEIGHT)
create_heatmap(top_pooled_counts_asvs, top_indepe_counts_asvs, top_pseudo_counts_asvs, "log_count", FALSE)
dev.off()
pdf(file="log_counts_asv_bottom_independant_nw.pdf", width=WIDTH, height=HEIGHT)
create_heatmap(bot_pooled_counts_asvs, bot_indepe_counts_asvs, bot_pseudo_counts_asvs, "log_count", FALSE)
dev.off()
pdf(file="log_counts_asv_all_independant_nw.pdf", width=WIDTH, height=HEIGHT)
create_heatmap(all_pooled_counts_asvs, all_indepe_counts_asvs, all_pseudo_counts_asvs, "log_count", FALSE)
dev.off()

#==============================================================================
# Rowley Shoals data
#==============================================================================

# Make all long dfs
pooled_presnc_species <- asv_counts_longer(RS_pooled, RSSPLITPAT, RSGREPPAT, "species", TRUE) %>% mutate(Option = 'Pooled')
indepe_presnc_species <- asv_counts_longer(RS_indepe, RSSPLITPAT, RSGREPPAT, "species", TRUE) %>% mutate(Option = 'Independent')
pseudo_presnc_species <- asv_counts_longer(RS_pseudo, RSSPLITPAT, RSGREPPAT, "species", TRUE) %>% mutate(Option = 'Pseudo')
pooled_presnc_asvs    <- asv_counts_longer(RS_pooled, RSSPLITPAT, RSGREPPAT, "ASV_sequence", TRUE) %>% mutate(Option = 'Pooled')
indepe_presnc_asvs    <- asv_counts_longer(RS_indepe, RSSPLITPAT, RSGREPPAT, "ASV_sequence", TRUE) %>% mutate(Option = 'Independent')
pseudo_presnc_asvs    <- asv_counts_longer(RS_pseudo, RSSPLITPAT, RSGREPPAT, "ASV_sequence", TRUE) %>% mutate(Option = 'Pseudo')
pooled_counts_species <- asv_counts_longer(RS_pooled, RSSPLITPAT, RSGREPPAT, "species", FALSE) %>% mutate(Option = 'Pooled')
indepe_counts_species <- asv_counts_longer(RS_indepe, RSSPLITPAT, RSGREPPAT, "species", FALSE) %>% mutate(Option = 'Independent')
pseudo_counts_species <- asv_counts_longer(RS_pseudo, RSSPLITPAT, RSGREPPAT, "species", FALSE) %>% mutate(Option = 'Pseudo')
pooled_counts_asvs    <- asv_counts_longer(RS_pooled, RSSPLITPAT, RSGREPPAT, "ASV_sequence", FALSE) %>% mutate(Option = 'Pooled')
indepe_counts_asvs    <- asv_counts_longer(RS_indepe, RSSPLITPAT, RSGREPPAT, "ASV_sequence", FALSE) %>% mutate(Option = 'Independent')
pseudo_counts_asvs    <- asv_counts_longer(RS_pseudo, RSSPLITPAT, RSGREPPAT, "ASV_sequence", FALSE) %>% mutate(Option = 'Pseudo')

# We need to start with the top/bottom independent dfs; some of these will be normalized
top_indepe_presnc_species <- subset_asvs(indepe_presnc_species, "top", FALSE, N)
bot_indepe_presnc_species <- subset_asvs(indepe_presnc_species, "bottom", FALSE, N)
all_indepe_presnc_species <- subset_asvs(indepe_presnc_species, "all", FALSE, N)
top_indepe_presnc_asvs    <- subset_asvs(indepe_presnc_asvs, "top", FALSE, N)
bot_indepe_presnc_asvs    <- subset_asvs(indepe_presnc_asvs, "bottom", FALSE, N)
all_indepe_presnc_asvs    <- subset_asvs(indepe_presnc_asvs, "all", FALSE, N)
top_indepe_counts_species <- subset_asvs(indepe_counts_species, "top", TRUE, N)
bot_indepe_counts_species <- subset_asvs(indepe_counts_species, "bottom", TRUE, N)
all_indepe_counts_species <- subset_asvs(indepe_counts_species, "all", TRUE, N)
top_indepe_counts_asvs    <- subset_asvs(indepe_counts_asvs, "top", TRUE, N)
bot_indepe_counts_asvs    <- subset_asvs(indepe_counts_asvs, "bottom", TRUE, N)
all_indepe_counts_asvs    <- subset_asvs(indepe_counts_asvs, "all", TRUE, N)

# We need to filter pooled and pseudo to match independent
top_pooled_presnc_species <- filter(pooled_presnc_species, asv %in% top_indepe_presnc_species$asv)
top_pseudo_presnc_species <- filter(pseudo_presnc_species, asv %in% top_indepe_presnc_species$asv)
bot_pooled_presnc_species <- filter(pooled_presnc_species, asv %in% bot_indepe_presnc_species$asv)
bot_pseudo_presnc_species <- filter(pseudo_presnc_species, asv %in% bot_indepe_presnc_species$asv)
all_pooled_presnc_species <- filter(pooled_presnc_species, asv %in% all_indepe_presnc_species$asv)
all_pseudo_presnc_species <- filter(pseudo_presnc_species, asv %in% all_indepe_presnc_species$asv)
top_pooled_presnc_asvs    <- filter(pooled_presnc_asvs, asv %in% top_indepe_presnc_asvs$asv)
top_pseudo_presnc_asvs    <- filter(pseudo_presnc_asvs, asv %in% top_indepe_presnc_asvs$asv)
bot_pooled_presnc_asvs    <- filter(pooled_presnc_asvs, asv %in% bot_indepe_presnc_asvs$asv)
bot_pseudo_presnc_asvs    <- filter(pseudo_presnc_asvs, asv %in% bot_indepe_presnc_asvs$asv)
all_pooled_presnc_asvs    <- filter(pooled_presnc_asvs, asv %in% all_indepe_presnc_asvs$asv)
all_pseudo_presnc_asvs    <- filter(pseudo_presnc_asvs, asv %in% all_indepe_presnc_asvs$asv)
top_pooled_counts_species <- filter(pooled_counts_species, asv %in% top_indepe_counts_species$asv)
top_pseudo_counts_species <- filter(pseudo_counts_species, asv %in% top_indepe_counts_species$asv)
bot_pooled_counts_species <- filter(pooled_counts_species, asv %in% bot_indepe_counts_species$asv)
bot_pseudo_counts_species <- filter(pseudo_counts_species, asv %in% bot_indepe_counts_species$asv)
all_pooled_counts_species <- filter(pooled_counts_species, asv %in% all_indepe_counts_species$asv)
all_pseudo_counts_species <- filter(pseudo_counts_species, asv %in% all_indepe_counts_species$asv)
top_pooled_counts_asvs    <- filter(pooled_counts_asvs, asv %in% top_indepe_counts_asvs$asv)
top_pseudo_counts_asvs    <- filter(pseudo_counts_asvs, asv %in% top_indepe_counts_asvs$asv)
bot_pooled_counts_asvs    <- filter(pooled_counts_asvs, asv %in% bot_indepe_counts_asvs$asv)
bot_pseudo_counts_asvs    <- filter(pseudo_counts_asvs, asv %in% bot_indepe_counts_asvs$asv)
all_pooled_counts_asvs    <- filter(pooled_counts_asvs, asv %in% all_indepe_counts_asvs$asv)
all_pseudo_counts_asvs    <- filter(pseudo_counts_asvs, asv %in% all_indepe_counts_asvs$asv)

# We need to normalize the pooled and pseudo counts dfs; independent already got normalized
top_pooled_counts_species <- normalize_asvs(top_pooled_counts_species)
top_pseudo_counts_species <- normalize_asvs(top_pseudo_counts_species)
bot_pooled_counts_species <- normalize_asvs(bot_pooled_counts_species)
bot_pseudo_counts_species <- normalize_asvs(bot_pseudo_counts_species)
all_pooled_counts_species <- normalize_asvs(all_pooled_counts_species)
all_pseudo_counts_species <- normalize_asvs(all_pseudo_counts_species)
top_pooled_counts_asvs    <- normalize_asvs(top_pooled_counts_asvs)
top_pseudo_counts_asvs    <- normalize_asvs(top_pseudo_counts_asvs)
bot_pooled_counts_asvs    <- normalize_asvs(bot_pooled_counts_asvs)
bot_pseudo_counts_asvs    <- normalize_asvs(bot_pseudo_counts_asvs)
all_pooled_counts_asvs    <- normalize_asvs(all_pooled_counts_asvs)
all_pseudo_counts_asvs    <- normalize_asvs(all_pseudo_counts_asvs)

# Now create the heatmaps
pdf(file="pres_absence_species_top_independant_rs.pdf", width=WIDTH, height=HEIGHT)
create_heatmap(top_pooled_presnc_species, top_indepe_presnc_species, top_pseudo_presnc_species, "count", TRUE)
dev.off()
pdf(file="pres_absence_species_bottom_independant_rs.pdf", width=WIDTH, height=HEIGHT)
create_heatmap(bot_pooled_presnc_species, bot_indepe_presnc_species, bot_pseudo_presnc_species, "count", TRUE)
dev.off()
pdf(file="pres_absence_species_all_independant_rs.pdf", width=WIDTH, height=HEIGHT)
create_heatmap(all_pooled_presnc_species, all_indepe_presnc_species, all_pseudo_presnc_species, "count", TRUE)
dev.off()
pdf(file="pres_absence_asv_top_independant_rs.pdf", width=WIDTH, height=HEIGHT)
rs_top <- create_heatmap(top_pooled_presnc_asvs, top_indepe_presnc_asvs, top_pseudo_presnc_asvs, "count", FALSE)
dev.off()
pdf(file="pres_absence_asv_bottom_independant_rs.pdf", width=WIDTH, height=HEIGHT)
rs_bot <- create_heatmap(bot_pooled_presnc_asvs, bot_indepe_presnc_asvs, bot_pseudo_presnc_asvs, "count", FALSE)
dev.off()
pdf(file="pres_absence_asv_all_independant_rs.pdf", width=WIDTH, height=HEIGHT)
create_heatmap(all_pooled_presnc_asvs, all_indepe_presnc_asvs, all_pseudo_presnc_asvs, "count", FALSE)
dev.off()
pdf(file="log_counts_species_top_independant_rs.pdf", width=WIDTH, height=HEIGHT)
create_heatmap(top_pooled_counts_species, top_indepe_counts_species, top_pseudo_counts_species, "log_count", TRUE)
dev.off()
pdf(file="log_counts_species_bottom_independant_rs.pdf", width=WIDTH, height=HEIGHT)
create_heatmap(bot_pooled_counts_species, bot_indepe_counts_species, bot_pseudo_counts_species, "log_count", TRUE)
dev.off()
pdf(file="log_counts_species_all_independant_rs.pdf", width=WIDTH, height=HEIGHT)
create_heatmap(all_pooled_counts_species, all_indepe_counts_species, all_pseudo_counts_species, "log_count", TRUE)
dev.off()
pdf(file="log_counts_asv_top_independant_rs.pdf", width=WIDTH, height=HEIGHT)
create_heatmap(top_pooled_counts_asvs, top_indepe_counts_asvs, top_pseudo_counts_asvs, "log_count", FALSE)
dev.off()
pdf(file="log_counts_asv_bottom_independant_rs.pdf", width=WIDTH, height=HEIGHT)
create_heatmap(bot_pooled_counts_asvs, bot_indepe_counts_asvs, bot_pseudo_counts_asvs, "log_count", FALSE)
dev.off()
pdf(file="log_counts_asv_all_independant_rs.pdf", width=WIDTH, height=HEIGHT)
create_heatmap(all_pooled_counts_asvs, all_indepe_counts_asvs, all_pseudo_counts_asvs, "log_count", FALSE)
dev.off()

plot_grid(cocos_top + ggtitle("Cocos Islands") + ylab("Top 20 ASVs") + theme(axis.text.x = element_blank()) + theme(axis.text.y = element_blank()),
          cocos_bot + ylab("Bottom 20 ASVs") + theme(axis.text.x = element_blank()) + theme(axis.text.y = element_blank()),
          ncol = 1)

plot_grid(rs_top + ggtitle("Rowley Shoals") + ylab("Top 20 ASVs") + theme(axis.text.x = element_blank()) + theme(axis.text.y = element_blank()),
          rs_bot + ylab("Bottom 20 ASVs") + theme(axis.text.x = element_blank()) + theme(axis.text.y = element_blank()),
          ncol = 1)

plot_grid(nw_top + ggtitle("NW WA") + ylab("Top 20 ASVs") + theme(axis.text.x = element_blank()) + theme(axis.text.y = element_blank()),
          nw_bot + ylab("Bottom 20 ASVs") + theme(axis.text.x = element_blank()) + theme(axis.text.y = element_blank()),
          ncol = 1)
