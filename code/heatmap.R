# ASV count heatmaps for paper

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
CCCUTOFF   <- 500
RSSPLITPAT <- "_(?!.*_)"
RSGREPPAT  <- "_"
RSCUTOFF   <- 20000
NWSPLITPAT <- "(?<=\\d)(?=[[:alpha:]])"
NWGREPPAT  <- "[A-Za-z]$"
NWCUTOFF   <- 100000
TOPN       <- 20
BOTTOMN    <- 20


#==============================================================================
# Functions
#==============================================================================

split_string <- function(x, pattern) {
  split_x <- strsplit(x, pattern, perl=TRUE)
  return(split_x[[1]][1])
}

add_zero_one <- function(x) {
  ifelse(x > 0, 1, 0)
}

asv_counts_longer <- function(phyloseq, split_pattern, grep_pattern){
  # Extract data tables 
  otu <- rownames_to_column(as.data.frame(otu_table(phyloseq)@.Data, var = "ASV")) %>% 
    dplyr::rename(ASV = rowname) %>%
    select(-(contains(c('WC', 'BC', 'EB', 'DI', 'Ext', 'FC'))))
  
  # CODE to use LCAs instead of ASVs
  tax <- rownames_to_column(as.data.frame(tax_table(phyloseq)@.Data, var = "ASV")) %>% 
    dplyr::rename(ASV = rowname)
  merged_df <- merge(otu, tax[, c("ASV", "ASV_sequence")], by = "ASV", all.x = TRUE)
  otu$ASV <- merged_df$ASV_sequence
  otu <- otu[!(otu$ASV %in% c("dropped", NA)), ]
  
  # Get the site names from the sample names
  sample_names <- names(otu)
  sample_names <- sample_names[-grep('ASV', sample_names)]
  sites        <- unique(unlist(lapply(sample_names, split_string, pattern=split_pattern)))
  #sites <- sample_names
  
  # Create a new df with the sites instead of samples
  new_otu <- data.frame(matrix(NA, nrow = nrow(otu), ncol = length(sites)))
  colnames(new_otu) <- sites
  for(i in sites){
    otu_sample <- otu[,grep(paste0('^', i, grep_pattern), names(otu))]
    #otu_sample <- data.frame(lapply(otu_sample, add_zero_one)) 
    new_otu[,i] <- rowSums(otu_sample)
    #new_otu[,i] <- otu[,i]
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
  #for (asv in asvs) {
  #  trans_otu[[`asv`]]  <- as.integer(trans_otu[[asv]])
  #}
  for (i in 1:(ncol(trans_otu))) {
    trans_otu[[i]] <- as.integer(trans_otu[[i]])
  }
  
  trans_otu$sample    <- rownames(trans_otu)
  
  # Reformat to long format
  trans_otu %>% 
    pivot_longer(c(asvs)) %>%
    dplyr::rename(asv = 'name',
                  count = 'value') -> otu_long
  return(otu_long)
}



create_heatmap <- function(pooled, independent, pseudo){
  rbind(pooled, independent, pseudo) %>%
    #mutate(asv_short = str_trunc(asv, width = 60, side = "right")) %>%
    ggplot(aes(x=sample, y=asv, fill = log_count)) +
    geom_tile() +
    theme_classic(base_size = 8) +
    facet_wrap(~Option) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), axis.text.y = element_blank()) +
    scale_fill_gradient(low = "#8FBCBB", high = "#B48EAD")
}

remove_outliers <- function(df, cutoff) {
  df <-df[df[["count"]] <= cutoff, ]
  
  return(df)
}

normalize_top_n_asvs <- function(df, n) {
  grouped_df <- df %>%
    group_by(asv) %>%
    summarize(total_count = sum(count))
  
  top_asvs <- grouped_df %>%
    arrange(desc(total_count)) %>%
    head(n) %>%
    select(asv)
  
  top_df <- df %>%
    filter(asv %in% top_asvs$asv)
  #return(top_df)
  top_df_log <- top_df %>%
    mutate(log_count = log(count + 1))
  
  return(top_df_log)
}

normalize_bottom_n_asvs <- function(df, n) {
  grouped_df <- df %>%
    group_by(asv) %>%
    summarize(total_count = sum(count))
  
  bottom_asvs <- grouped_df %>%
    arrange(desc(total_count)) %>%
    tail(n) %>%
    select(asv)
  
  bottom_df <- df %>%
    filter(asv %in% bottom_asvs$asv)
  #return(bottom_df)
  bottom_df_log <- bottom_df %>%
    mutate(log_count = log(count + 1))
  
  return(bottom_df_log)
}

normalize_all_asvs <- function(df) {
  grouped_df <- df %>%
    group_by(asv) %>%
    summarize(total_count = sum(count))
  
  all_asvs <- grouped_df %>%
    arrange(desc(total_count)) %>%
    select(asv)
  
  all_df <- df %>%
    filter(asv %in% all_asvs$asv)
  #return(all_df)
  all_df_log <- all_df %>%
    mutate(log_count = log(count + 1))
  
  return(all_df_log)
}


#==============================================================================
# Import data
#==============================================================================

Cocos_independent <- readRDS('data/phyloseq_objects/trim_fix/CoCosV10I_16S_phyloseq_nt_FALSE_decontaminated.rds')
Cocos_pooled      <- readRDS('data/phyloseq_objects/trim_fix/CoCosV10I_16S_phyloseq_nt_TRUE_decontaminated.rds')
Cocos_pseudo      <- readRDS('data/phyloseq_objects/trim_fix/CoCosV10I_16S_phyloseq_nt_pseudo_decontaminated.rds')

NW_pooled         <- readRDS('data/phyloseq_objects/Pool_TRUE_KWEST_16S_phyloseq_nt_decontaminated.rds')
NW_independent    <- readRDS('data/phyloseq_objects/Pool_FALSE_KWEST_16S_phyloseq_nt_decontaminated.rds')
NW_pseudo         <- readRDS('data/phyloseq_objects/Pool_pseudo_KWEST_16S_phyloseq_nt_decontaminated.rds')
NW_pseudo         <- subset_samples(NW_pseudo, sample_names(NW_pseudo) != "71a")

RS_pooled         <- readRDS('data/phyloseq_objects/trim_fix/RS21AUG_16S_phyloseq_nt_TRUE_decontam.rds')
RS_independent    <- readRDS('data/phyloseq_objects/trim_fix/RS21AUG_16S_phyloseq_nt_FALSE_decontam.rds')
RS_pseudo         <- readRDS('data/phyloseq_objects/trim_fix/RS21AUG_16S_phyloseq_nt_pseudo_decontam.rds')
RS_pooled         <- subset_samples(RS_pooled, sample_names(RS_pooled) != "RS1_ME_S4_1_2_")
RS_independent    <- subset_samples(RS_independent, sample_names(RS_independent) != "RS1_ME_S4_1_2_")
RS_pseudo         <- subset_samples(RS_pseudo, sample_names(RS_pseudo) != "RS1_ME_S4_1_2_")


#==============================================================================
# CococsI
#==============================================================================

pooled      <- asv_counts_longer(Cocos_pooled, CCSPLITPAT, CCGREPPAT) %>% mutate(Option = 'Pooled')
independent <- asv_counts_longer(Cocos_independent, CCSPLITPAT, CCGREPPAT) %>% mutate(Option = 'Independent')
pseudo      <- asv_counts_longer(Cocos_pseudo, CCSPLITPAT, CCGREPPAT) %>% mutate(Option = 'Pseudo')
#pooled      <- remove_outliers(pooled, CCCUTOFF)
#independent <- remove_outliers(independent, CCCUTOFF)
#pseudo      <- remove_outliers(pseudo, CCCUTOFF)
#top_pooled         <- normalize_top_n_asvs(pooled, TOPN)
top_independent    <- normalize_top_n_asvs(independent, TOPN)
#top_pseudo         <- normalize_top_n_asvs(pseudo, TOPN)
#bottom_pooled      <- normalize_bottom_n_asvs(pooled, BOTTOMN)
bottom_independent <- normalize_bottom_n_asvs(independent, BOTTOMN)
#bottom_pseudo      <- normalize_bottom_n_asvs(pseudo, BOTTOMN)
#all_pooled         <- normalize_all_asvs(pooled)
all_independent    <- normalize_all_asvs(independent)
#all_pseudo         <- normalize_all_asvs(pseudo)
top_pooled <- filter(pooled, asv %in% top_independent$asv)
top_pooled         <- normalize_top_n_asvs(top_pooled, TOPN)
top_pseudo <- filter(pseudo, asv %in% top_independent$asv)
top_pseudo         <- normalize_top_n_asvs(top_pseudo, TOPN)
bottom_pooled <- filter(pooled, asv %in% bottom_independent$asv)
bottom_pooled <- normalize_bottom_n_asvs(bottom_pooled, BOTTOMN)
bottom_pseudo <- filter(pseudo, asv %in% bottom_independent$asv)
bottom_pseudo <- normalize_bottom_n_asvs(bottom_pseudo, BOTTOMN)
all_pooled <- filter(pooled, asv %in% all_independent$asv)
all_pooled    <- normalize_all_asvs(all_pooled)
all_pseudo <- filter(pseudo, asv %in% all_independent$asv)
all_pseudo    <- normalize_all_asvs(all_pseudo)


pdf(file="log_counts_asv_top_independant_cocos.pdf", width=WIDTH, height=HEIGHT)
create_heatmap(top_pooled, top_independent, top_pseudo)
dev.off()
pdf(file="log_counts_asv_bottom_independant_cocos.pdf", width=WIDTH, height=HEIGHT)
create_heatmap(bottom_pooled, bottom_independent, bottom_pseudo)
dev.off()
pdf(file="log_counts_asv_all_independant_cocos.pdf", width=WIDTH, height=HEIGHT)
create_heatmap(all_pooled, all_independent, all_pseudo)
dev.off()


#==============================================================================
# North west WA data
#==============================================================================

pooled      <- asv_counts_longer(NW_pooled, NWSPLITPAT, NWGREPPAT) %>% mutate(Option = 'Pooled')
independent <- asv_counts_longer(NW_independent, NWSPLITPAT, NWGREPPAT) %>% mutate(Option = 'Independent')
pseudo      <- asv_counts_longer(NW_pseudo, NWSPLITPAT, NWGREPPAT) %>% mutate(Option = 'Pseudo')
#pooled      <- remove_outliers(pooled, NWCUTOFF)
#independent <- remove_outliers(independent, NWCUTOFF)
#pseudo      <- remove_outliers(pseudo, NWCUTOFF)
#top_pooled         <- normalize_top_n_asvs(pooled, TOPN)
top_independent    <- normalize_top_n_asvs(independent, TOPN)
#top_pseudo         <- normalize_top_n_asvs(pseudo, TOPN)
#bottom_pooled      <- normalize_bottom_n_asvs(pooled, BOTTOMN)
bottom_independent <- normalize_bottom_n_asvs(independent, BOTTOMN)
#bottom_pseudo      <- normalize_bottom_n_asvs(pseudo, BOTTOMN)
#all_pooled         <- normalize_all_asvs(pooled)
all_independent    <- normalize_all_asvs(independent)
#all_pseudo         <- normalize_all_asvs(pseudo)
top_pooled <- filter(pooled, asv %in% top_independent$asv)
top_pooled         <- normalize_top_n_asvs(top_pooled, TOPN)
top_pseudo <- filter(pseudo, asv %in% top_independent$asv)
top_pseudo         <- normalize_top_n_asvs(top_pseudo, TOPN)
bottom_pooled <- filter(pooled, asv %in% bottom_independent$asv)
bottom_pooled <- normalize_bottom_n_asvs(bottom_pooled, BOTTOMN)
bottom_pseudo <- filter(pseudo, asv %in% bottom_independent$asv)
bottom_pseudo <- normalize_bottom_n_asvs(bottom_pseudo, BOTTOMN)
all_pooled <- filter(pooled, asv %in% all_independent$asv)
all_pooled    <- normalize_all_asvs(all_pooled)
all_pseudo <- filter(pseudo, asv %in% all_independent$asv)
all_pseudo    <- normalize_all_asvs(all_pseudo)


pdf(file="log_counts_asv_top_independant_nw.pdf", width=WIDTH, height=HEIGHT)
create_heatmap(top_pooled, top_independent, top_pseudo)
dev.off()
pdf(file="log_counts_asv_bottom_independant_nw.pdf", width=WIDTH, height=HEIGHT)
create_heatmap(bottom_pooled, bottom_independent, bottom_pseudo)
dev.off()
pdf(file="log_counts_asv_all_independant_nw.pdf", width=WIDTH, height=HEIGHT)
create_heatmap(all_pooled, all_independent, all_pseudo)
dev.off()


#==============================================================================
# Rowley Shoals data
#==============================================================================

pooled      <- asv_counts_longer(RS_pooled, RSSPLITPAT, RSGREPPAT) %>% mutate(Option = 'Pooled')
independent <- asv_counts_longer(RS_independent, RSSPLITPAT, RSGREPPAT) %>% mutate(Option = 'Independent')
pseudo      <- asv_counts_longer(RS_pseudo, RSSPLITPAT, RSGREPPAT) %>% mutate(Option = 'Pseudo')
#pooled      <- remove_outliers(pooled, RSCUTOFF)
#independent <- remove_outliers(independent, RSCUTOFF)
#pseudo      <- remove_outliers(pseudo, RSCUTOFF)
#top_pooled         <- normalize_top_n_asvs(pooled, TOPN)
top_independent    <- normalize_top_n_asvs(independent, TOPN)
#top_pseudo         <- normalize_top_n_asvs(pseudo, TOPN)
#bottom_pooled      <- normalize_bottom_n_asvs(pooled, BOTTOMN)
bottom_independent <- normalize_bottom_n_asvs(independent, BOTTOMN)
#bottom_pseudo      <- normalize_bottom_n_asvs(pseudo, BOTTOMN)
#all_pooled         <- normalize_all_asvs(pooled)
all_independent    <- normalize_all_asvs(independent)
#all_pseudo         <- normalize_all_asvs(pseudo)
top_pooled <- filter(pooled, asv %in% top_independent$asv)
top_pooled         <- normalize_top_n_asvs(top_pooled, TOPN)
top_pseudo <- filter(pseudo, asv %in% top_independent$asv)
top_pseudo         <- normalize_top_n_asvs(top_pseudo, TOPN)
bottom_pooled <- filter(pooled, asv %in% bottom_independent$asv)
bottom_pooled <- normalize_bottom_n_asvs(bottom_pooled, BOTTOMN)
bottom_pseudo <- filter(pseudo, asv %in% bottom_independent$asv)
bottom_pseudo <- normalize_bottom_n_asvs(bottom_pseudo, BOTTOMN)
all_pooled <- filter(pooled, asv %in% all_independent$asv)
all_pooled    <- normalize_all_asvs(all_pooled)
all_pseudo <- filter(pseudo, asv %in% all_independent$asv)
all_pseudo    <- normalize_all_asvs(all_pseudo)

pdf(file="log_counts_asv_top_independant_rs.pdf", width=WIDTH, height=HEIGHT)
create_heatmap(top_pooled, top_independent, top_pseudo)
dev.off()
pdf(file="log_counts_asv_bottom_independant_rs.pdf", width=WIDTH, height=HEIGHT)
create_heatmap(bottom_pooled, bottom_independent, bottom_pseudo)
dev.off()
pdf(file="log_counts_asv_all_independant_rs.pdf", width=WIDTH, height=HEIGHT)
create_heatmap(all_pooled, all_independent, all_pseudo)
dev.off()

