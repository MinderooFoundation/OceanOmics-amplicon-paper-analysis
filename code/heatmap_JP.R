# ASV count heatmaps for paper

#==============================================================================
# Setup
#==============================================================================

library(tidyverse)
library(tidyr)
library(phyloseq)
library(cowplot)
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

asv_counts_longer <- function(phyloseq, split_pattern, grep_pattern){
  # Extract data tables 
  otu <- rownames_to_column(as.data.frame(otu_table(phyloseq)@.Data, var = "ASV")) %>% 
    dplyr::rename(ASV = rowname) %>%
    select(-(contains(c('WC', 'BC', 'EB', 'DI', 'Ext', 'FC'))))
  
  # Get the site names from the sample names
  sample_names <- names(otu)
  sample_names <- sample_names[-grep('ASV', sample_names)]
  sites        <- unique(unlist(lapply(sample_names, split_string, pattern=split_pattern)))
  
  # Create a new df with the sites instead of samples
  new_otu <- data.frame(matrix(NA, nrow = nrow(otu), ncol = length(sites)))
  colnames(new_otu) <- sites
  for(i in sites){
    otu_sample <- otu[,grep(paste0('^', i, grep_pattern), names(otu))]
    new_otu[,i] <- rowSums(otu_sample)
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
  for (asv in asvs) {
    trans_otu[[asv]]  <- as.integer(trans_otu[[asv]])
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
    ggplot(aes(x=sample, y=asv, fill = log_count)) +
    geom_tile() +
    theme_classic(base_size = 12) +
    facet_wrap(~Option) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
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
  
  all_df_log <- all_df %>%
    mutate(log_count = log(count + 1))
  
  return(all_df_log)
}


#==============================================================================
# Import data
#==============================================================================

Cocos_indepe <- readRDS('data/phyloseq_objects/trim_fix/CoCosV10I_16S_phyloseq_nt_FALSE_decontaminated.rds')
Cocos_pooled      <- readRDS('data/phyloseq_objects/trim_fix/CoCosV10I_16S_phyloseq_nt_TRUE_decontaminated.rds')
Cocos_pseudo      <- readRDS('data/phyloseq_objects/trim_fix/CoCosV10I_16S_phyloseq_nt_pseudo_decontaminated.rds')

NW_pooled         <- readRDS('data/phyloseq_objects/Pool_TRUE_KWEST_16S_phyloseq_nt_decontaminated.rds')
NW_indepe    <- readRDS('data/phyloseq_objects/Pool_FALSE_KWEST_16S_phyloseq_nt_decontaminated.rds')
NW_pseudo         <- readRDS('data/phyloseq_objects/Pool_pseudo_KWEST_16S_phyloseq_nt_decontaminated.rds')
NW_pseudo         <- subset_samples(NW_pseudo, sample_names(NW_pseudo) != "71a")

RS_pooled         <- readRDS('data/phyloseq_objects/trim_fix/RS21AUG_16S_phyloseq_nt_TRUE_decontam.rds')
RS_inde    <- readRDS('data/phyloseq_objects/trim_fix/RS21AUG_16S_phyloseq_nt_FALSE_decontam.rds')
RS_pseudo         <- readRDS('data/phyloseq_objects/trim_fix/RS21AUG_16S_phyloseq_nt_pseudo_decontam.rds')
RS_pooled         <- subset_samples(RS_pooled, sample_names(RS_pooled) != "RS1_ME_S4_1_2_")
RS_indepe    <- subset_samples(RS_independent, sample_names(RS_independent) != "RS1_ME_S4_1_2_")
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
top_pooled         <- normalize_top_n_asvs(pooled, TOPN)
top_independent    <- normalize_top_n_asvs(independent, TOPN)
top_pseudo         <- normalize_top_n_asvs(pseudo, TOPN)
bottom_pooled      <- normalize_bottom_n_asvs(pooled, BOTTOMN)
bottom_independent <- normalize_bottom_n_asvs(independent, BOTTOMN)
bottom_pseudo      <- normalize_bottom_n_asvs(pseudo, BOTTOMN)
all_pooled         <- normalize_all_asvs(pooled)
all_independent    <- normalize_all_asvs(independent)
all_pseudo         <- normalize_all_asvs(pseudo)

# pdf(file="top_asv_counts_cocos.pdf", width=WIDTH, height=HEIGHT)
cocos_top <- create_heatmap(top_pooled, top_independent, top_pseudo)
# dev.off()
# pdf(file="bottom_asv_counts_cocos.pdf", width=WIDTH, height=HEIGHT)
cocos_bottom <- create_heatmap(bottom_pooled, bottom_independent, bottom_pseudo)
# dev.off()
# pdf(file="all_asv_counts_cocos.pdf", width=WIDTH, height=HEIGHT)
# create_heatmap(all_pooled, all_independent, all_pseudo)
# dev.off()

#==============================================================================
# North west WA data
#==============================================================================

pooled      <- asv_counts_longer(NW_pooled, NWSPLITPAT, NWGREPPAT) %>% mutate(Option = 'Pooled')
independent <- asv_counts_longer(NW_independent, NWSPLITPAT, NWGREPPAT) %>% mutate(Option = 'Independent')
pseudo      <- asv_counts_longer(NW_pseudo, NWSPLITPAT, NWGREPPAT) %>% mutate(Option = 'Pseudo')
#pooled      <- remove_outliers(pooled, NWCUTOFF)
#independent <- remove_outliers(independent, NWCUTOFF)
#pseudo      <- remove_outliers(pseudo, NWCUTOFF)
top_pooled         <- normalize_top_n_asvs(pooled, TOPN)
top_independent    <- normalize_top_n_asvs(independent, TOPN)
top_pseudo         <- normalize_top_n_asvs(pseudo, TOPN)
bottom_pooled      <- normalize_bottom_n_asvs(pooled, BOTTOMN)
bottom_independent <- normalize_bottom_n_asvs(independent, BOTTOMN)
bottom_pseudo      <- normalize_bottom_n_asvs(pseudo, BOTTOMN)
all_pooled         <- normalize_all_asvs(pooled)
all_independent    <- normalize_all_asvs(independent)
all_pseudo         <- normalize_all_asvs(pseudo)

# pdf(file="top_asv_counts_nw.pdf", width=WIDTH, height=HEIGHT)
nw_top <- create_heatmap(top_pooled, top_independent, top_pseudo)
# dev.off()
# pdf(file="bottom_asv_counts_nw.pdf", width=WIDTH, height=HEIGHT)
nw_bottom <- create_heatmap(bottom_pooled, bottom_independent, bottom_pseudo)
# dev.off()
# pdf(file="all_asv_counts_nw.pdf", width=WIDTH, height=HEIGHT)
# create_heatmap(all_pooled, all_independent, all_pseudo)
# dev.off()


#==============================================================================
# Rowley Shoals data
#==============================================================================

pooled      <- asv_counts_longer(RS_pooled, RSSPLITPAT, RSGREPPAT) %>% mutate(Option = 'Pooled')
independent <- asv_counts_longer(RS_independent, RSSPLITPAT, RSGREPPAT) %>% mutate(Option = 'Independent')
pseudo      <- asv_counts_longer(RS_pseudo, RSSPLITPAT, RSGREPPAT) %>% mutate(Option = 'Pseudo')
#pooled      <- remove_outliers(pooled, RSCUTOFF)
#independent <- remove_outliers(independent, RSCUTOFF)
#pseudo      <- remove_outliers(pseudo, RSCUTOFF)
top_pooled         <- normalize_top_n_asvs(pooled, TOPN)
top_independent    <- normalize_top_n_asvs(independent, TOPN)
top_pseudo         <- normalize_top_n_asvs(pseudo, TOPN)
bottom_pooled      <- normalize_bottom_n_asvs(pooled, BOTTOMN)
bottom_independent <- normalize_bottom_n_asvs(independent, BOTTOMN)
bottom_pseudo      <- normalize_bottom_n_asvs(pseudo, BOTTOMN)
all_pooled         <- normalize_all_asvs(pooled)
all_independent    <- normalize_all_asvs(independent)
all_pseudo         <- normalize_all_asvs(pseudo)

# pdf(file="top_asv_counts_rs.pdf", width=WIDTH, height=HEIGHT)
rs_top <- create_heatmap(top_pooled, top_independent, top_pseudo)
# dev.off()
# pdf(file="bottom_asv_counts_rs.pdf", width=WIDTH, height=HEIGHT)
rs_bottom <- create_heatmap(bottom_pooled, bottom_independent, bottom_pseudo)
# dev.off()
# pdf(file="all_asv_counts_rs.pdf", width=WIDTH, height=HEIGHT)
# create_heatmap(all_pooled, all_independent, all_pseudo)
# dev.off()

plot_grid(cocos_top + ggtitle("Cocos Islands") + ylab("Top 20 ASVs") + theme(axis.text.x = element_blank()) + theme(axis.text.y = element_blank()),
          cocos_bottom + ylab("Bottom 20 ASVs") + theme(axis.text.x = element_blank()) + theme(axis.text.y = element_blank()),
          ncol = 1)

plot_grid(rs_top + ggtitle("Rowley Shoals") + ylab("Top 20 ASVs") + theme(axis.text.x = element_blank()) + theme(axis.text.y = element_blank()),
          rs_bottom + ylab("Bottom 20 ASVs") + theme(axis.text.x = element_blank()) + theme(axis.text.y = element_blank()),
          ncol = 1)

plot_grid(nw_top + ggtitle("NW WA") + ylab("Top 20 ASVs") + theme(axis.text.x = element_blank()) + theme(axis.text.y = element_blank()),
          nw_bottom + ylab("Bottom 20 ASVs") + theme(axis.text.x = element_blank()) + theme(axis.text.y = element_blank()),
          ncol = 1)
