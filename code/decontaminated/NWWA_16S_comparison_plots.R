# This is a test script to check if we can do the pool comparison with
# a phyloseq object only

library(phyloseq)
library(tidyverse)
library(ggpubr)


RS_false  <- readRDS('data/phyloseq_objects_decontaminated/Pool_FALSE_KWEST_16S_phyloseq_nt_decontaminated.rds')
RS_true  <- readRDS('data/phyloseq_objects_decontaminated/Pool_TRUE_KWEST_16S_phyloseq_nt_decontaminated.rds')
RS_pseudo <- readRDS('data/phyloseq_objects_decontaminated/Pool_pseudo_KWEST_16S_phyloseq_nt_decontaminated.rds')


RS_false_seqtab <- t(otu_table(RS_false))
RS_true_seqtab  <- t(otu_table(RS_true))
RS_pseudo_seqtab <- t(otu_table(RS_pseudo))



nsam <- dim(RS_true_seqtab)[1]
df.obs <- data.frame(observed=c(rowSums(RS_false_seqtab>0), rowSums(RS_pseudo_seqtab>0), rowSums(RS_true_seqtab>0)),
                     Mode=rep(c("independent", "pseudo", "pooled"), each=nsam),
                     rank=rank(rowSums(RS_true_seqtab>0))) #%>%

mode_comp <- ggplot(data = df.obs, aes(x = rank, y = observed, color=Mode)) + 
  geom_point(alpha=0.5, size = 3) +
  xlab("Samples") + 
  ylab("Observed ASVs") +
  ggtitle("Number of ASVs per sample") +
  theme_bw(base_size = 18) +
  geom_smooth(method = "lm") #+


#--------------------------------------------------------------------------------------------------------------------------------------------------
# TAXA
#--------------------------------------------------------------------------------------------------------------------------------------------------

# Extract the taxa table and reformat it so that the rows are samples and columns are taxa
reformat_tax_table <- function(phyloseq_object, taxa_column = c('domain', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'LCA')){
  
  taxa_table <- merge(as.data.frame(tax_table(phyloseq_object)), as.data.frame(otu_table(phyloseq_object)),
                      by = 'row.names', all = TRUE) %>%
    group_by(!!sym(taxa_column)) %>% 
    summarize_at(vars(names(as.data.frame(otu_table(phyloseq_object)))), list(sum=sum)) -> taxa_table
  
  taxa_table             <- as.data.frame((taxa_table)) 
  rownames(taxa_table)   <- taxa_table[,taxa_column]
  taxa_table[,taxa_column] <- NULL
  taxa_table             <- t(taxa_table)
  taxa_table
  
}

# LCA
RS_false_tax  <- reformat_tax_table(RS_false, "LCA")
RS_true_tax   <- reformat_tax_table(RS_true, "LCA")
RS_pseudo_tax <- reformat_tax_table(RS_pseudo, "LCA")



#=============================================================
# PLOT
#=============================================================
df.obs_lca <- data.frame(observed=c(rowSums(RS_false_tax>0), rowSums(RS_pseudo_tax>0), rowSums(RS_true_tax>0)),
                         Mode=rep(c("independent", "pseudo", "pooled"), each=nsam),
                         rank=rank(rowSums(RS_true_tax>0))) 

mode_comp_taxa <- ggplot(data = df.obs_lca, aes(x = rank, y = observed, color=Mode)) + 
  geom_point(alpha=0.5, size = 3) +
  xlab("Samples") + 
  ylab("Observed unique LCA") +
  ggtitle("Number of unique LCA level taxa per sample") +
  theme_bw(base_size = 18) +
  geom_smooth(method = "lm") 


# Species
RS_false_tax  <- reformat_tax_table(RS_false, "species")
RS_true_tax   <- reformat_tax_table(RS_true, "species")
RS_pseudo_tax <- reformat_tax_table(RS_pseudo, "species")

#=============================================================
# PLOT
#=============================================================
df.obs_species <- data.frame(observed=c(rowSums(RS_false_tax>0), rowSums(RS_pseudo_tax>0), rowSums(RS_true_tax>0)),
                             Mode=rep(c("independent", "pseudo", "pooled"), each=nsam),
                             rank=rank(rowSums(RS_true_tax>0)), times=4) 

mode_comp_species <- ggplot(data = df.obs_species, aes(x = rank, y = observed, color=Mode)) + 
  geom_point(alpha=0.5, size = 3) +
  xlab("Samples") + 
  ylab("Observed unique species") +
  ggtitle("Number of unique species level taxa per sample") +
  theme_bw(base_size = 18) +
  geom_smooth(method = "lm") #+


ggarrange(mode_comp,mode_comp_taxa, mode_comp_species, ncol = 3 , common.legend = TRUE, legend = "top")


#===================================================================================================
# Venn Diagram


individual  <- readRDS('data/phyloseq_objects_decontaminated/Pool_FALSE_KWEST_16S_16S_phyloseq_nt_decontaminated.rds')
pooled      <- readRDS('data/phyloseq_objects_decontaminated/Pool_TRUE_KWEST_16S_16S_phyloseq_nt_decontaminated.rds')
pseudo      <- readRDS('data/phyloseq_objects_decontaminated/Pool_pseudo_KWEST_16S_16S_phyloseq_nt_decontaminated.rds')

taxa_pooled     <- unique(as.character(tax_table(pooled)@.Data[,"species"]))
taxa_individual <- unique(as.character(tax_table(individual)@.Data[,"species"]))
taxa_pseudo     <- unique(as.character(tax_table(pseudo)@.Data[,"species"]))


list_venn <- list(`Independent Analysis` = taxa_individual,
                  `Pooled Analysis` = taxa_pooled,
                  `Pseudo Analysis` = taxa_pseudo)

VENN <- ggvenn(list_venn, c("Independent Analysis", "Pooled Analysis", "Pseudo Analysis"),
       fill_color =c('#F8766D', '#00BA38','#619CFF')) +
  labs(title = "Overlap of unique species",
       subtitle = "Comparison of the three different dada2 'pool=' modes",
       caption = "Data source: North West WA (16S), West K, et al.")









