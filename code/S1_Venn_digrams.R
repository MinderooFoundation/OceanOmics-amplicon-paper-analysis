#===================================================================================================
# Venn Diagram: Overlap of identified species

# Load libraries
library(ggvenn)
library(cowplot)

individual  <- readRDS('data/phyloseq_objects/CoCosV10I_16S_phyloseq_nt_FALSE_decontaminated.rds')
pooled  <- readRDS('data/phyloseq_objects/CoCosV10I_16S_phyloseq_nt_TRUE_decontaminated.rds')
pseudo <- readRDS('data/phyloseq_objects/CoCosV10I_16S_phyloseq_nt_pseudo_decontaminated.rds')


taxa_pooled     <- unique(as.character(tax_table(pooled)@.Data[,"species"]))
taxa_individual <- unique(as.character(tax_table(individual)@.Data[,"species"]))
taxa_pseudo     <- unique(as.character(tax_table(pseudo)@.Data[,"species"]))


list_venn <- list(`Independent` = taxa_individual,
                  `Pooled` = taxa_pooled,
                  `Pseudo` = taxa_pseudo)

VENN_cocos <- ggvenn(list_venn, c("Independent", "Pooled", "Pseudo"),
               fill_color =c('#F8766D', '#00BA38','#619CFF'))
  # labs(title = "Overlap of unique species",
       # subtitle = "Comparison of the three different dada2 'pool=' modes",
       # caption = "Data source: CoCos Island transect (16S), December 2022, Minderoo Foundation OceanOmics")
VENN_cocos


#===================================================================================================
# Rowley Shoals data


individual  <- readRDS('data/phyloseq_objects/RS21AUG_16S_phyloseq_nt_FALSE_decontam.rds')
pooled      <- readRDS('data/phyloseq_objects/RS21AUG_16S_phyloseq_nt_TRUE_decontam.rds')
pseudo      <- readRDS('data/phyloseq_objects/RS21AUG_16S_phyloseq_nt_pseudo_decontam.rds')

taxa_pooled     <- unique(as.character(tax_table(pooled)@.Data[,"species"]))
taxa_individual <- unique(as.character(tax_table(individual)@.Data[,"species"]))
taxa_pseudo     <- unique(as.character(tax_table(pseudo)@.Data[,"species"]))


list_venn <- list(`Independent` = taxa_individual,
                  `Pooled` = taxa_pooled,
                  `Pseudo` = taxa_pseudo)

VENN_RS <- ggvenn(list_venn, c("Independent", "Pooled", "Pseudo"),
               fill_color =c('#F8766D', '#00BA38','#619CFF')) 
  # labs(title = "Overlap of unique species",
  #      subtitle = "Comparison of the three different dada2 'pool=' modes",
  #      caption = "Data source: Rowley Shoals (16S), August 2021, Minderoo Foundation OceanOmics")
VENN_RS

#===================================================================================================
# North-West WA data


individual  <- readRDS('data/phyloseq_objects/NWWA_16S_phyloseq_nt_FALSE_decontam.rds')
pooled      <- readRDS('data/phyloseq_objects/NWWA_16S_phyloseq_nt_TRUE_decontam.rds')
pseudo      <- readRDS('data/phyloseq_objects/NWWA_16S_phyloseq_nt_pseudo_decontam.rds')

taxa_pooled     <- unique(as.character(tax_table(pooled)@.Data[,"species"]))
taxa_individual <- unique(as.character(tax_table(individual)@.Data[,"species"]))
taxa_pseudo     <- unique(as.character(tax_table(pseudo)@.Data[,"species"]))


list_venn <- list(`Independent` = taxa_individual,
                  `Pooled` = taxa_pooled,
                  `Pseudo` = taxa_pseudo)

VENN_NW <- ggvenn(list_venn, c("Independent", "Pooled", "Pseudo"),
               fill_color =c('#F8766D', '#00BA38','#619CFF')) 
  # labs(title = "Overlap of unique species",
  #      subtitle = "Comparison of the three different dada2 'pool=' modes",
  #      caption = "Data source: North West WA (16S), West K, et al.")
VENN_NW


plot_grid(VENN_cocos,
          VENN_RS,
          VENN_NW,
          labels = 'AUTO', nrow = 1, label_size = 20)

