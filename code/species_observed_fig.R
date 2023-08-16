## comparison of pooling methods on species

### Figures for paper
library(phyloseq)
library(tidyverse)
library(ggvenn)
library(cowplot)
library(ggsignif)

# CoCos_false  <- readRDS('data/phyloseq_objects/CoCosV10I_16S_phyloseq_nt_FALSE.rds')
# CoCos_true  <- readRDS('data/phyloseq_objects/CoCosV10I_16S_phyloseq_nt_TRUE.rds')
# CoCos_pseudo <- readRDS('data/phyloseq_objects/CoCosV10I_16S_phyloseq_nt_pseudo.rds')
# 
# tax_tab_false <- as.data.frame(CoCos_false@tax_table)
# length(unique((tax_tab_false$species)))
# length(unique((tax_tab_false$LCA)))
# 
# # species_filtered <- tax_tab_false %>%
# #   filter(!(species %in% tax_tab_false$species[grep(NA, (tax_tab_false$species))])) %>%
# #   filter(!(species %in% tax_tab_false$species[grep("dropped", (tax_tab_false$species))]))
# # length(species_filtered$species)
# # 
# # species_filtered <- tax_tab_false %>%
# #   filter(!(genus %in% tax_tab_false$genus[grep(NA, (tax_tab_false$genus))])) %>%
# #   filter(!(genus %in% tax_tab_false$genus[grep("dropped", (tax_tab_false$genus))]))
# # length(species_filtered$genus)
# # 
# 
# tax_tab_pseudo <- as.data.frame(CoCos_pseudo@tax_table)
# length(unique((tax_tab_pseudo$species)))
# length(unique((tax_tab_pseudo$LCA)))
# 
# tax_tab_true <- as.data.frame(CoCos_true@tax_table)
# length(unique((tax_tab_true$species)))
# length(unique((tax_tab_true$LCA)))
# 
# # species_filtered <- tax_tab_true %>%
# #   filter(!(species %in% tax_tab_true$species[grep(NA, (tax_tab_false$species))])) %>%
# #   filter(!(species %in% tax_tab_true$species[grep("dropped", (tax_tab_false$species))]))
# # length(species_filtered$species)
# # 
# # species_filtered <- tax_tab_true %>%
# #   filter(!(genus %in% tax_tab_true$genus[grep(NA, (tax_tab_false$genus))])) %>%
# #   filter(!(genus %in% tax_tab_true$genus[grep("dropped", (tax_tab_false$genus))]))
# # length(species_filtered$genus)
# 
# 
# diff_species <- setdiff(tax_tab_true$species, tax_tab_false$species)
# diff_species
# # [1] "Desmodema polystictum" - heaps in one samples - probable    "Sillago schomburgkii" -  potentially
# # [3] "Ectreposebastes imus" - likely regional but deep water      "Trachinocephalus myops"    - likely
# # [5] "Diretmus argenteus"   - likely      "Magnisudis atlantica"  -one sightings in indian ocean   so potential
# # [7] "Upeneichthys stotti"   - unlikely this species only SWWA but could be close rel    "Bolinichthys distofax"  - global tropical distribution but no recordings near cocos
# # [9] "Notoscopelus caudispinosus" - very rare - likely "Decapterus maruadsi"  - abundant - likely
# # [11] "Kali indica"   - rarer   - possible          "Lethrinus genivittatus"    - rarer - likely
# # [13] "Vinciguerria nimbaria" - very rare   - very likely
# 
# ##### New results
# # [1] "Pempheris mangula"           "Pempheris sp. AWFS-F16-0594"
# 
# rowSums(CoCos_true@otu_table)
# 
# setdiff(unique(tax_tab_true$family), unique(tax_tab_false$family))
# setdiff(unique(tax_tab_true$LCA), unique(tax_tab_false$LCA))
# ##### New results
# # [1] "Pempheris mangula"           "Pempheris sp. AWFS-F16-0594"


### Decontaminated results

# par(mfrow=c(2,2))
# CoCos_false  <- readRDS('data/phyloseq_objects/CoCosV10I_16S_phyloseq_nt_FALSE_decontaminated.rds')
# CoCos_true  <- readRDS('data/phyloseq_objects/CoCosV10I_16S_phyloseq_nt_TRUE_decontaminated.rds')
# CoCos_pseudo <- readRDS('data/phyloseq_objects/CoCosV10I_16S_phyloseq_nt_pseudo_decontaminated.rds')

CoCos_false  <- readRDS('data/phyloseq_objects/trim_fix/CoCosV10I_16S_phyloseq_nt_FALSE_decontaminated.rds')
CoCos_true  <- readRDS('data/phyloseq_objects/trim_fix/CoCosV10I_16S_phyloseq_nt_TRUE_decontaminated.rds')
CoCos_pseudo <- readRDS('data/phyloseq_objects/trim_fix/CoCosV10I_16S_phyloseq_nt_pseudo_decontaminated.rds')


# Remove the controls for Cococs dataets
# pooled_no_controls = subset_samples(CoCos_pseudo, sample_names(CoCos_pseudo) !=
#                                       rownames(sample_data(CoCos_pseudo))[grep("WC", rownames(sample_data(CoCos_pseudo)))])
# samples_to_keep <- rownames(sample_data(CoCos_pseudo))[-grep("WC", rownames(sample_data(CoCos_pseudo)))]
# samples_to_keep <- samples_to_keep[-grep("BC", samples_to_keep)]
# samples_to_keep <- samples_to_keep[-grep("EB", samples_to_keep)]
# samples_to_keep <- samples_to_keep[-grep("DI", samples_to_keep)]
# CoCos_pseudo_nocontrols <- prune_samples(samples_to_keep, CoCos_pseudo)
# 
# samples_to_keep <- rownames(sample_data(CoCos_true))[-grep("WC", rownames(sample_data(CoCos_true)))]
# samples_to_keep <- samples_to_keep[-grep("BC", samples_to_keep)]
# samples_to_keep <- samples_to_keep[-grep("EB", samples_to_keep)]
# samples_to_keep <- samples_to_keep[-grep("DI", samples_to_keep)]
# CoCos_true_nocontrols <- prune_samples(samples_to_keep, CoCos_true)
# 
# samples_to_keep <- rownames(sample_data(CoCos_false))[-grep("WC", rownames(sample_data(CoCos_false)))]
# samples_to_keep <- samples_to_keep[-grep("BC", samples_to_keep)]
# samples_to_keep <- samples_to_keep[-grep("EB", samples_to_keep)]
# samples_to_keep <- samples_to_keep[-grep("DI", samples_to_keep)]
# CoCos_false_nocontrols <- prune_samples(samples_to_keep, CoCos_false)

# Remove the sample that are not in all results
samples_to_keep <- Reduce(intersect, list(rownames(sample_data(CoCos_false)), rownames(sample_data(CoCos_pseudo)), rownames(sample_data(CoCos_true))))

CoCos_false <- prune_samples(samples_to_keep, CoCos_false)
CoCos_true <- prune_samples(samples_to_keep, CoCos_true)
CoCos_pseudo <- prune_samples(samples_to_keep, CoCos_pseudo)

CoCos_false_seqtab <- t(otu_table(CoCos_false))
CoCos_true_seqtab  <- t(otu_table(CoCos_true))
CoCos_pseudo_seqtab <- t(otu_table(CoCos_pseudo))

#--------------------------------------------------------------------------------------------------------------------------------------------------
# ASVs
#--------------------------------------------------------------------------------------------------------------------------------------------------

nsam <- dim(CoCos_true_seqtab)[1]
df.obs <- data.frame(observed=c(rowSums(CoCos_false_seqtab>0), rowSums(CoCos_pseudo_seqtab>0), rowSums(CoCos_true_seqtab>0)),
                     Mode=rep(c("independent", "pseudo", "pooled"), each=nsam),
                     rank=rank(rowSums(CoCos_true_seqtab>0))) #%>%

mode_comp_CocosI_detam <- ggplot(data = df.obs, aes(x = rank, y = observed, color=Mode)) + 
  geom_point(alpha=0.5, size = 3) +
  scale_color_brewer(palette="Dark2") +
  xlab(NULL) + 
  ylab("Observed ASVs") +
  ggtitle("Cocos Islands") +
  theme_bw(base_size = 12) +
  geom_smooth(method = "lm") #+

mode_signif_cocos <- ggplot(df.obs, aes(x = Mode, y = observed, color = Mode)) +
  geom_boxplot() + # using `ggsignif` to display comparison of interest
  geom_signif(
    comparisons = list(c("independent", "pooled")),
    map_signif_level = FALSE, vjust = 3) +
  geom_signif(
    comparisons = list(c("pseudo", "pooled")),
    map_signif_level = FALSE, vjust = 3) +
  # geom_signif(
  #   comparisons = list(c("pseudo", "independent")),
  #   map_signif_level = TRUE) +
  scale_color_brewer(palette="Dark2") +
  xlab(NULL) + 
  ylab("Observed ASVs") +
  ggtitle("Cocos Islands") +
  theme_bw(base_size = 12)
  # scale_color_manual(values = c("#D55E00", "#00FFFF", "#FF00FF"))
  # geom_jitter()


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

#=============================================================
# Look at differences in species
#=============================================================

tax_tab_false <- as.data.frame(CoCos_false@tax_table)
length(unique((tax_tab_false$species)))
length(unique((tax_tab_false$LCA)))

tax_tab_pseudo <- as.data.frame(CoCos_pseudo@tax_table)
length(unique((tax_tab_pseudo$species)))
length(unique((tax_tab_pseudo$LCA)))

tax_tab_true <- as.data.frame(CoCos_true@tax_table)
length(unique((tax_tab_true$species)))
length(unique((tax_tab_true$LCA)))

diff_species <- setdiff(tax_tab_true$species, tax_tab_false$species)
diff_species

### New results
# [1] "Bolinichthys distofax"      "Notoscopelus caudispinosus" "Kali indica"               
# [4] "Lethrinus genivittatus"     "Vinciguerria nimbaria"      "Desmodema polystictum"     
# [7] "Trachinocephalus myops"     "Diretmus argenteus"         "Magnisudis atlantica"      
# [10] "Upeneichthys stotti"        "Kyphosus bigibbus"       

rowSums(CoCos_true@otu_table)

setdiff(tax_tab_false$species, tax_tab_true$species)

setdiff(unique(tax_tab_true$family), unique(tax_tab_false$family))
# [1] "Lethrinidae"    "Phosichthyidae" "Trachipteridae" "Pempheridae"    "Synodontidae"  
# [6] "Diretmidae"     "Paralepididae"  "Istiophoridae" 
setdiff(unique(tax_tab_true$LCA), unique(tax_tab_false$LCA))

#=============================================================
# PLOT species
#=============================================================

CoCos_false_tax  <- reformat_tax_table(CoCos_false, "species")
CoCos_true_tax   <- reformat_tax_table(CoCos_true, "species")
CoCos_pseudo_tax <- reformat_tax_table(CoCos_pseudo, "species")

nsam <- dim(CoCos_true_tax)[1]
df.obs_species <- data.frame(observed=c(rowSums(CoCos_false_tax>0), rowSums(CoCos_pseudo_tax>0), rowSums(CoCos_true_tax>0)),
                             Mode=rep(c("independent", "pseudo", "pooled"), each=nsam),
                             rank=rank(rowSums(CoCos_true_tax>0)), times=4) 

Cocos_species <- ggplot(data = df.obs_species, aes(x = rank, y = observed, color=Mode)) + 
  geom_point(alpha=0.5, size = 3) +
  scale_color_brewer(palette="Dark2") +
  xlab("Samples") + 
  ylab("Observed unique species") +
  # ggtitle("Cocos Islands") +
  theme_bw(base_size = 12) +
  geom_smooth(method = "lm") +
  theme(plot.title = element_text(hjust = 0.5))

mode_signif_cocos_species <- ggplot(df.obs_species, aes(x = Mode, y = observed, color = Mode)) +
  geom_boxplot() + # using `ggsignif` to display comparison of interest
  geom_signif(
    comparisons = list(c("independent", "pooled")),
    map_signif_level = FALSE, vjust = 3) +
  geom_signif(
    comparisons = list(c("pseudo", "pooled")),
    map_signif_level = FALSE, vjust = 3) +
  # geom_signif(
  #   comparisons = list(c("pseudo", "independent")),
  #   map_signif_level = TRUE) +
  scale_color_brewer(palette="Dark2") +
  xlab(NULL) + 
  ylab("Observed unique species") +
  ggtitle("Cocos Islands") +
  theme_bw(base_size = 12)

## add Rowleys data

# RS_false  <- readRDS('data/phyloseq_objects/RS21AUG_16S_phyloseq_nt_FALSE_decontaminated.rds')
# RS_true  <- readRDS('data/phyloseq_objects/RS21AUG_16S_phyloseq_nt_true_decontaminated.rds')
# RS_pseudo <- readRDS('data/phyloseq_objects/RS21AUG_16S_phyloseq_nt_pseudo_decontaminated.rds')

RS_false  <- readRDS('data/phyloseq_objects/trim_fix/RS21AUG_16S_phyloseq_nt_FALSE_decontam.rds')
RS_true  <- readRDS('data/phyloseq_objects/trim_fix/RS21AUG_16S_phyloseq_nt_TRUE_decontam.rds')
RS_pseudo <- readRDS('data/phyloseq_objects/trim_fix/RS21AUG_16S_phyloseq_nt_pseudo_decontam.rds')

# RS_false  <- readRDS('data/phyloseq_objects/RS21AUG_16S_phyloseq_nt_FALSE_pangea_decontaminated.rds')
# RS_true  <- readRDS('data/phyloseq_objects/RS21AUG_16S_phyloseq_nt_TRUE_pangea_decontaminated.rds')
# RS_pseudo  <- readRDS('data/phyloseq_objects/RS21AUG_16S_phyloseq_nt_pseudo_pangea_decontaminated.rds')
#--------------------------------------------------------------------------------------------------------------------------------------------------
# ASVs
#--------------------------------------------------------------------------------------------------------------------------------------------------
RS_false_seqtab <- t(otu_table(RS_false))
RS_true_seqtab  <- t(otu_table(RS_true))
RS_pseudo_seqtab <- t(otu_table(RS_pseudo))


RS_false_seqtab <- RS_false_seqtab[!(row.names(RS_false_seqtab) %in% "RS1_ME_S4_1_2_"),]
RS_true_seqtab <- RS_true_seqtab[!(row.names(RS_true_seqtab) %in% "RS1_ME_S4_1_2_"),]
RS_pseudo_seqtab <- RS_pseudo_seqtab[!(row.names(RS_pseudo_seqtab) %in% "RS1_ME_S4_1_2_"),]


nsam <- dim(RS_true_seqtab)[1]
df.obs <- data.frame(observed=c(rowSums(RS_false_seqtab>0), rowSums(RS_pseudo_seqtab>0), rowSums(RS_true_seqtab>0)),
                     Mode=rep(c("independent", "pseudo", "pooled"), each=nsam),
                     rank=rank(rowSums(RS_true_seqtab>0))) #%>%

mode_comp_RS16S_detam <- ggplot(data = df.obs, aes(x = rank, y = observed, color=Mode)) + 
  geom_point(alpha=0.5, size = 3) +
  scale_color_brewer(palette="Dark2") +
  xlab(NULL) + 
  ylab(NULL) +
  ggtitle("Rowley Shoals") +
  theme_bw(base_size = 12) +
  geom_smooth(method = "lm") #+

mode_signif_RS <- ggplot(df.obs, aes(x = Mode, y = observed, color = Mode)) +
  geom_boxplot() + # using `ggsignif` to display comparison of interest
  geom_signif(
    comparisons = list(c("independent", "pooled")),
    map_signif_level = FALSE, vjust = 3) +
  geom_signif(
    comparisons = list(c("pseudo", "pooled")),
    map_signif_level = FALSE, vjust = 3) +
  # geom_signif(
  #   comparisons = list(c("pseudo", "independent")),
  #   map_signif_level = TRUE) +
  scale_color_brewer(palette="Dark2") +
  xlab(NULL) + 
  ylab(NULL) +
  ggtitle("Rowley Shoals") +
  theme_bw(base_size = 12)

# Remove the controls for Rowley Shoals dataset
# pooled_no_controls = subset_samples(RS_pseudo, sample_names(RS_pseudo) != 
#                                       rownames(sample_data(RS_pseudo))[grep("WC", rownames(sample_data(RS_pseudo)))])
# samples_to_keep <- rownames(sample_data(RS_pseudo))[-grep("WC", rownames(sample_data(RS_pseudo)))]
# samples_to_keep <- samples_to_keep[-grep("BC", samples_to_keep)]
# RS_pseudo_nocontrols <- prune_samples(samples_to_keep, RS_pseudo)
# 
# samples_to_keep <- rownames(sample_data(RS_false))[-grep("WC", rownames(sample_data(RS_false)))]
# samples_to_keep <- samples_to_keep[-grep("BC", samples_to_keep)]
# RS_false_nocontrols <- prune_samples(samples_to_keep, RS_false)
# 
# samples_to_keep <- rownames(sample_data(RS_true))[-grep("WC", rownames(sample_data(RS_true)))]
# samples_to_keep <- samples_to_keep[-grep("BC", samples_to_keep)]
# RS_true_nocontrols <- prune_samples(samples_to_keep, RS_true)

tax_tab_false <- as.data.frame(RS_false@tax_table)
length(unique((tax_tab_false$species))) #
length(unique((tax_tab_false$LCA))) #

tax_tab_pseudo <- as.data.frame(RS_pseudo@tax_table)
length(unique((tax_tab_pseudo$species))) #
length(unique((tax_tab_pseudo$LCA))) #

tax_tab_true <- as.data.frame(RS_true@tax_table)
length(unique((tax_tab_true$species))) #
length(unique((tax_tab_true$LCA))) #

setdiff(tax_tab_true$species, tax_tab_false$species)

rowSums(RS_true@otu_table)

setdiff(tax_tab_false$species, tax_tab_true$species) # 0

rowSums(RS_false@otu_table)

setdiff(unique(tax_tab_true$family), unique(tax_tab_false$family))
# [1] "Bramidae"         "Stomiidae"        "Syngnathidae"     "Triglidae"       
# [5] "Moringuidae"      "Pomacanthidae"    "Gempylidae"       "Nomeidae"        
# [9] "Platycephalidae"  "Mugilidae"        "Carcharhinidae"   "Leiognathidae"   
# [13] "Terapontidae"     "Tetrarogidae"     "Potamotrygonidae" "Sphyraenidae"   
setdiff(unique(tax_tab_false$family), unique(tax_tab_true$family)) #0

#=============================================================
# PLOT species
#=============================================================
# Species
RS_false_tax  <- reformat_tax_table(RS_false, "species")
RS_true_tax   <- reformat_tax_table(RS_true, "species")
RS_pseudo_tax <- reformat_tax_table(RS_pseudo, "species")

# Remove the sample that is not in the site specific analysis
RS_false_tax <- RS_false_tax[!(row.names(RS_false_tax) %in% "RS1_ME_S4_1_2__sum"),]
RS_true_tax <- RS_true_tax[!(row.names(RS_true_tax) %in% "RS1_ME_S4_1_2__sum"),]
RS_pseudo_tax <- RS_pseudo_tax[!(row.names(RS_pseudo_tax) %in% "RS1_ME_S4_1_2__sum"),]
#=============================================================
# PLOT
#=============================================================
nsam <- dim(RS_true_tax)[1]
df.obs_species <- data.frame(observed=c(rowSums(RS_false_tax>0), rowSums(RS_pseudo_tax>0), rowSums(RS_true_tax>0)),
                             Mode=rep(c("independent", "pseudo", "pooled"), each=nsam),
                             rank=rank(rowSums(RS_true_tax>0)), times=4) 

RS16S_species <- ggplot(data = df.obs_species, aes(x = rank, y = observed, color=Mode)) + 
  geom_point(alpha=0.5, size = 3) +
  scale_color_brewer(palette="Dark2") +
  xlab("Samples") + 
  ylab(NULL) +
  # ggtitle("Rowley Shoals") +
  theme_bw(base_size = 12) +
  geom_smooth(method = "lm") +
  theme(plot.title = element_text(hjust = 0.5))

mode_signif_RS_species <- ggplot(df.obs_species, aes(x = Mode, y = observed, color = Mode)) +
  geom_boxplot() + # using `ggsignif` to display comparison of interest
  geom_signif(
    comparisons = list(c("independent", "pooled")),
    map_signif_level = FALSE, vjust = 3) +
  geom_signif(
    comparisons = list(c("pseudo", "pooled")),
    map_signif_level = FALSE, vjust = 3) +
  # geom_signif(
  #   comparisons = list(c("pseudo", "independent")),
  #   map_signif_level = TRUE) +
  scale_color_brewer(palette="Dark2") +
  xlab(NULL) + 
  ylab(NULL) +
  ggtitle(NULL) +
  theme_bw(base_size = 12)


## add NWWA data
NW_false  <- readRDS('data/phyloseq_objects/Pool_FALSE_KWEST_16S_phyloseq_nt_decontaminated.rds')
NW_true  <- readRDS('data/phyloseq_objects/Pool_TRUE_KWEST_16S_phyloseq_nt_decontaminated.rds')
NW_pseudo <- readRDS('data/phyloseq_objects/Pool_pseudo_KWEST_16S_phyloseq_nt_decontaminated.rds')

# NW_false  <- readRDS('data/phyloseq_objects/Pool_FALSE_KWEST_16S_16S_phyloseq_nt.rds')
# NW_true  <- readRDS('data/phyloseq_objects/Pool_TRUE_KWEST_16S_16S_phyloseq_nt.rds')
# NW_pseudo <- readRDS('data/phyloseq_objects/Pool_pseudo_KWEST_16S_16S_phyloseq_nt.rds')

# Remove the controls for Northwest dataset
# pooled_no_controls = subset_samples(NW_pseudo, sample_names(NW_pseudo) != 
#                                       rownames(sample_data(NW_pseudo))[grep("FC", rownames(sample_data(NW_pseudo)))])
# samples_to_keep <- rownames(sample_data(NW_pseudo))[-grep("FC", rownames(sample_data(NW_pseudo)))]
# samples_to_keep <- samples_to_keep[-grep("Extbl", samples_to_keep)]
# NW_pseudo_nocontrols <- prune_samples(samples_to_keep, NW_pseudo)
# 
# samples_to_keep <- rownames(sample_data(NW_false))[-grep("FC", rownames(sample_data(NW_false)))]
# samples_to_keep <- samples_to_keep[-grep("Extbl", samples_to_keep)]
# NW_false_nocontrols <- prune_samples(samples_to_keep, NW_false)
# 
# samples_to_keep <- rownames(sample_data(NW_true))[-grep("FC", rownames(sample_data(NW_true)))]
# samples_to_keep <- samples_to_keep[-grep("Extbl", samples_to_keep)]
# NW_true_nocontrols <- prune_samples(samples_to_keep, NW_true)
tax_tab_false <- as.data.frame(NW_false@tax_table)
length(unique((tax_tab_false$species))) 
length(unique((tax_tab_false$LCA)))

tax_tab_pseudo <- as.data.frame(NW_pseudo@tax_table)
length(unique((tax_tab_pseudo$species)))
length(unique((tax_tab_pseudo$LCA)))

tax_tab_true <- as.data.frame(NW_true@tax_table)
length(unique((tax_tab_true$species)))
length(unique((tax_tab_true$LCA))) 

ind <- setdiff(tax_tab_false$species, tax_tab_true$species)
# [1] "Trachinocephalus myops"         "Decapterus akaadsi"            
# [3] "Platax batavianus"              "Paramonacanthus choirocephalus"
# [5] "Torquigener whitleyi"           "Eleutheronema tetradactylum"   
# [7] "Pomacentrus milleri"            "Ulua mentalis"                 
# [9] "Pseudorhombus jenynsii"
pseudo <- setdiff(tax_tab_pseudo$species, tax_tab_true$species)
# [1] "Platax batavianus"              "Paramonacanthus choirocephalus"
# [3] "Epinephelus corallicola"        "Torquigener whitleyi"          
# [5] "Selar crumenophthalmus"         "Ulua mentalis"                 
# [7] "Pegasus sp. SAA-2015"           "Pseudorhombus jenynsii"        
# [9] "Priacanthus tayenus"  
setdiff(ind, pseudo)
# [1] "Trachinocephalus myops"      "Decapterus akaadsi"         
# [3] "Eleutheronema tetradactylum" "Pomacentrus milleri"  
setdiff(pseudo, ind)
# [1] "Epinephelus corallicola" "Selar crumenophthalmus"  "Pegasus sp. SAA-2015"   
# [4] "Priacanthus tayenus" 

rowSums(NW_false@otu_table)

setdiff(tax_tab_false$family, tax_tab_true$family)
# [1] "Synodontidae"
setdiff(tax_tab_true$family, tax_tab_false$family)

# Remove the sample that are not in all resul
samples_to_keep <- Reduce(intersect, list(rownames(sample_data(NW_false)), rownames(sample_data(NW_pseudo)), rownames(sample_data(NW_true))))

NW_false <- prune_samples(samples_to_keep, NW_false)
NW_true <- prune_samples(samples_to_keep, NW_true)
NW_pseudo <- prune_samples(samples_to_keep, NW_pseudo)

#=============================================================
# ASVs
#=============================================================
NW_false_seqtab <- t(otu_table(NW_false))
NW_true_seqtab  <- t(otu_table(NW_true))
NW_pseudo_seqtab <- t(otu_table(NW_pseudo))

nsam <- dim(NW_true_seqtab)[1]
df.obs <- data.frame(observed=c(rowSums(NW_false_seqtab>0), rowSums(NW_pseudo_seqtab>0), rowSums(NW_true_seqtab>0)),
                     Mode=rep(c("independent", "pseudo", "pooled"), each=nsam),
                     rank=rank(rowSums(NW_true_seqtab>0))) #%>%

mode_comp_NWWA_detam <- ggplot(data = df.obs, aes(x = rank, y = observed, color=Mode)) + 
  geom_point(alpha=0.5, size = 3) +
  scale_color_brewer(palette="Dark2") +
  xlab(NULL) + 
  ylab(NULL) +
  ggtitle("NW WA") +
  theme_bw(base_size = 12) +
  geom_smooth(method = "lm")



mode_signif_NW <- ggplot(df.obs, aes(x = Mode, y = observed, color = Mode)) +
  geom_boxplot() + # using `ggsignif` to display comparison of interest
  geom_signif(
    comparisons = list(c("independent", "pooled")),
    map_signif_level = FALSE, vjust = 3) +
  geom_signif(
    comparisons = list(c("pseudo", "pooled")),
    map_signif_level = FALSE, vjust = 3) +
  # geom_signif(
  #   comparisons = list(c("pseudo", "independent")),
  #   map_signif_level = TRUE) +
  scale_color_brewer(palette="Dark2") +
  xlab(NULL) + 
  ylab(NULL) +
  ggtitle("NW WA") +
  theme_bw(base_size = 12)





# Species
NW_false_tax  <- reformat_tax_table(NW_false, "species")
NW_true_tax   <- reformat_tax_table(NW_true, "species")
NW_pseudo_tax <- reformat_tax_table(NW_pseudo, "species")

#=============================================================
# PLOT
#=============================================================
nsam <- dim(NW_true_tax)[1]
df.obs_species <- data.frame(observed=c(rowSums(NW_false_tax>0), rowSums(NW_pseudo_tax>0), rowSums(NW_true_tax>0)),
                             Mode=rep(c("independent", "pseudo", "pooled"), each=nsam),
                             rank=rank(rowSums(NW_true_tax>0)), times=4) 

NWWA_species <- ggplot(data = df.obs_species, aes(x = rank, y = observed, color=Mode)) + 
  geom_point(alpha=0.5, size = 3) +
  scale_color_brewer(palette="Dark2") +
  xlab("Samples") + 
  ylab(NULL) +
  # ggtitle("North-west") +
  theme_bw(base_size = 12) +
  geom_smooth(method = "lm") +
  theme(plot.title = element_text(hjust = 0.5))


mode_signif_NW_species <- ggplot(df.obs_species, aes(x = Mode, y = observed, color = Mode)) +
  geom_boxplot() + # using `ggsignif` to display comparison of interest
  geom_signif(
    comparisons = list(c("independent", "pooled")),
    map_signif_level = FALSE, vjust = 3) +
  geom_signif(
    comparisons = list(c("pseudo", "pooled")),
    map_signif_level = FALSE, vjust = 3) +
  # geom_signif(
  #   comparisons = list(c("pseudo", "independent")),
  #   map_signif_level = TRUE) +
  scale_color_brewer(palette="Dark2") +
  xlab(NULL) + 
  ylab(NULL) +
  ggtitle(NULL) +
  theme_bw(base_size = 12)



# # Create plot grid
# plot <- plot_grid(Cocos_species + theme(legend.position="none"),
#                   RS16S_species + theme(legend.position="none"),
#                   NWWA_species + theme(legend.position="none"),
#                   Cocos_family + theme(legend.position="bottom") + labs(fill = "Replicates") + ylab("Families") + xlab("Sites") + theme(axis.ticks = element_blank()),
#                   RS_family + theme(legend.position="bottom") + labs(fill = "Replicates") + xlab("Sites") + theme(axis.ticks = element_blank()),
#                   NW_family + theme(legend.position="bottom") + labs(fill = "Replicates") + xlab("Sites") + theme(axis.ticks = element_blank()),
#                   # nrow = 1,
#                   labels = 'AUTO',
#                   label_size = 14,
#                   hjust = -1.1,
#                   vjust = 2)
# 
# plot

# plot <- plot_grid(Cocos_species + theme(legend.position="none"),
#                   RS16S_species + theme(legend.position="none"),
#                   NWWA_species + theme(legend.position="none"),
#                   nrow = 1,
#                   labels = 'AUTO',
#                   label_size = 14,
#                   hjust = -1.1,
#                   vjust = 2)
# 
# plot
# 
# # add in legend
# legend <- get_legend(
#   RS16S_species +
#     guides(color = guide_legend(nrow = 1)) +
#     theme(legend.position = "bottom") +
#     theme(legend.text = element_text(size = 12))
# )

# plot_1 <- plot_grid(plot, legend, ncol = 1, rel_heights = c(1, .1))
# plot_2 <- plot_grid(Cocos_family + theme(legend.position="bottom") + labs(fill = "Replicates") + ylab("Families") + xlab("Sites") + theme(axis.ticks = element_blank()),
#                   RS_family + theme(legend.position="bottom") + labs(fill = "Replicates") + xlab("Sites") + theme(axis.ticks = element_blank()),
#                   NW_family + theme(legend.position="bottom") + labs(fill = "Replicates") + xlab("Sites") + theme(axis.ticks = element_blank()),
#                   nrow = 1,
#                   labels = c('D', 'E', 'F'),
#                   label_size = 14,
#                   hjust = -1.1,
#                   vjust = 2)
# 
# plot_2
# plot_grid(plot_1, plot_2, ncol = 1, rel_heights = c(0.5, 1))

# Create plot grid
plot <- plot_grid(mode_comp_CocosI_detam + theme(legend.position="none") + theme(axis.title=element_text(size=12)),
                         mode_comp_RS16S_detam + theme(legend.position="none"),
                         mode_comp_NWWA_detam + theme(legend.position="none"),
                         Cocos_species + theme(legend.position="none") + theme(axis.title=element_text(size=12)),
                         RS16S_species + theme(legend.position="none"),
                         NWWA_species + theme(legend.position="none"),
                         # nrow = 1,
                         labels = 'AUTO',
                         label_size = 14)
                         # hjust = -1.1,
                         # vjust = 2)

# add in legend
legend <- get_legend(
  mode_comp_RS16S_detam +
    guides(color = guide_legend(nrow = 1)) +
    theme(legend.position = "bottom")
)

plot_grid(plot, legend, ncol = 1, rel_heights = c(1, .1))





# Create plot grid
plot <- plot_grid(mode_signif_cocos + theme(legend.position="none") + theme(axis.title=element_text(size=12)) + theme(axis.text.x = element_blank()),
                  mode_signif_RS + theme(legend.position="none")  + theme(axis.text.x = element_blank()),
                  mode_signif_NW + theme(legend.position="none")  + theme(axis.text.x = element_blank()),
                  mode_signif_cocos_species + theme(legend.position="none") + theme(axis.title=element_text(size=12)) + ggtitle(NULL)  + theme(axis.text.x = element_blank()),
                  mode_signif_RS_species + theme(legend.position="none")  + theme(axis.text.x = element_blank()),
                  mode_signif_NW_species + theme(legend.position="none")  + theme(axis.text.x = element_blank()),
                  # nrow = 1,
                  labels = 'AUTO',
                  label_size = 14)
# hjust = -1.1,
# vjust = 2)

# add in legend
legend <- get_legend(
  mode_signif_cocos +
    guides(color = guide_legend(nrow = 1)) +
    theme(legend.position = "bottom")
)

plot_grid(plot, legend, ncol = 1, rel_heights = c(1, .1))



