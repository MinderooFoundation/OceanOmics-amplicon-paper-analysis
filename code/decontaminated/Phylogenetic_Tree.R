# Creat phylogenetic circle plot
# https://yulab-smu.top/treedata-book/chapter10.html

library(tidyverse)
library(ggtree)
library(rotl)
library(rfishbase)
library(ggnewscale)
library(phyloseq)
library(ggtreeExtra) # Bioconductor

RS_false  <- readRDS('data/phyloseq_objects_decontaminated/RSV5_16S_FALSE_phyloseq_nt_decontaminated.rds')
RS_true  <- readRDS('data/phyloseq_objects_decontaminated/RSV5_16S_TRUE_phyloseq_nt_decontaminated.rds')
RS_pseudo <- readRDS('data/phyloseq_objects_decontaminated/RSV5_16S_pseudo_phyloseq_nt_decontaminated.rds')

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
RS_false_tax  <- reformat_tax_table(RS_false, "species")
RS_true_tax   <- reformat_tax_table(RS_true, "species")
RS_pseudo_tax <- reformat_tax_table(RS_pseudo, "species")


# Remove the sample that is not in the site specific analysis
RS_false_tax <- RS_false_tax[!(row.names(RS_false_tax) %in% "RS1_ME_S4_1(2)_sum"),]
RS_true_tax <- RS_true_tax[!(row.names(RS_true_tax) %in% "RS1_ME_S4_1(2)_sum"),]
RS_pseudo_tax <- RS_pseudo_tax[!(row.names(RS_pseudo_tax) %in% "RS1_ME_S4_1(2)_sum"),]


# The whole section below on preparing the data can be moved out of this script?
if(!file.exists('data/fishbase/fishbase_species.RDS')) {
  fishbase_species <- fb_tbl('species')
  fishbase_families <- fb_tbl('families')
  fishbase_synonyms <- fb_tbl('synonyms')
  saveRDS(object = fishbase_species, file = 'data/phylo_tree/fishbase_species.RDS')
  saveRDS(object = fishbase_families, file='data/phylo_tree/fishbase_families.RDS')
  saveRDS(object = fishbase_synonyms, file='data/phylo_tree/fishbase_synonyms.RDS')
} else {
  fishbase_species <- readRDS('data/phylo_tree/fishbase_species.RDS')
  fishbase_families <- readRDS('data/phylo_tree/fishbase/fishbase_families.RDS')
  fishbase_synonyms <- readRDS('data/phylo_tree/fishbase/fishbase_synonyms.RDS')
}

#=============================================================================================
# Data cleanup

synomyms_lookup <- fishbase_synonyms %>% 
  unite('SynSpecies',  SynGenus:SynSpecies, sep= ' ') %>% 
  filter(Status != 'misapplied name') %>% 
  left_join(fishbase_species, by='SpecCode') %>% 
  select(Genus, Species, SynSpecies) %>%
  unite('Species', Genus:Species, sep=' ')

synomyms_lookup <- synomyms_lookup %>%
  filter(Species != 'NA NA')

  fishbase_species_to_family <- fishbase_species |> 
  left_join(fishbase_families |> 
              select(FamCode, Family)) |> 
  select(Genus, Species, Family) |> 
  unite(col = 'Species', Genus:Species, sep = ' ')

# GBIF
# Read in the GBIF data, get the fishbase entries and filter the GBIF data by fishbase 
gbif             <- read_tsv('data/phylo_tree/0198426-220831081235567.csv')

fishbase_species <- fishbase_species %>% 
  unite('Species2', Genus:Species, sep = ' ')

gbif <- gbif %>% 
  filter(species %in% fishbase_species$Species2)

# Remove Elasmobranchii and hagfishes (there's two sightings of a hagfish in GBIF for RS!)

gbif <- gbif %>% #filter( (class != 'Elasmobranchii') %>% replace_na(TRUE)) %>% 
  filter( (!order %in% c('Chimaeriformes', 'Myxiniformes')) %>% replace_na(TRUE))

#=============================================================================================
# TREE
# Reading in the guide tree that Philipp prepared
#tr <- readRDS('data/phylo_tree/phyloreturns.RDS')
tr <- treeio::read.newick('data/phylo_tree/betancour.tre')
# Extract the taxa table and reformat it so that the rows are samples and columns are taxa
get_taxa_names <- function(phyloseq_object, taxa_column = c('domain', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'LCA')){
  
  taxa_table <- merge(as.data.frame(tax_table(phyloseq_object)), as.data.frame(otu_table(phyloseq_object)),
                      by = 'row.names', all = TRUE) %>%
    group_by(!!sym(taxa_column)) %>% 
    summarize_at(vars(names(as.data.frame(otu_table(phyloseq_object)))), list(sum=sum)) -> taxa_table
  
  taxa_table             <- as.data.frame((taxa_table)) 
  rownames(taxa_table)   <- taxa_table[,taxa_column]
  taxa_table[,taxa_column] <- NULL
  taxa_table             <- t(taxa_table)
  unique(colnames(taxa_table))
  
}


# Species level
TREE_TAX <- unique(unlist(lapply(sub("_[^_]+$", "", unlist(lapply(tr$tip.label, function(x) sub("^[^_]*_", "", x)))), function(x) str_replace_all(x, '_', ' '))))





annot_table <- as.data.frame(tibble(
  Node = TREE_TAX, 
  RS_pseudo = TREE_TAX %in% unique(colnames(RS_pseudo_tax)), 
  RS_false = TREE_TAX %in% unique(colnames(RS_false_tax)),
  RS_true = TREE_TAX %in% unique(colnames(RS_true_tax)),
  GBIF = TREE_TAX %in% unique(gbif$species)))

row.names(annot_table) <- annot_table$Node
annot_table$Node <- NULL

counts_of_falses <- rowSums(annot_table == FALSE)
bad_names <- names(counts_of_falses[counts_of_falses==3])
good_tr <- ape::drop.tip(phy = tr, tip = bad_names)


all_taxa <- rownames(annot_table[which(rowSums(annot_table) > 0), ])



tr$tip.label <- unlist(lapply(sub("_[^_]+$", "", unlist(lapply(tr$tip.label, function(x) sub("^[^_]*_", "", x)))), function(x) str_replace_all(x, '_', ' ')))
good_tr <- ape::drop.tip(phy = tr, tip = setdiff(rownames(annot_table), all_taxa))


circ <- ggtree(good_tr, layout = "circular")
circ2 <- rotate_tree(circ, 90)

# The tree lables were still the "old ones" with "_" and the code at the end. so replace them here as well
circ2$data <- circ2$data %>%
  mutate(label = (unlist(lapply(sub("_[^_]+$", "", unlist(lapply(circ2$data$label, function(x) sub("^[^_]*_", "", x)))), function(x) str_replace_all(x, '_', ' ')))))

good_tr %>%
  coord_polar()
# Plot the heatmap tree
p1 <- 
gheatmap(
  circ2,
  annot_table["RS_false"],
  width = .1,
  offset = 0.4,
  custom_column_labels = "Individual", #c('Pseudo', 'Individual', 'Pooled', 'GBIF'),
  font.size = 4,
  color = 'black'
) +
  scale_fill_manual(values = c("grey", "#1B9E77")) + 
  theme(legend.position = 'none')

circ2 +
  geom_tile(data = annot_table, aes(x, y, fill = value), 
                      width = width, color = color, inherit.aes = FALSE)


gheatmap(
  circ2,
  annot_table[c("RS_true", "RS_false", "GBIF")],
  width = .3,
  offset = 0.4,
  custom_column_labels = c('Pooled', 'Individual', 'GBIF'),
  font.size = 4,
  color = 'black'
) +
  scale_fill_manual(values = c("grey", "#1B9E77")) + 
  theme(legend.position = 'none') +
  geom_cladelab()

circ2 + geom_tile(fill = annot_table["RS_true"])

