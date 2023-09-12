#-----------------------------------------------------------------------------------------------------------------
# Decontamination of all data: this takes the phyloseq objects and 

# Load libraries

library(cowplot)
library(readxl)

LoadandCite(pkgs=c("repmis", "knitr", "tinytex", "phyloseq", "ggplot2", 
                   "zCompositions", "propr", "compositions", "ggfortify", 
                   "ALDEx2", "EnhancedVolcano", "plyr", "microViz", "decontam", 
                   "patchwork", "readxl", "DivNet", "dplyr", "cowplot"))

## Cocos Islands

## Load data
metadata           <- data.frame(read_excel(path = "data/metadata/V10_CKI_eDNA_metadata_P1.xlsx"))
samples            <- list()
samples$cki_ind    <- readRDS(file = "data/phyloseq_objects/CoCosV10I_16S_phyloseq_nt_FALSE.rds")
samples$cki_pooled <- readRDS(file = "data/phyloseq_objects/CoCosV10I_16S_phyloseq_nt_TRUE.rds")
samples$cki_pseudo <- readRDS(file = "data/phyloseq_objects/CoCosV10I_16S_phyloseq_nt_pseudo.rds")

length(rownames(samples$cki_ind@otu_table))
length(rownames(samples$cki_pseudo@otu_table))
length(rownames(samples$cki_pooled@otu_table))

## Format and filter data for ASVs that occur in more than 3 samples/replicates
samples <- lapply(samples, FUN = function(y) {sample_data(y)$Sample <- gsub(pattern = "V10_CKI_", replacement = "", x = rownames(sample_data(y))); return(y)})
samples <- lapply(samples, FUN = function(x) tax_mutate(x, LCA_ASV = paste0(unname(tax_table(x)@.Data[, "LCA"]), " (", rownames(tax_table(x)), ")")))
samples <- lapply(samples, FUN = function(x) subset_samples(x, !(replicate_id %in% c("BC", "DI"))))
samples <- lapply(samples, FUN = function(x) {sample_data(x)$Control <- ifelse(sample_data(x)$replicate_id %in% c("EB", "WC"), TRUE, FALSE); return(x)})
samples <- lapply(samples, FUN = function(x) filter_taxa(x, flist = function(y) sum(y >= 1) >= 3, prune = TRUE))
samples <- lapply(samples, FUN = function(x) {colnames(sample_data(x))[2] <- "sampling_method"; return(x)})
samples <- lapply(samples, FUN = function(x) {sample_data(x)$EEZ <- metadata$EEZ[match(rownames(sample_data(x)), metadata$sample_id)]; return(x)})
samples <- lapply(samples, FUN = function(x) {sample_data(x)$Project <- metadata$Project[match(rownames(sample_data(x)), metadata$sample_id)]; return(x)})
samples <- lapply(samples, FUN = function(x) {sample_data(x)$Group <- metadata$Group[match(rownames(sample_data(x)), metadata$sample_id)]; return(x)})

length(rownames(samples$cki_ind@otu_table)) 
length(rownames(samples$cki_pseudo@otu_table)) 
length(rownames(samples$cki_pooled@otu_table)) 

## Decontamination of data
samples_contam <- lapply(samples, FUN = function(x) isContaminant(x, method = "prevalence", neg = "Control", threshold = 0.2))
samples_contam <- lapply(samples_contam, FUN = function(x) {x$Prevalence <- cut(x$prev, breaks = c(3, 5, 10, max(x$prev)), labels = c(">=3 and <=5", ">5 and <=10", ">10"), include.lowest = TRUE); return(x)})
samples <- mapply(names(samples), FUN = function(x) tax_mutate(samples[[x]], Contamination = as.character(samples_contam[[x]]$contaminant)))

p1 <- ggplot(samples_contam$cki_ind, aes(x = p, fill = Prevalence)) + geom_histogram(colour = "black", breaks = seq(from = 0, to = 1, by = 0.05), alpha = 0.2) + geom_vline(xintercept = 0.2, colour = "red") + xlab("Probability") + ylab("Count") + scale_y_continuous(expand = c(0, 0)) + scale_x_continuous(expand = c(0, 0))
p2 <- ggplot(samples_contam$cki_pooled, aes(x = p, fill = Prevalence)) + geom_histogram(colour = "black", breaks = seq(from = 0, to = 1, by = 0.05), alpha = 0.2) + geom_vline(xintercept = 0.2, colour = "red") + xlab("Probability") + ylab("Count") + scale_y_continuous(expand = c(0, 0)) + scale_x_continuous(expand = c(0, 0))
p3 <- ggplot(samples_contam$cki_pseudo, aes(x = p, fill = Prevalence)) + geom_histogram(colour = "black", breaks = seq(from = 0, to = 1, by = 0.05), alpha = 0.2) + geom_vline(xintercept = 0.2, colour = "red") + xlab("Probability") + ylab("Count") + scale_y_continuous(expand = c(0, 0)) + scale_x_continuous(expand = c(0, 0))

p1 + p2 + p3 + plot_layout(guides = "collect", ncol = 2) + plot_annotation(tag_levels = "A")
rm(p1, p2, p3)

# Prevalence plots looks the same with 0.2 threshold
samples <- lapply(samples, FUN = function(x) subset_taxa(x, Contamination == FALSE))
samples <- lapply(samples, FUN = function(x) subset_samples(x, !(replicate_id %in% c("EB", "WC"))))
samples <- lapply(samples, FUN = function(x) subset_samples(x, !is.na(EEZ)))

samples <- lapply(samples, FUN = function(y) prune_samples(samples = sample_sums(y) >= 1000, x = y))
samples <- lapply(samples, FUN = function(x) {sample_data(x)$site_id <- factor(as.numeric(sample_data(x)$site_id)); return(x)})

length(rownames(samples$cki_ind@otu_table)) 
length(rownames(samples$cki_pseudo@otu_table)) 
length(rownames(samples$cki_pooled@otu_table)) 

saveRDS(samples$cki_ind, file = "data/phyloseq_objects/CoCosV10I_16S_phyloseq_nt_FALSE_decontaminated.rds")
saveRDS(samples$cki_pooled, file = "data/phyloseq_objects/CoCosV10I_16S_phyloseq_nt_TRUE_decontaminated.rds")
saveRDS(samples$cki_pseudo, file = "data/phyloseq_objects/CoCosV10I_16S_phyloseq_nt_pseudo_decontaminated.rds")




## Rowley Shoals

## Load data
metadata <- read.csv(file = "data/metadata/RSV5_metadata.csv")

samples <- list()
samples$rs_ind <- readRDS(file = "data/phyloseq_objects/RS21AUG_16S_phyloseq_nt_FALSE.rds")
samples$rs_pooled <- readRDS(file = "data/phyloseq_objects/RS21AUG_16S_phyloseq_nt_TRUE.rds")
samples$rs_pseudo <- readRDS(file = "data/phyloseq_objects/RS21AUG_16S_phyloseq_nt_pseudo.rds")

length(rownames(samples$rs_ind@otu_table)) 
length(rownames(samples$rs_pseudo@otu_table))
length(rownames(samples$rs_pooled@otu_table))

## Format and filter data for ASVs that occur in more than 3 samples/replicates
samples <- lapply(samples, FUN = function(x) {sample_data(x)$Sample <- rownames(sample_data(x)); return(x)})
samples <- lapply(samples, FUN = function(x) tax_mutate(x, LCA_ASV = paste0(unname(tax_table(x)@.Data[, "LCA"]), " (", rownames(tax_table(x)), ")")))
samples <- lapply(samples, FUN = function(x) subset_samples(x, Replicate.ID != "BC"))
samples <- lapply(samples, FUN = function(x) {sample_data(x)$Control <- ifelse(sample_data(x)$Replicate.ID == "WC", TRUE, FALSE); return(x)})
samples <- lapply(samples, FUN = function(x) filter_taxa(x, flist = function(y) sum(y >= 1) >= 3, prune = TRUE))

length(rownames(samples$rs_ind@otu_table)) 
length(rownames(samples$rs_pseudo@otu_table))
length(rownames(samples$rs_pooled@otu_table))

## Decontamination of data
samples_contam <- lapply(samples, FUN = function(x) isContaminant(x, method = "prevalence", neg = "Control", threshold = 0.5))
samples_contam <- lapply(samples_contam, FUN = function(x) {x$Prevalence <- cut(x$prev, breaks = c(3, 5, 10, max(x$prev)), labels = c(">=3 and <=5", ">5 and <=10", ">10"), include.lowest = TRUE); return(x)})

# Prevalence plots 
p1 <- ggplot(samples_contam$rs_ind, aes(x = p, fill = Prevalence)) + geom_histogram(colour = "black", breaks = seq(from = 0, to = 1, by = 0.05), alpha = 0.5) + geom_vline(xintercept = 0.5, colour = "red") + xlab("Probability") + ylab("Count") + scale_y_continuous(expand = c(0, 0)) + scale_x_continuous(expand = c(0, 0))
p2 <- ggplot(samples_contam$rs_pooled, aes(x = p, fill = Prevalence)) + geom_histogram(colour = "black", breaks = seq(from = 0, to = 1, by = 0.05), alpha = 0.5) + geom_vline(xintercept = 0.5, colour = "red") + xlab("Probability") + ylab("Count") + scale_y_continuous(expand = c(0, 0)) + scale_x_continuous(expand = c(0, 0))
p3 <- ggplot(samples_contam$rs_pseudo, aes(x = p, fill = Prevalence)) + geom_histogram(colour = "black", breaks = seq(from = 0, to = 1, by = 0.05), alpha = 0.5) + geom_vline(xintercept = 0.5, colour = "red") + xlab("Probability") + ylab("Count") + scale_y_continuous(expand = c(0, 0)) + scale_x_continuous(expand = c(0, 0))

p1 + p2 + p3 + plot_layout(guides = "collect", ncol = 2) + plot_annotation(tag_levels = "A")
rm(p1, p2, p3)

samples <- mapply(names(samples), FUN = function(x) tax_mutate(samples[[x]], Contamination = as.character(samples_contam[[x]]$contaminant)))
samples <- lapply(samples, FUN = function(x) subset_taxa(x, Contamination == FALSE))
samples <- lapply(samples, FUN = function(x) subset_samples(x, Replicate.ID != "WC"))
samples <- lapply(samples, FUN = function(x) subset_samples(x, Replicate.ID != "NTC"))
samples <- lapply(samples, FUN = function(y) prune_samples(samples = sample_sums(y) >= 1000, x = y))

length(rownames(samples$rs_ind@otu_table)) 
length(rownames(samples$rs_pseudo@otu_table))
length(rownames(samples$rs_pooled@otu_table))

saveRDS(samples$rs_ind, file = "data/phyloseq_objects/RS21AUG_16S_phyloseq_nt_FALSE_decontam.rds")
saveRDS(samples$rs_pooled, file = "data/phyloseq_objects/RS21AUG_16S_phyloseq_nt_TRUE_decontam.rds")
saveRDS(samples$rs_pseudo, file = "data/phyloseq_objects/RS21AUG_16S_phyloseq_nt_pseudo_decontam.rds")



## North-west Western Australia (NW WA)

## Load data and create metadata info for this dataset
samples <- list()
samples$west_ind <- readRDS(file = "data/phyloseq_objects/NWWA_16S_phyloseq_nt_FALSE.rds")
samples$west_pooled <- readRDS(file = "data/phyloseq_objects/NWWA_16S_phyloseq_nt_TRUE.rds")
samples$west_pseudo <- readRDS(file = "data/phyloseq_objects/NWWA_16S_phyloseq_nt_pseudo.rds")

## Format and filter data for ASVs that occur in more than 3 samples/replicates
samples <- lapply(samples, FUN = function(x) {sample_data(x)$Sample <- rownames(sample_data(x)); return(x)})
samples <- lapply(samples, FUN = function(x) tax_mutate(x, LCA_ASV = paste0(unname(tax_table(x)@.Data[, "LCA"]), " (", rownames(tax_table(x)), ")")))
samples <- lapply(samples, FUN = function(y) {sample_data(y)$Control <- grepl(pattern = "Extbl|FC", x = sample_data(y)$Sample); return(y)})
samples <- lapply(samples, FUN = function(x) filter_taxa(x, flist = function(y) sum(y >= 1) >= 3, prune = TRUE))
samples <- lapply(samples, FUN = function(x) {sample_data(x)$Bioregion <- ifelse(sample_data(x)$site %in% 1:7, "Canning", "Kimberley"); return(x)})
samples <- lapply(samples, FUN = function(x) {sample_data(x)$Subregion <- ifelse(sample_data(x)$site %in% 1:7, "Dampier Peninsula", 
                                                                                 ifelse(sample_data(x)$site %in% c(8:44, 67:71), "South Kimberley", "North Kimberley")); return(x)})
samples <- lapply(samples, FUN = function(x) {sample_data(x)$Habitat <- ifelse(sample_data(x)$site %in% c(1:12, 26:29, 34:37, 42:48, 52, 55:63), "Inshore",
                                                                               ifelse(sample_data(x)$site %in% c(13:17, 21:25, 49:51, 53:54, 64:65), "Coastal", "Nearshore Estuarine")); return(x)})
samples <- lapply(samples, FUN = function(x) {sample_data(x)$Habitat[sample_data(x)$site == 66] <- "Midshelf"; return(x)})

## Decontamination of data
samples_contam <- lapply(samples, FUN = function(x) isContaminant(x, method = "prevalence", neg = "Control", threshold = 0.3))
samples_contam <- lapply(samples_contam, FUN = function(x) {x$Prevalence <- cut(x$prev, breaks = c(3, 5, 10, max(x$prev)), labels = c(">=3 and <=5", ">5 and <=10", ">10"), include.lowest = TRUE); return(x)})

# Prevalence plots 
p1 <- ggplot(samples_contam$west_ind, aes(x = p, fill = Prevalence)) + geom_histogram(colour = "black", breaks = seq(from = 0, to = 1, by = 0.05), alpha = 0.5) + geom_vline(xintercept = 0.5, colour = "red") + xlab("Probability") + ylab("Count") + scale_y_continuous(expand = c(0, 0)) + scale_x_continuous(expand = c(0, 0))
p2 <- ggplot(samples_contam$west_pooled, aes(x = p, fill = Prevalence)) + geom_histogram(colour = "black", breaks = seq(from = 0, to = 1, by = 0.05), alpha = 0.5) + geom_vline(xintercept = 0.5, colour = "red") + xlab("Probability") + ylab("Count") + scale_y_continuous(expand = c(0, 0)) + scale_x_continuous(expand = c(0, 0))
p3 <- ggplot(samples_contam$west_pseudo, aes(x = p, fill = Prevalence)) + geom_histogram(colour = "black", breaks = seq(from = 0, to = 1, by = 0.05), alpha = 0.5) + geom_vline(xintercept = 0.5, colour = "red") + xlab("Probability") + ylab("Count") + scale_y_continuous(expand = c(0, 0)) + scale_x_continuous(expand = c(0, 0))

p1 + p2 + p3 + plot_layout(guides = "collect", ncol = 2) + plot_annotation(tag_levels = "A")
rm(p1, p2, p3)

samples <- mapply(names(samples), FUN = function(x) tax_mutate(samples[[x]], Contamination = as.character(samples_contam[[x]]$contaminant)))
samples <- lapply(samples, FUN = function(x) subset_taxa(x, Contamination == FALSE))
samples <- lapply(samples, FUN = function(x) subset_samples(x, Control != TRUE))
samples <- lapply(samples, FUN = function(y) prune_samples(samples = sample_sums(y) >= 1000, x = y))

length(rownames(samples$west_ind@otu_table)) 
length(rownames(samples$west_pseudo@otu_table))
length(rownames(samples$west_pooled@otu_table))

saveRDS(samples$west_ind, file = "data/phyloseq_objects/NWWA_16S_phyloseq_nt_FALSE_decontam.rds")
saveRDS(samples$west_pooled, file = "data/phyloseq_objects/NWWA_16S_phyloseq_nt_TRUE_decontam.rds")
saveRDS(samples$west_pseudo, file = "data/phyloseq_objects/NWWA_16S_phyloseq_nt_pseudo_decontam.rds")










