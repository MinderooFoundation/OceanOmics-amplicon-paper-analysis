# Recreate Denise's plots
# install.packages("remotes")   ## run this line if you do not already have remotes installed
# library(remotes)
# remotes::install_github("adw96/breakaway")
# remotes::install_github("adw96/DivNet")
# library(DivNet)
# # if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# # BiocManager::install(c("phyloseq", "microbiome", "ComplexHeatmap"), update = FALSE)
# # install.packages(
# #   "microViz",
# #   repos = c(davidbarnett = "https://david-barnett.r-universe.dev", getOption("repos"))
# # )
# library(microViz)
# # if (!require("BiocManager", quietly = TRUE))
# #   install.packages("BiocManager")
# # BiocManager::install("decontam")
# library(decontam)
# library(zCompositions)
# # devtools::install_github("tpq/propr")
# library("propr")
# BiocManager::install("ALDEx2")
# BiocManager::install("EnhancedVolcano")

# Load libraries
library("repmis")
LoadandCite(pkgs=c("repmis", "knitr", "tinytex", "phyloseq", "ggplot2", 
                   "zCompositions", "propr", "compositions", "ggfortify", 
                   "ALDEx2", "EnhancedVolcano", "plyr", "microViz", "decontam", 
                   "patchwork", "readxl", "DivNet", "dplyr"), file = "report/packages.bib")




## Load data and create metadata info for this dataset
samples <- list()
samples$west_ind <- readRDS(file = "data/phyloseq_objects/Pool_FALSE_KWEST_16S_16S_phyloseq_nt.rds")
samples$west_pooled <- readRDS(file = "data/phyloseq_objects/Pool_TRUE_KWEST_16S_16S_phyloseq_nt.rds")
samples$west_pseudo <- readRDS(file = "data/phyloseq_objects/Pool_pseudo_KWEST_16S_16S_phyloseq_nt.rds")

## Create bioplot function
pca_biplot <- function(pca_data, phyloseq, colour, colour_label, shape, shape_label, title) {
  long_arrows <- names(which(abs(pca_data$rotation[ ,"PC1"]) > 0.1 | abs(pca_data$rotation[ ,"PC2"]) > 0.1))
  loadings_label <- tax_table(phyloseq)@.Data[ , "LCA_ASV"]
  loadings_label[!names(loadings_label) %in% long_arrows] <- ""
  
  biplot <- autoplot(pca_data,
                     data = sample_data(phyloseq),
                     colour = colour, shape = shape,
                     label = TRUE, label.label = sample_data(phyloseq)$Sample,
                     label.size = 4, label.vjust = -0.5,
                     loadings = TRUE, loadings.colour = "gray80", loadings.label = TRUE,
                     loadings.label.label = unname(loadings_label),
                     loadings.label.size = 3, loadings.label.colour = "gray30", loadings.label.repel = TRUE) + theme_bw() + guides(colour = guide_legend(title = colour_label), shape = guide_legend(title = shape_label)) + theme(legend.title = element_text(size = 12), legend.text = element_text(size = 11)) + ggtitle(label = title)
  
  return(biplot)
}

pca_biplot_unlabelled <- function(pca_data, phyloseq, size, colour, colour_label, shape, shape_label, title) {
  long_arrows <- names(which(abs(pca_data$rotation[ ,"PC1"]) > 0.1 | abs(pca_data$rotation[ ,"PC2"]) > 0.1))
  loadings_label <- tax_table(phyloseq)@.Data[ , "LCA_ASV"]
  loadings_label[!names(loadings_label) %in% long_arrows] <- ""
  
  biplot <- autoplot(pca_data,
                     data = sample_data(phyloseq),
                     colour = colour, size = size, shape = shape,
                     loadings = TRUE, loadings.colour = "gray80", loadings.label = TRUE,
                     loadings.label.label = unname(loadings_label),
                     loadings.label.size = 3, loadings.label.colour = "gray30", loadings.label.repel = TRUE) + theme_bw() + guides(colour = guide_legend(title = colour_label), shape = guide_legend(title = shape_label)) + theme(legend.title = element_text(size = 12), legend.text = element_text(size = 11)) + ggtitle(label = title)
  
  return(biplot)
}


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
samples <- mapply(names(samples), FUN = function(x) tax_mutate(samples[[x]], Contamination = as.character(samples_contam[[x]]$contaminant)))
samples <- lapply(samples, FUN = function(x) subset_taxa(x, Contamination == FALSE))
samples <- lapply(samples, FUN = function(x) subset_samples(x, Control != TRUE))
samples <- lapply(samples, FUN = function(y) prune_samples(samples = sample_sums(y) >= 1000, x = y))

## Transform zeros
samples_zimp <- lapply(samples, FUN = function(x) cmultRepl(t(otu_table(x)), method = "CZM", output = "p-counts"))

## Transform data with a centred log ratio (CLR)
samples_clr <- lapply(samples_zimp, FUN = function(x) propr(counts = x))

## Perform PCA on the CLR transformed data and plot on covariance biplots
samples_pca <- lapply(samples_clr, FUN = function(x) prcomp(x@logratio))

## Biplots
### Independent
NW_pca_independent <- plot(pca_biplot_unlabelled(pca_data = samples_pca$west_ind, phyloseq = samples$west_ind, size = 3, colour = "Subregion", colour_label = "Sub region", shape = "Habitat", shape_label = "Habitat", title = "Independent"))

### Pooled
NW_pca_pooled <- plot(pca_biplot_unlabelled(pca_data = samples_pca$west_pooled, phyloseq = samples$west_pooled, size = 3, colour = "Subregion", colour_label = "Sub region", shape = "Habitat", shape_label = "Habitat", title = "Pooled"))

### Pseudo
NW_pca_pseudo <- plot(pca_biplot_unlabelled(pca_data = samples_pca$west_pseudo, phyloseq = samples$west_pseudo, size = 3, colour = "Subregion", colour_label = "Sub region", shape = "Habitat", shape_label = "Habitat", title = "Pseudo"))

### Create figure
plot1 <- plot_grid(NW_pca_independent + theme(legend.position="none"), 
                  # NW_pca_pseudo + theme(legend.position="none"),
                  NW_pca_pooled + theme(legend.position="none"), 
                  nrow = 1)

# add in legend
legend <- get_legend(
  NW_pca_independent +
    guides(color = guide_legend(nrow = 1)) +
    theme(legend.position = "bottom") +
    theme(legend.text = element_text(size = 12))
)

plot1 <- plot_grid(plot1, legend, rel_heights = c(1, .1), ncol = 1)
plot1







## Divnet analysis
samples <- lapply(samples, FUN = function(x) {sample_data(x)$Site <- paste(sample_data(x)$Subregion, sample_data(x)$Habitat, sample_data(x)$site, sep = "_"); return(x)})

divnet_alpha_NW <- lapply(samples, FUN = function(x) divnet(W = x, X = "Site", base = names(sort(taxa_sums(x))[round(length(taxa_sums(x))/2)])))
saveRDS(divnet_alpha_NW, file = "data/divnet_alpha_NW.rds")
# divnet_alpha <- readRDS(file = "data/divnet_alpha_NW.rds")

### by habitat
NW_alpha_independent <- data.frame(divnet_alpha$west_ind$shannon %>% summary, sample_data(samples$west_ind)) %>%
  distinct(estimate, error, lower, upper, Site, Habitat) %>%
  ggplot(aes(x = Site, y = estimate, col = Habitat)) + geom_point() + geom_segment(aes(x = .data[["Site"]], xend = .data[["Site"]], y = .data[["lower"]], yend = .data[["upper"]])) + xlab("Site") + ylab("Shannon entropy estimate") + theme_bw() + theme(axis.text.x = element_text(angle=90, hjust = 1))

NW_alpha_pooled <- data.frame(divnet_alpha$west_pooled$shannon %>% summary, sample_data(samples$west_pooled)) %>%
  distinct(estimate, error, lower, upper, Site, Habitat) %>%
  ggplot(aes(x = Site, y = estimate, col = Habitat)) + geom_point() + geom_segment(aes(x = .data[["Site"]], xend = .data[["Site"]], y = .data[["lower"]], yend = .data[["upper"]])) + xlab("Site") + ylab("Shannon entropy estimate") + theme_bw() + theme(axis.text.x = element_text(angle=90, hjust = 1))

NW_alpha_pseudo <- data.frame(divnet_alpha$west_pseudo$shannon %>% summary, sample_data(samples$west_pseudo)) %>%
  distinct(estimate, error, lower, upper, Site, Habitat) %>%
  ggplot(aes(x = Site, y = estimate, col = Habitat)) + geom_point() + geom_segment(aes(x = .data[["Site"]], xend = .data[["Site"]], y = .data[["lower"]], yend = .data[["upper"]])) + xlab("Site") + ylab("Shannon entropy estimate") + theme_bw() + theme(axis.text.x = element_text(angle=90, hjust = 1))

### by subregion
NW_alpha_independent <- data.frame(divnet_alpha$west_ind$shannon %>% summary, sample_data(samples$west_ind)) %>%
  distinct(estimate, error, lower, upper, Site, Subregion) %>%
  ggplot(aes(x = Site, y = estimate, col = Subregion)) + geom_point() + geom_segment(aes(x = .data[["Site"]], xend = .data[["Site"]], y = .data[["lower"]], yend = .data[["upper"]])) + xlab("Site") + ylab("Shannon entropy estimate") + theme_bw() + theme(axis.text.x = element_text(angle=90, hjust = 1))

NW_alpha_pooled <- data.frame(divnet_alpha$west_pooled$shannon %>% summary, sample_data(samples$west_pooled)) %>%
  distinct(estimate, error, lower, upper, Site, Subregion) %>%
  ggplot(aes(x = Site, y = estimate, col = Subregion)) + geom_point() + geom_segment(aes(x = .data[["Site"]], xend = .data[["Site"]], y = .data[["lower"]], yend = .data[["upper"]])) + xlab("Site") + ylab("Shannon entropy estimate") + theme_bw() + theme(axis.text.x = element_text(angle=90, hjust = 1))

NW_alpha_pseudo <- data.frame(divnet_alpha$west_pseudo$shannon %>% summary, sample_data(samples$west_pseudo)) %>%
  distinct(estimate, error, lower, upper, Site, Subregiont) %>%
  ggplot(aes(x = Site, y = estimate, col = Subregion)) + geom_point() + geom_segment(aes(x = .data[["Site"]], xend = .data[["Site"]], y = .data[["lower"]], yend = .data[["upper"]])) + xlab("Site") + ylab("Shannon entropy estimate") + theme_bw() + theme(axis.text.x = element_text(angle=90, hjust = 1))

### Create figure
plot <- plot_grid(NW_alpha_independent + theme(legend.position="none") + theme(axis.text.x = element_blank()), 
                  NW_alpha_pooled + theme(legend.position="none") + ylab(NULL)  + theme(axis.text.x = element_blank()), 
                  nrow = 1)
plot
# add in legend
legend <- get_legend(
  NW_alpha_independent +
    guides(color = guide_legend(nrow = 1)) +
    theme(legend.position = "bottom") +
    theme(legend.text = element_text(size = 12))
)

NW_alpha_subregion <- plot_grid(plot, legend, rel_heights = c(1, .1), ncol = 1)

plot <- plot_grid(NW_alpha_independent + theme(legend.position="none"), 
                  NW_alpha_pooled + theme(legend.position="none") + ylab(NULL), 
                  nrow = 1)
plot

# ### Create figure
# plot <- plot_grid(NW_pca_independent + theme(legend.position="none"), 
#                   NW_pca_pooled + theme(legend.position="none"),
#                   NW_alpha_independent + theme(legend.position="none") + ggtitle(NULL) + xlab(NULL), 
#                   NW_alpha_pooled + theme(legend.position="none") + ylab(NULL) + ggtitle(NULL) + xlab(NULL))
# # nrow = 1)
# plot

# add in legend
legend <- get_legend(
  NW_alpha_independent +
    guides(color = guide_legend(nrow = 1)) +
    theme(legend.position = "bottom") +
    theme(legend.text = element_text(size = 12))
)

NW_alpha_habitat <- plot_grid(plot, legend, rel_heights = c(1, .1), ncol = 1)

alphas <- plot_grid(NW_alpha_subregion, NW_alpha_habitat, ncol = 1, rel_heights = c(1.1, 2))

plot_grid(plot1, alphas, ncol = 1)

