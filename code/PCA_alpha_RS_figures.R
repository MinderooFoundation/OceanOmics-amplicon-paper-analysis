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
library(tidyverse)
library(ggvenn)
library(cowplot)
library("repmis")
LoadandCite(pkgs=c("repmis", "knitr", "tinytex", "phyloseq", "ggplot2", 
                   "zCompositions", "propr", "compositions", "ggfortify", 
                   "ALDEx2", "EnhancedVolcano", "plyr", "microViz", "decontam", 
                   "patchwork", "readxl", "DivNet", "dplyr"), file = "report/packages.bib")

## Load data
metadata <- read.csv(file = "data/metadata/RSV5_metadata.csv")

samples <- list()
# samples$rs_ind <- readRDS(file = "data/phyloseq_objects/RS21AUG_16S_phyloseq_nt_FALSE.rds")
# samples$rs_pooled <- readRDS(file = "data/phyloseq_objects/RS21AUG_16S_phyloseq_nt_TRUE.rds")
# samples$rs_pseudo <- readRDS(file = "data/phyloseq_objects/RS21AUG_16S_phyloseq_nt_pseudo.rds")
# samples$rs_ind <- readRDS(file = "data/phyloseq_objects/RS21AUG_16S_phyloseq_nt_FALSE_pangaea.rds")
# samples$rs_pooled <- readRDS(file = "data/phyloseq_objects/RS21AUG_16S_phyloseq_nt_TRUE_pangaea.rds")
# samples$rs_pseudo <- readRDS(file = "data/phyloseq_objects/RS21AUG_16S_phyloseq_nt_pseudo_pangaea.rds")

samples$rs_ind <- readRDS(file = "data/phyloseq_objects/trim_fix/RS21AUG_16S_phyloseq_nt_FALSE.rds")
samples$rs_pooled <- readRDS(file = "data/phyloseq_objects/trim_fix/RS21AUG_16S_phyloseq_nt_TRUE.rds")
samples$rs_pseudo <- readRDS(file = "data/phyloseq_objects/trim_fix/RS21AUG_16S_phyloseq_nt_pseudo.rds")


length(rownames(samples$rs_ind@otu_table)) #1167
length(rownames(samples$rs_pseudo@otu_table)) #1284
length(rownames(samples$rs_pooled@otu_table)) #1256

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
samples <- lapply(samples, FUN = function(x) subset_samples(x, Replicate.ID != "BC"))
samples <- lapply(samples, FUN = function(x) {sample_data(x)$Control <- ifelse(sample_data(x)$Replicate.ID == "WC", TRUE, FALSE); return(x)})
samples <- lapply(samples, FUN = function(x) filter_taxa(x, flist = function(y) sum(y >= 1) >= 3, prune = TRUE))

length(rownames(samples$rs_ind@otu_table)) #625
length(rownames(samples$rs_pseudo@otu_table)) #857
length(rownames(samples$rs_pooled@otu_table)) #1166

## Decontamination of data
samples_contam <- lapply(samples, FUN = function(x) isContaminant(x, method = "prevalence", neg = "Control", threshold = 0.5))
samples_contam <- lapply(samples_contam, FUN = function(x) {x$Prevalence <- cut(x$prev, breaks = c(3, 5, 10, max(x$prev)), labels = c(">=3 and <=5", ">5 and <=10", ">10"), include.lowest = TRUE); return(x)})

# Prevalence plots 
# p1 <- ggplot(samples_contam$rs_ind, aes(x = p, fill = Prevalence)) + geom_histogram(colour = "black", breaks = seq(from = 0, to = 1, by = 0.05), alpha = 0.5) + geom_vline(xintercept = 0.5, colour = "red") + xlab("Probability") + ylab("Count") + scale_y_continuous(expand = c(0, 0)) + scale_x_continuous(expand = c(0, 0))
# p2 <- ggplot(samples_contam$rs_pooled, aes(x = p, fill = Prevalence)) + geom_histogram(colour = "black", breaks = seq(from = 0, to = 1, by = 0.05), alpha = 0.5) + geom_vline(xintercept = 0.5, colour = "red") + xlab("Probability") + ylab("Count") + scale_y_continuous(expand = c(0, 0)) + scale_x_continuous(expand = c(0, 0))
# p3 <- ggplot(samples_contam$rs_pseudo, aes(x = p, fill = Prevalence)) + geom_histogram(colour = "black", breaks = seq(from = 0, to = 1, by = 0.05), alpha = 0.5) + geom_vline(xintercept = 0.5, colour = "red") + xlab("Probability") + ylab("Count") + scale_y_continuous(expand = c(0, 0)) + scale_x_continuous(expand = c(0, 0))
# 
# p1 + p2 + p3 + plot_layout(guides = "collect", ncol = 2) + plot_annotation(tag_levels = "A")
# rm(p1, p2, p3)

samples <- mapply(names(samples), FUN = function(x) tax_mutate(samples[[x]], Contamination = as.character(samples_contam[[x]]$contaminant)))
samples <- lapply(samples, FUN = function(x) subset_taxa(x, Contamination == FALSE))
samples <- lapply(samples, FUN = function(x) subset_samples(x, Replicate.ID != "WC"))
samples <- lapply(samples, FUN = function(x) subset_samples(x, Replicate.ID != "NTC"))
samples <- lapply(samples, FUN = function(y) prune_samples(samples = sample_sums(y) >= 1000, x = y))

length(rownames(samples$rs_ind@otu_table)) #624
length(rownames(samples$rs_pseudo@otu_table)) #829
length(rownames(samples$rs_pooled@otu_table)) #1132

saveRDS(samples$rs_ind, file = "data/phyloseq_objects/trim_fix/RS21AUG_16S_phyloseq_nt_FALSE_decontam.rds")
saveRDS(samples$rs_pooled, file = "data/phyloseq_objects/trim_fix/RS21AUG_16S_phyloseq_nt_TRUE_decontam.rds")
saveRDS(samples$rs_pseudo, file = "data/phyloseq_objects/trim_fix/RS21AUG_16S_phyloseq_nt_pseudo_decontam.rds")


## Transform zeros
samples_zimp <- lapply(samples, FUN = function(x) cmultRepl(t(otu_table(x)), method = "CZM", output = "p-counts"))

## Transform data with a centred log ratio (CLR)
samples_clr <- lapply(samples_zimp, FUN = function(x) propr(counts = x))

## Perform PCA on the CLR transformed data and plot on covariance biplots
samples_pca <- lapply(samples_clr, FUN = function(x) prcomp(x@logratio))

## Biplots
### Independent
RS_pca_independent <- plot(pca_biplot_unlabelled(pca_data = samples_pca$rs_ind, phyloseq = samples$rs_ind, size = 3, colour = "Atoll", colour_label = "Atoll", shape = "Environment", shape_label = "Environment", title = "Independent"))

### Pooled
RS_pca_pooled <- plot(pca_biplot_unlabelled(pca_data = samples_pca$rs_pooled, phyloseq = samples$rs_pooled, size = 3, colour = "Atoll", colour_label = "Atoll", shape = "Environment", shape_label = "Environment", title = "Pooled"))

### Pseudo
RS_pca_pseudo <- plot(pca_biplot_unlabelled(pca_data = samples_pca$rs_pseudo, phyloseq = samples$rs_pseudo, size = 3, colour = "Atoll", colour_label = "Atoll", shape = "Environment", shape_label = "Environment", title = "Pseudo"))

### Create figure
plot <- plot_grid(RS_pca_independent + theme(legend.position="none"), 
                  # RS_pca_pseudo + theme(legend.position="none"),
                  RS_pca_pooled + theme(legend.position="none"), 
                  nrow = 1)

# add in legend
legend <- get_legend(
  RS_pca_independent +
    guides(color = guide_legend(ncol = 1)) +
    theme(legend.position = "right") +
    theme(legend.text = element_text(size = 12))
)

plot_grid(plot, legend, rel_widths = c(2, .2))


## Divnet analysis
samples <- lapply(samples, FUN = function(x) {sample_data(x)$Site <- paste(sample_data(x)$Atoll, sample_data(x)$Environment, sample_data(x)$Sites, sep = "_"); return(x)})

divnet_alpha <- lapply(samples, FUN = function(x) divnet(W = x, X = "Site", base = names(sort(taxa_sums(x))[round(length(taxa_sums(x))/2)])))
saveRDS(divnet_alpha, file = "data/Divenet_alpha_RS.rds")

RS_alpha_independent_atoll <- data.frame(divnet_alpha$rs_ind$shannon %>% summary, sample_data(samples$rs_ind)) %>%
  distinct(estimate, error, lower, upper, Site, Atoll) %>%
  ggplot(aes(x = Site, y = estimate, col = Atoll)) + geom_point() + geom_segment(aes(x = .data[["Site"]], xend = .data[["Site"]], y = .data[["lower"]], yend = .data[["upper"]])) + xlab("Site") + ylab("Shannon entropy estimate") + theme_bw() + theme(axis.text.x = element_text(angle=90))
saveRDS(RS_alpha_independent_atoll, file = "data/RS_alpha_independent_atoll.rds")

RS_alpha_pooled_atoll <- data.frame(divnet_alpha$rs_pooled$shannon %>% summary, sample_data(samples$rs_pooled)) %>%
  distinct(estimate, error, lower, upper, Site, Atoll) %>%
  ggplot(aes(x = Site, y = estimate, col = Atoll)) + geom_point() + geom_segment(aes(x = .data[["Site"]], xend = .data[["Site"]], y = .data[["lower"]], yend = .data[["upper"]])) + xlab("Site") + ylab("Shannon entropy estimate") + theme_bw() + theme(axis.text.x = element_text(angle=90))
saveRDS(RS_alpha_pooled_atoll, file = "data/RS_alpha_pooled_atoll.rds")

RS_alpha_independent_env <- data.frame(divnet_alpha$rs_ind$shannon %>% summary, sample_data(samples$rs_ind)) %>%
  distinct(estimate, error, lower, upper, Site, Environment) %>%
  ggplot(aes(x = Site, y = estimate, col = Environment)) + geom_point() + geom_segment(aes(x = .data[["Site"]], xend = .data[["Site"]], y = .data[["lower"]], yend = .data[["upper"]])) + xlab("Site") + ylab("Shannon entropy estimate") + theme_bw() + theme(axis.text.x = element_text(angle=90))
saveRDS(RS_alpha_independent_env, file = "data/RS_alpha_independent_env.rds")

RS_alpha_pooled_env <- data.frame(divnet_alpha$rs_pooled$shannon %>% summary, sample_data(samples$rs_pooled)) %>%
  distinct(estimate, error, lower, upper, Site, Environment) %>%
  ggplot(aes(x = Site, y = estimate, col = Environment)) + geom_point() + geom_segment(aes(x = .data[["Site"]], xend = .data[["Site"]], y = .data[["lower"]], yend = .data[["upper"]])) + xlab("Site") + ylab("Shannon entropy estimate") + theme_bw() + theme(axis.text.x = element_text(angle=90))
saveRDS(RS_alpha_pooled_env, file = "data/RS_alpha_pooled_env.rds")

# RS_alpha_pseudo <- data.frame(divnet_alpha$rs_pseudo$shannon %>% summary, sample_data(samples$rs_pseudo)) %>%
#   distinct(estimate, error, lower, upper, Site, Environment) %>%
#   ggplot(aes(x = Site, y = estimate, col = Environment)) + geom_point() + geom_segment(aes(x = .data[["Site"]], xend = .data[["Site"]], y = .data[["lower"]], yend = .data[["upper"]])) + xlab("Site") + ylab("Shannon entropy estimate") + theme_bw() + theme(axis.text.x = element_text(angle=90))

### Create figure
plot <- plot_grid(RS_pca_independent + theme(legend.position="none"), 
                  # RS_pca_pseudo + theme(legend.position="none"),
                  RS_pca_pooled + theme(legend.position="none"), 
                  nrow = 1)

# add in legend
legend <- get_legend(
  RS_pca_independent +
    guides(color = guide_legend(nrow = 1)) +
    theme(legend.position = "bottom") +
    theme(legend.text = element_text(size = 12))
)

plot_grid(plot, legend, rel_heights = c(2, .2), ncol = 1)

plot <- plot_grid(RS_alpha_independent_atoll + theme(legend.position="none") + theme(axis.text.x = element_blank()), 
                  RS_alpha_pooled_atoll + theme(legend.position="none")  + theme(axis.text.x = element_blank()) + ylab(NULL), 
                  RS_alpha_independent_env + theme(legend.position="none") + ggtitle(NULL), 
                  RS_alpha_pooled_env + theme(legend.position="none") + ylab(NULL) + ggtitle(NULL))
                  # nrow = 1)
plot

# add in legend
legend <- get_legend(
  RS_alpha_independent +
    guides(color = guide_legend(nrow = 1)) +
    theme(legend.position = "bottom")
)

plot_grid(plot, legend, ncol = 1, rel_heights = c(1, .1))

### Create figure
plot <- plot_grid(RS_alpha_independent_atoll + theme(legend.position="none") + theme(axis.text.x = element_blank()), 
                  RS_alpha_pooled_atoll + theme(legend.position="none") + ylab(NULL)  + theme(axis.text.x = element_blank()), 
                  nrow = 1)
plot
# add in legend
legend <- get_legend(
  RS_alpha_independent_atoll +
    guides(color = guide_legend(nrow = 1)) +
    theme(legend.position = "bottom") +
    theme(legend.text = element_text(size = 12))
)

RS_alpha_atoll <- plot_grid(plot, legend, rel_heights = c(1, .1), ncol = 1)

plot <- plot_grid(RS_alpha_independent_env + theme(legend.position="none"), 
                  RS_alpha_pooled_env + theme(legend.position="none") + ylab(NULL), 
                  nrow = 1)
plot
# add in legend
legend <- get_legend(
  RS_alpha_independent_env +
    guides(color = guide_legend(nrow = 1)) +
    theme(legend.position = "bottom") +
    theme(legend.text = element_text(size = 12))
)

RS_alpha_env <- plot_grid(plot, legend, rel_heights = c(1, .1), ncol = 1)

alphas <- plot_grid(RS_alpha_atoll, RS_alpha_env, ncol = 1, rel_heights = c(1.2, 2))
alphas
