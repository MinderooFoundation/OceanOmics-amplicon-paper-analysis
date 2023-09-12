#-----------------------------------------------------------------------------------------------------------------
# Generate beta diversity PCA plots 

# Load libraries
LoadandCite(pkgs=c("repmis", "knitr", "tinytex", "phyloseq", "ggplot2", 
                   "zCompositions", "propr", "compositions", "ggfortify", 
                   "ALDEx2", "EnhancedVolcano", "plyr", "microViz", "decontam", 
                   "patchwork", "readxl", "DivNet", "dplyr", "cowplot"))

## Load data
metadata           <- data.frame(read_excel(path = "data/metadata/V10_CKI_eDNA_metadata_P1.xlsx"))
samples            <- list()
samples$cki_ind    <- readRDS(file = "data/phyloseq_objects/CoCosV10I_16S_phyloseq_nt_FALSE.rds")
samples$cki_pooled <- readRDS(file = "data/phyloseq_objects/CoCosV10I_16S_phyloseq_nt_TRUE.rds")
samples$cki_pseudo <- readRDS(file = "data/phyloseq_objects/CoCosV10I_16S_phyloseq_nt_pseudo.rds")

length(rownames(samples$cki_ind@otu_table))
length(rownames(samples$cki_pseudo@otu_table))
length(rownames(samples$cki_pooled@otu_table))

## Create bioplot function
pca_biplot <- function(pca_data, phyloseq, colour, colour_label, shape, shape_label, title) {
  long_arrows <- names(which(abs(pca_data$rotation[ ,"PC1"]) > 0.1 | abs(pca_data$rotation[ ,"PC2"]) > 0.1))
  loadings_label <- tax_table(phyloseq)@.Data[ , "LCA_ASV"]
  loadings_label[!names(loadings_label) %in% long_arrows] <- ""
  
  biplot <- autoplot(pca_data,
                     data = sample_data(phyloseq),
                     colour = colour, shape = shape,
                     label = TRUE, label.label = sample_data(phyloseq)$Sample,
                     label.size = 3, label.vjust = -0.5,
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
                     loadings.label.size = 5, loadings.label.colour = "gray30", loadings.label.repel = TRUE) + theme_bw() + guides(colour = guide_legend(title = colour_label), shape = guide_legend(title = shape_label)) + theme(legend.title = element_text(size = 12), legend.text = element_text(size = 11)) + ggtitle(label = title)
  
  return(biplot)
}


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

# saveRDS(samples$cki_ind, file = "data/phyloseq_objects/CoCosV10I_16S_phyloseq_nt_FALSE_decontaminated.rds")
# saveRDS(samples$cki_pooled, file = "data/phyloseq_objects/CoCosV10I_16S_phyloseq_nt_TRUE_decontaminated.rds")
# saveRDS(samples$cki_pseudo, file = "data/phyloseq_objects/CoCosV10I_16S_phyloseq_nt_pseudo_decontaminated.rds")

## Transform zeros
samples_zimp <- lapply(samples, FUN = function(x) cmultRepl(t(otu_table(x)), method = "CZM", output = "p-counts"))

## Transform data with a centred log ratio (CLR)
samples_clr <- lapply(samples_zimp, FUN = function(x) propr(counts = x))

## Perform PCA on the CLR transformed data and plot on covariance biplots
samples_pca <- lapply(samples_clr, FUN = function(x) prcomp(x@logratio))

## Biplots
### Independent
Cocos_pca_independent <- plot(pca_biplot_unlabelled(pca_data = samples_pca$cki_ind, phyloseq = samples$cki_ind, size = 3, colour = "site_id", colour_label = "Site", shape = "EEZ", shape_label = "EEZ", title = "Independent"))

### Pooled
Cocos_pca_pooled <- plot(pca_biplot_unlabelled(pca_data = samples_pca$cki_pooled, phyloseq = samples$cki_pooled, size = 3, colour = "site_id", colour_label = "Site", shape = "EEZ", shape_label = "EEZ", title = "Pooled"))

### Pseudo
Cocos_pca_pseudo <- plot(pca_biplot_unlabelled(pca_data = samples_pca$cki_pseudo, phyloseq = samples$cki_pseudo, size = 3, colour = "site_id", colour_label = "Site", shape = "EEZ", shape_label = "EEZ", title = "Pseudo"))
Cocos_pca_pseudo + scale_color_viridis_d(direction = -1, option = "H")
Cocos_pca_pseudo + scale_color_viridis_d(direction = -1, option = "D")
Cocos_pca_pseudo + scale_color_d3(palette = "category20")
Cocos_pca_pseudo + scale_color_tableau(palette = "Classic 20")

### Create figure
plot_pca <- plot_grid(Cocos_pca_independent + theme(legend.position="none")  + scale_color_viridis_d(direction = -1, option = "D"), 
          # Cocos_pca_pseudo + theme(legend.position="none"),
          Cocos_pca_pooled + theme(legend.position="none") + scale_color_viridis_d(direction = -1, option = "D"),
          nrow = 1)
plot_pca
# add in legend
legend <- get_legend(
  Cocos_pca_independent +
    guides(color = guide_legend(ncol = 1)) +
    theme(legend.position = "right") +
    theme(legend.text = element_text(size = 12)) +
    scale_color_viridis_d(direction = -1, option = "D")
)

plot_grid(plot_pca, legend, rel_widths = c(2, .2))








## Divnet analysis
# Remove the sample that are not in all results
samples_to_keep <- Reduce(intersect, list(rownames(sample_data(samples$cki_ind)), rownames(sample_data(samples$cki_pseudo)), rownames(sample_data(samples$cki_pooled))))
samples$cki_ind <- prune_samples(samples_to_keep, samples$cki_ind)
samples$cki_pseudo <- prune_samples(samples_to_keep, samples$cki_pseudo)
samples$cki_pooled <- prune_samples(samples_to_keep, samples$cki_pooled)

divnet_alpha <- lapply(samples, FUN = function(x) divnet(W = x, X = "site_id", base = names(sort(taxa_sums(x))[round(length(taxa_sums(x))/2)])))
divnet_alpha_cki_ind <- divnet_alpha$cki_ind$shannon
divnet_alpha_cki_pseudo <- divnet_alpha$cki_pseudo$shannon
divnet_alpha_cki_pooled <- divnet_alpha$cki_pooled$shannon

saveRDS(divnet_alpha_cki_ind, file = "data/divnet_alpha_cki_ind.rds")
saveRDS(divnet_alpha_cki_pseudo, file = "data/divnet_alpha_cki_pseudo.rds")
saveRDS(divnet_alpha_cki_pooled, file = "data/divnet_alpha_cki_pooled.rds")

divnet_alpha_cki_ind <- readRDS(file = "data/divnet_alpha_cki_ind.rds")
divnet_alpha_cki_pseudo <- readRDS(file = "data/divnet_alpha_cki_pseudo.rds")
divnet_alpha_cki_pooled <- readRDS(file = "data/divnet_alpha_cki_pooled.rds")

Cocos_alpha_independent <- data.frame(divnet_alpha_cki_ind %>% summary, sample_data(samples$cki_ind)) %>%
  distinct(estimate, error, lower, upper, site_id, EEZ) %>%
  ggplot(aes(x = site_id, y = estimate, col = site_id)) + geom_point(aes(shape = EEZ), size=2) + scale_color_viridis_d(direction = -1, option = "D") + geom_segment(aes(x = .data[["site_id"]], xend = .data[["site_id"]], y = .data[["lower"]], yend = .data[["upper"]])) + guides(colour = FALSE) + xlab("Site identifier") + ylab("Shannon entropy estimate") + theme_bw()

Cocos_alpha_pooled <- data.frame(divnet_alpha_cki_pooled %>% summary, sample_data(samples$cki_pooled)) %>%
  distinct(estimate, error, lower, upper, site_id, EEZ) %>%
  ggplot(aes(x = site_id, y = estimate, col = site_id)) + geom_point(aes(shape = EEZ), size=2) + scale_color_viridis_d(direction = -1, option = "D") + geom_segment(aes(x = .data[["site_id"]], xend = .data[["site_id"]], y = .data[["lower"]], yend = .data[["upper"]])) + guides(colour = FALSE) + xlab("Site identifier") + ylab("Shannon entropy estimate") + theme_bw()

Cocos_alpha_pseudo <- data.frame(divnet_alpha$cki_pseudo$shannon %>% summary, sample_data(samples$cki_pseudo)) %>%
  distinct(estimate, error, lower, upper, site_id) %>%
  ggplot(aes(x = site_id, y = estimate, col = site_id)) + geom_point() + geom_segment(aes(x = .data[["site_id"]], xend = .data[["site_id"]], y = .data[["lower"]], yend = .data[["upper"]])) + guides(colour = FALSE) + xlab("Site identifier") + ylab("Shannon entropy estimate") + theme_bw()


### Create figure
plot_cocos <- plot_grid(Cocos_pca_independent + theme(legend.position="none") + scale_color_viridis_d(direction = -1, option = "D"), 
                  Cocos_pca_pooled + theme(legend.position="none") + scale_color_viridis_d(direction = -1, option = "D"),
                  Cocos_alpha_independent + theme(legend.position="none") + ggtitle(NULL) + coord_cartesian(ylim = c(0, 5)) + scale_color_viridis_d(direction = -1, option = "D"), 
                  Cocos_alpha_pooled + theme(legend.position="none") + ylab(NULL) + ggtitle(NULL) + coord_cartesian(ylim = c(0, 5)) + scale_color_viridis_d(direction = -1, option = "D"),
                  rel_heights = c(2, 1), labels = c("A","B","C","D"), label_size = 14)
                  # nrow = 1)
plot_cocos

# add in legend
legend <- get_legend(
  Cocos_alpha_independent +
    guides(color = guide_legend(ncol = 1)) +
    theme(legend.position = "right") +
    theme(legend.text = element_text(size = 12)) +
    theme(legend.title = element_text("Site ID"))
)

plot_grid(plot_cocos, legend, rel_widths = c(2, .2), align = "hv", axis = "tblr")
