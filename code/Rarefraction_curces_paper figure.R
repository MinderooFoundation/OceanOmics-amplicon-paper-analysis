# Load libraries
library(ggplot2)
library(vegan)
library(tidyverse)
library(phyloseq)
library(cowplot)
library(RColorBrewer)

# Read in ggrare function
ggrare <- function(physeq_object, step = 10, label = NULL, color = NULL, plot = TRUE, parallel = FALSE, se = TRUE) {
  
  x <- methods::as(phyloseq::otu_table(physeq_object), "matrix")
  if (phyloseq::taxa_are_rows(physeq_object)) { x <- t(x) }
  
  ## This script is adapted from vegan `rarecurve` function
  tot <- rowSums(x)
  S <- rowSums(x > 0)
  nr <- nrow(x)
  
  rarefun <- function(i) {
    cat(paste("rarefying sample", rownames(x)[i]), sep = "\n")
    n <- seq(1, tot[i], by = step)
    if (n[length(n)] != tot[i]) {
      n <- c(n, tot[i])
    }
    y <- vegan::rarefy(x[i, ,drop = FALSE], n, se = se)
    if (nrow(y) != 1) {
      rownames(y) <- c(".S", ".se")
      return(data.frame(t(y), Size = n, Sample = rownames(x)[i]))
    } else {
      return(data.frame(.S = y[1, ], Size = n, Sample = rownames(x)[i]))
    }
  }
  if (parallel) {
    out <- parallel::mclapply(seq_len(nr), rarefun, mc.preschedule = FALSE)
  } else {
    out <- lapply(seq_len(nr), rarefun)
  }
  df <- do.call(rbind, out)
  
  # Get sample data
  if (!is.null(phyloseq::sample_data(physeq_object, FALSE))) {
    sdf <- methods::as(phyloseq::sample_data(physeq_object), "data.frame")
    sdf$Sample <- rownames(sdf)
    data <- merge(df, sdf, by = "Sample")
    labels <- data.frame(x = tot, y = S, Sample = rownames(x))
    labels <- merge(labels, sdf, by = "Sample")
  }
  
  # Add, any custom-supplied plot-mapped variables
  if ( length(color) > 1 ) {
    data$color <- color
    names(data)[names(data) == "color"] <- deparse(substitute(color))
    color <- deparse(substitute(color))
  }
  
  if ( length(label) > 1 ) {
    labels$label <- label
    names(labels)[names(labels) == "label"] <- deparse(substitute(label))
    label <- deparse(substitute(label))
  }
  
  p <- ggplot2::ggplot(data = data,
                       ggplot2::aes_string(x = "Size",
                                           y = ".S",
                                           group = "Sample",
                                           color = color))
  
  p <- p + ggplot2::labs(x = "Sequence Sample Size", y = "Species Richness")
  
  if (!is.null(label)) {
    p <- p + ggplot2::geom_text(data = labels,
                                ggplot2::aes_string(x = "x",
                                                    y = "y",
                                                    label = label,
                                                    color = color),
                                size = 4, hjust = 0)
  }
  
  p <- p + ggplot2::geom_line()
  if (se) { ## add standard error if available
    p <- p +
      ggplot2::geom_ribbon(ggplot2::aes_string(ymin = ".S - .se",
                                               ymax = ".S + .se",
                                               color = NULL,
                                               fill = color),
                           alpha = 0.2)
  }
  if (plot) {
    plot(p)
  }
  invisible(p)
}

# Read in phyloseq objects
CoCos16S_false  <- readRDS('data/phyloseq_objects/trim_fix/CoCosV10I_16S_phyloseq_nt_FALSE_decontaminated.rds')
CoCos16S_true  <- readRDS('data/phyloseq_objects/trim_fix/CoCosV10I_16S_phyloseq_nt_TRUE_decontaminated.rds')
CoCos16S_pseudo <- readRDS('data/phyloseq_objects/trim_fix/CoCosV10I_16S_phyloseq_nt_pseudo_decontaminated.rds')

NW_false  <- readRDS('data/phyloseq_objects/Pool_FALSE_KWEST_16S_phyloseq_nt_decontaminated.rds')
NW_true  <- readRDS('data/phyloseq_objects/Pool_TRUE_KWEST_16S_phyloseq_nt_decontaminated.rds')
NW_pseudo <- readRDS('data/phyloseq_objects/CoCosV10I_16S_phyloseq_nt_pseudo_decontaminated.rds')

RS16S_false  <- readRDS('data/phyloseq_objects/trim_fix/RS21AUG_16S_phyloseq_nt_FALSE_decontam.rds')
RS16S_true  <- readRDS('data/phyloseq_objects/trim_fix/RS21AUG_16S_phyloseq_nt_TRUE_decontam.rds')
RS16S_pseudo <- readRDS('data/phyloseq_objects/trim_fix/RS21AUG_16S_phyloseq_nt_pseudo_decontam.rds')

# # create plots for CococsI
# Cocos_false_curve <- ggrare(CoCos16S_false, step = 10, label = NULL, mycolor='#1B9E77',
#                             plot = TRUE, parallel = FALSE, se = TRUE)

Cocos_false_curve <- ggrare(CoCos16S_false, step = 10, label = NULL,
                            plot = TRUE, parallel = FALSE, se = TRUE)

saveRDS(Cocos_false_curve, file = "data/rarefaction_objects/Cocos_false_curve_plot.rds")
Cocos_false_curve <- readRDS(file = "data/rarefaction_objects/Cocos_false_curve_plot.rds")

# Cocos_true_curve <- ggrare(CoCos16S_true, step = 10, label = NULL,  mycolor='#D95F02',
#                             plot = TRUE, parallel = FALSE, se = TRUE)

Cocos_true_curve <- ggrare(CoCos16S_true, step = 10, label = NULL,
                           plot = TRUE, parallel = FALSE, se = TRUE)

saveRDS(Cocos_true_curve, file = "data/rarefaction_objects/Cocos_true_curve_plot.rds")
Cocos_true_curve <- readRDS(file = "data/rarefaction_objects/Cocos_true_curve_plot.rds")

# Testing other metadata options analysis Cocos data set
Cocos_false_curve_site <- ggrare(CoCos16S_false, step = 10, label = NULL, color = "site_id",
                                 plot = TRUE, parallel = FALSE, se = TRUE)

saveRDS(Cocos_false_curve_site, file = "data/rarefaction_objects/Cocos_false_curve_plot_site.rds")

Cocos_false_curve_method <- ggrare(CoCos16S_false, step = 10, label = NULL, color = "sampling_method...20",
                                   plot = TRUE, parallel = FALSE, se = TRUE)

saveRDS(Cocos_false_curve_method, file = "data/rarefaction_objects/Cocos_false_curve_plot_method.rds")

Cocos_true_curve_site <- ggrare(CoCos16S_true, step = 10, label = NULL, color = "site_id",
                           plot = TRUE, parallel = FALSE, se = TRUE)

saveRDS(Cocos_true_curve_site, file = "data/rarefaction_objects/Cocos_true_curve_plot_site.rds")

Cocos_true_curve_method <- ggrare(CoCos16S_true, step = 10, label = NULL, color = "sampling_method...20",
                                          plot = TRUE, parallel = FALSE, se = TRUE)

saveRDS(Cocos_true_curve_method, file = "data/rarefaction_objects/Cocos_true_curve_plot_method.rds")



# Plots for Rowley Shoals dataset
# RS_false_curve <- ggrare(RS16S_false, step = 10, label = NULL,  mycolor='#1B9E77',
#                          plot = TRUE, parallel = FALSE, se = TRUE)
RS_false_curve <- ggrare(RS16S_false, step = 10, label = NULL,
                           plot = TRUE, parallel = FALSE, se = TRUE)

saveRDS(RS_false_curve, file = "data/rarefaction_objects/RS_false_curve_plot.rds")
RS_false_curve <- readRDS(file = "data/rarefaction_objects/RS_false_curve_plot.rds")
# test <- readRDS("RS_false_curve_plot.rds")

# RS_true_curve <- ggrare(RS16S_true, step = 10, label = NULL, mycolor='#D95F02',
#                         plot = TRUE, parallel = FALSE, se = TRUE)
RS_true_curve <- ggrare(RS16S_true, step = 10, label = NULL,
                        plot = TRUE, parallel = FALSE, se = TRUE)

saveRDS(RS_true_curve, file = "RS_true_curve_plot.rds")
RS_true_curve <- readRDS(file = "data/rarefaction_objects/RS_true_curve_plot.rds")

# Testing metadata variables
RS_false_curve_atoll <- ggrare(RS16S_false, step = 10, label = NULL, color = "Atoll",
                         plot = TRUE, parallel = FALSE, se = TRUE)

saveRDS(RS_false_curve_atoll, file = "RS_false_curve_plot_atoll.rds")

RS_true_curve_atoll <- ggrare(RS16S_true, step = 10, label = NULL, color = "Atoll",
                        plot = TRUE, parallel = FALSE, se = TRUE)

saveRDS(RS_true_curve_atoll, file = "RS_true_curve_plot_atoll.rds")




# Plot curves for Northwest data set
# NW_false_curve <- ggrare(NW_false, step = 10, label = NULL,  mycolor='#1B9E77',
#                         plot = TRUE, parallel = FALSE, se = TRUE)
NW_false_curve <- ggrare(NW_false, step = 10, label = NULL,
                         plot = TRUE, parallel = FALSE, se = TRUE)

saveRDS(NW_false_curve, file = "data/phyloseq_objects/NW_false_curve_plot.rds")
NW_false_curve <- readRDS(file = "data/rarefaction_objects/NW_false_curve_plot.rds")

# NW_true_curve <- ggrare(NW_true, step = 10, label = NULL, mycolor='#D95F02',
#                         plot = TRUE, parallel = FALSE, se = TRUE)
NW_true_curve <- ggrare(NW_true, step = 10, label = NULL,
                        plot = TRUE, parallel = FALSE, se = TRUE)

saveRDS(NW_true_curve, file = "data/phyloseq_objects/NW_true_curve_plot.rds")
NW_true_curve <- readRDS(file = "data/rarefaction_objects/NW_true_curve_plot.rds")

## Investigate potential ecological differences
NW_false_curve_bioregion <- ggrare(NW_false, step = 10, label = NULL, color = "Bioregion",
                                  plot = TRUE, parallel = FALSE, se = TRUE)
saveRDS(NW_false_curve_bioregion, file = "data/phyloseq_objects/NW_false_curve_plot_bioregion.rds")

NW_true_curve_bioregion <- ggrare(NW_true, step = 10, label = NULL, color = "Bioregion",
                        plot = TRUE, parallel = FALSE, se = TRUE)
saveRDS(NW_true_curve_bioregion, file = "data/phyloseq_objects/NW_true_curve_plot_bioregion.rds")

NW_false_curve_subregion <- ggrare(NW_false, step = 10, label = NULL, color = "Subregion",
                                   plot = TRUE, parallel = FALSE, se = TRUE)
saveRDS(NW_false_curve_subregion, file = "NW_false_curve_plot_subregion.rds")

NW_true_curve_subregion <- ggrare(NW_true, step = 10, label = NULL, color = "Subregion",
                                  plot = TRUE, parallel = FALSE, se = TRUE)
saveRDS(NW_true_curve_subregion, file = "data/phyloseq_objects/NW_true_curve_plot_subregion.rds")

NW_false_curve_habitat <- ggrare(NW_false, step = 10, label = NULL, color = "Habitat",
                                   plot = TRUE, parallel = FALSE, se = TRUE)
saveRDS(NW_true_curve, file = "NW_true_curve_plot.rds")

NW_true_curve_habitat <- ggrare(NW_true, step = 10, label = NULL, color = "Habitat",
                                  plot = TRUE, parallel = FALSE, se = TRUE)
saveRDS(NW_true_curve, file = "NW_true_curve_plot.rds")


# Create plot grid with same axis
# plot <- plot_grid(Cocos_false_curve + theme(legend.position="none", line = ) + xlab(NULL)  + coord_cartesian(ylim = c(0, 500)),
#                   RS_false_curve + theme(legend.position="none") + xlab(NULL) + ylab(NULL) + coord_cartesian(ylim = c(0, 500)),
#                   NW_false_curve + theme(legend.position="none") + xlab(NULL) + ylab(NULL)  + coord_cartesian(ylim = c(0, 500)),
#                   Cocos_true_curve + theme(legend.position="none")  + coord_cartesian(ylim = c(0, 500)),
#                   RS_true_curve + theme(legend.position="none") + ylab(NULL) + coord_cartesian(ylim = c(0, 500)),
#                   NW_true_curve + theme(legend.position="none")+ ylab(NULL)  + coord_cartesian(ylim = c(0, 500)),
#                   labels = 'AUTO', label_size = 12)
# plot

# Create plot
plot <- plot_grid(Cocos_false_curve + theme(legend.position="none") + xlab(NULL)  + coord_cartesian(ylim = c(0, 80)) + ggtitle("Cocos Islands") + theme(plot.title = element_text(hjust = 0.5)),
                  RS_false_curve + theme(legend.position="none") + xlab(NULL) + ylab(NULL) + coord_cartesian(ylim = c(0, 500)) + ggtitle("Rowley Shoals") + theme(plot.title = element_text(hjust = 0.5)),
                  NW_false_curve + theme(legend.position="none") + xlab(NULL) + ylab(NULL)  + coord_cartesian(ylim = c(0, 80)) + ggtitle("NW WA") + theme(plot.title = element_text(hjust = 0.5)),
                  Cocos_true_curve + theme(legend.position="none"),
                  RS_true_curve + theme(legend.position="none") + ylab(NULL) + coord_cartesian(ylim = c(0, 500)),
                  NW_true_curve + theme(legend.position="none"),
                  labels = 'AUTO',
                  label_size = 12)
plot

# Cocos ecolgy plot
plot_cocos <- plot_grid(Cocos_false_curve_site + theme(legend.position="none"),
                  Cocos_true_curve_site,
                  Cocos_false_curve_method + theme(legend.position="none"),
                  Cocos_true_curve_method,
                  labels = 'AUTO')
plot_cocos


# NW WA ecolgy plot
plot_NW <- plot_grid(NW_false_curve_bioregion + theme(legend.position="none"),
                     NW_true_curve_bioregion,
                     NW_false_curve_subregion + theme(legend.position="none"),
                     NW_true_curve_subregion,
                     NW_false_curve_habitat + theme(legend.position="none"),
                     NW_true_curve_habitat,
                     labels = 'AUTO', nrow = 3)
plot_NW

## Colour scheme from other figure
# brewer.pal(n = 3, "Dark2")
# display.brewer.pal(n = 3, "Dark2")
## FASLE = #1B9E77
## TRUE = #D95F02
