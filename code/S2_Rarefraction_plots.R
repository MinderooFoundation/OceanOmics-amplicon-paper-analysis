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
CoCos16S_false  <- readRDS('data/phyloseq_objects/CoCosV10I_16S_phyloseq_nt_FALSE_decontaminated.rds')
CoCos16S_true  <- readRDS('data/phyloseq_objects/CoCosV10I_16S_phyloseq_nt_TRUE_decontaminated.rds')
CoCos16S_pseudo <- readRDS('data/phyloseq_objects/CoCosV10I_16S_phyloseq_nt_pseudo_decontaminated.rds')

NW_false  <- readRDS('data/phyloseq_objects/Pool_FALSE_KWEST_16S_phyloseq_nt_decontaminated.rds')
NW_true  <- readRDS('data/phyloseq_objects/Pool_TRUE_KWEST_16S_phyloseq_nt_decontaminated.rds')
NW_pseudo <- readRDS('data/phyloseq_objects/CoCosV10I_16S_phyloseq_nt_pseudo_decontaminated.rds')

RS16S_false  <- readRDS('data/phyloseq_objects/RS21AUG_16S_phyloseq_nt_FALSE_decontam.rds')
RS16S_true  <- readRDS('data/phyloseq_objects/RS21AUG_16S_phyloseq_nt_TRUE_decontam.rds')
RS16S_pseudo <- readRDS('data/phyloseq_objects/RS21AUG_16S_phyloseq_nt_pseudo_decontam.rds')

# # create plots for CococsI
Cocos_false_curve <- ggrare(CoCos16S_false, step = 10, label = NULL,
                            plot = TRUE, parallel = FALSE, se = TRUE)

saveRDS(Cocos_false_curve, file = "data/rarefaction_objects/Cocos_false_curve_plot.rds")
# Cocos_false_curve <- readRDS(file = "data/rarefaction_objects/Cocos_false_curve_plot.rds")

Cocos_true_curve <- ggrare(CoCos16S_true, step = 10, label = NULL,
                           plot = TRUE, parallel = FALSE, se = TRUE)

saveRDS(Cocos_true_curve, file = "data/rarefaction_objects/Cocos_true_curve_plot.rds")
# Cocos_true_curve <- readRDS(file = "data/rarefaction_objects/Cocos_true_curve_plot.rds")



# Plots for Rowley Shoals dataset
RS_false_curve <- ggrare(RS16S_false, step = 10, label = NULL,
                           plot = TRUE, parallel = FALSE, se = TRUE)

saveRDS(RS_false_curve, file = "data/rarefaction_objects/RS_false_curve_plot.rds")
# RS_false_curve <- readRDS(file = "data/rarefaction_objects/RS_false_curve_plot.rds")

RS_true_curve <- ggrare(RS16S_true, step = 10, label = NULL,
                        plot = TRUE, parallel = FALSE, se = TRUE)

saveRDS(RS_true_curve, file = "RS_true_curve_plot.rds")
# RS_true_curve <- readRDS(file = "data/rarefaction_objects/RS_true_curve_plot.rds")



# Plot curves for Northwest data set
NW_false_curve <- ggrare(NW_false, step = 10, label = NULL,
                         plot = TRUE, parallel = FALSE, se = TRUE)

saveRDS(NW_false_curve, file = "data/phyloseq_objects/NW_false_curve_plot.rds")
# NW_false_curve <- readRDS(file = "data/rarefaction_objects/NW_false_curve_plot.rds")

NW_true_curve <- ggrare(NW_true, step = 10, label = NULL,
                        plot = TRUE, parallel = FALSE, se = TRUE)

saveRDS(NW_true_curve, file = "data/phyloseq_objects/NW_true_curve_plot.rds")
# NW_true_curve <- readRDS(file = "data/rarefaction_objects/NW_true_curve_plot.rds")



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
