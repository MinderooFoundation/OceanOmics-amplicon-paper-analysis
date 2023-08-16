#=====================================================================
# Plot ASV distribution by habitat type to see if the pattern remains
# Seb Rauschert
#=====================================================================

library(tidyverse)
library(phyloseq)

RS_false  <- readRDS('data/phyloseq_objects_decontaminated/RSV5_16S_individual_phyloseq_nt_decontaminated.rds')
RS_true  <- readRDS('data/phyloseq_objects_decontaminated/RSV5_16S_pooled_phyloseq_nt_decontaminated.rds')
RS_pseudo <- readRDS('data/phyloseq_objects_decontaminated/RSV5_16S_pseudo_phyloseq_nt_decontaminated.rds')

#=====================================================================
# Subset data for habitate comparison plots: pool = TRUE
# Clerke Habitats

RS_Clerke_true        <- subset_samples(RS_true, Atoll == "Clerke")
RS_Clerke_Lagoon_true <- subset_samples(RS_Clerke_true, Environment == "Lagoon")
RS_Clerke_Slope_true  <- subset_samples(RS_Clerke_true, Environment == "Slope")

# Mermaid Habitats
RS_Mermaid_true        <- subset_samples(RS_true, Atoll == "Mermaid")
RS_Mermaid_Lagoon_true <- subset_samples(RS_Mermaid_true, Environment == "Lagoon")
RS_Mermaid_Slope_true  <- subset_samples(RS_Mermaid_true, Environment == "Slope")

# Imperious Habitats
RS_Imperious_true        <- subset_samples(RS_true, Atoll == "Imperious")
RS_Imperious_Lagoon_true <- subset_samples(RS_Imperious_true, Environment == "Lagoon")
RS_Imperious_Slope_true  <- subset_samples(RS_Imperious_true, Environment == "Slope")
#=====================================================================



#=====================================================================
# Subset data for habitate comparison plots: pool = FALSE
# Clerke Habitats

RS_Clerke_false        <- subset_samples(RS_false, Atoll == "Clerke")
RS_Clerke_Lagoon_false <- subset_samples(RS_Clerke_false, Environment == "Lagoon")
RS_Clerke_Slope_false  <- subset_samples(RS_Clerke_false, Environment == "Slope")

# Mermaid Habitats
RS_Mermaid_false        <- subset_samples(RS_false, Atoll == "Mermaid")
RS_Mermaid_Lagoon_false <- subset_samples(RS_Mermaid_false, Environment == "Lagoon")
RS_Mermaid_Slope_false  <- subset_samples(RS_Mermaid_false, Environment == "Slope")

# Imperious Habitats
RS_Imperious_false        <- subset_samples(RS_false, Atoll == "Imperious")
RS_Imperious_Lagoon_false <- subset_samples(RS_Imperious_false, Environment == "Lagoon")
RS_Imperious_Slope_false  <- subset_samples(RS_Imperious_false, Environment == "Slope")
#=====================================================================



#=====================================================================
# Subset data for habitate comparison plots: pool = "pseudo"
# Clerke Habitats

RS_Clerke_pseudo        <- subset_samples(RS_pseudo, Atoll == "Clerke")
RS_Clerke_Lagoon_pseudo <- subset_samples(RS_Clerke_pseudo, Environment == "Lagoon")
RS_Clerke_Slope_pseudo  <- subset_samples(RS_Clerke_pseudo, Environment == "Slope")

# Mermaid Habitats
RS_Mermaid_pseudo        <- subset_samples(RS_pseudo, Atoll == "Mermaid")
RS_Mermaid_Lagoon_pseudo <- subset_samples(RS_Mermaid_pseudo, Environment == "Lagoon")
RS_Mermaid_Slope_pseudo  <- subset_samples(RS_Mermaid_pseudo, Environment == "Slope")

# Imperious Habitats
RS_Imperious_pseudo        <- subset_samples(RS_pseudo, Atoll == "Imperious")
RS_Imperious_Lagoon_pseudo <- subset_samples(RS_Imperious_pseudo, Environment == "Lagoon")
RS_Imperious_Slope_pseudo  <- subset_samples(RS_Imperious_pseudo, Environment == "Slope")
#=====================================================================


# Plot the same thing, by habitat type
# FUNCTION for plotting
asv_plot<- function(pooled_data, 
         individual_data,
         pseudo_data, 
         title){

    RS_false_seqtab <- t(otu_table(individual_data))
    RS_true_seqtab  <- t(otu_table(pooled_data))
    RS_pseudo_seqtab <- t(otu_table(pseudo_data))
    
    
    RS_false_seqtab <- RS_false_seqtab[!(row.names(RS_false_seqtab) %in% "RS1_ME_S4_1(2)"),]
    RS_true_seqtab <- RS_true_seqtab[!(row.names(RS_true_seqtab) %in% "RS1_ME_S4_1(2)"),]
    RS_pseudo_seqtab <- RS_pseudo_seqtab[!(row.names(RS_pseudo_seqtab) %in% "RS1_ME_S4_1(2)"),]
    
    
    nsam <- dim(RS_true_seqtab)[1]
    df.obs <- data.frame(observed=c(rowSums(RS_false_seqtab>0), rowSums(RS_pseudo_seqtab>0), rowSums(RS_true_seqtab>0)),
                         Mode=rep(c("independent", "pseudo", "pooled"), each=nsam),
                         rank=rank(rowSums(RS_true_seqtab>0))) #%>%
    
    mode_comp <- ggplot(data = df.obs, aes(x = rank, y = observed, color=Mode)) + 
      geom_point(alpha=0.5, size = 3) +
      xlab("Samples") + 
      ylab("Observed ASVs") +
      ggtitle(title) +
      theme_bw(base_size = 18) +
      geom_smooth(method = "lm")
    
    mode_comp

}

# FUNCTION

# Clerke
CLERKE <- asv_plot(RS_Clerke_true, RS_Clerke_false, RS_Clerke_pseudo, "Clerke")

# Imperious
IMPERIOUS <- asv_plot(RS_Imperious_true, RS_Imperious_false, RS_Imperious_pseudo, "Imperious")

# Mermaid
MERMAID <- asv_plot(RS_Mermaid_true, RS_Mermaid_false, RS_Mermaid_pseudo, "Mermaid")

# ggpubr::ggarrange(CLERKE, IMPERIOUS, MERMAID, ncol = 3)


# Clerke Lagoon
CLERKE_lagoon <- asv_plot(RS_Clerke_Lagoon_true, RS_Clerke_Lagoon_false, RS_Clerke_Lagoon_pseudo, "Clerke Lagoon")

# Imperious Lagoon
IMPERIOUS_lagoon <- asv_plot(RS_Imperious_Lagoon_true, RS_Imperious_Lagoon_false, RS_Imperious_Lagoon_pseudo, "Imperious Lagoon")

# Mermaid Lagoon
MERMAID_lagoon <- asv_plot(RS_Mermaid_Lagoon_true, RS_Mermaid_Lagoon_false, RS_Mermaid_Lagoon_pseudo, "Mermaid Lagoon")


# Clerke Slope
CLERKE_slope <- asv_plot(RS_Clerke_Slope_true, RS_Clerke_Slope_false, RS_Clerke_Slope_pseudo, "Clerke Slope")

# Imperious Slope
IMPERIOUS_slope <- asv_plot(RS_Imperious_Slope_true, RS_Imperious_Slope_false, RS_Imperious_Slope_pseudo, "Imperious Slope")

# Mermaid Slope
MERMAID_slope <- asv_plot(RS_Mermaid_Slope_true, RS_Mermaid_Slope_false, RS_Mermaid_Slope_pseudo, "Mermaid Slope")

