# Install required R packages for the analysis

install.packages("tidyverse")
install.packages("tidyr")
install.packages("stringr")
install.packages("mapdata")
install.packages("sf")
install.packages("sp")
install.packages("rgdal")
install.packages("terra")
install.packages("viridis")
install.packages("remotes")   ## run this line if you do not already have remotes installed
remotes::install_github("adw96/breakaway")
remotes::install_github("adw96/DivNet")
install.packages("BiocManager")
BiocManager::install(c("phyloseq", "microbiome", "ComplexHeatmap"), update = FALSE)
install.packages(
  "microViz",
  repos = c(davidbarnett = "https://david-barnett.r-universe.dev", getOption("repos"))
)
BiocManager::install("phyloseq")
BiocManager::install("decontam")
devtools::install_github("tpq/propr")
BiocManager::install("ALDEx2")
BiocManager::install("EnhancedVolcano")
install.packages("vegan")
install.packages("repmis")
