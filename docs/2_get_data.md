# Get Data 

## Setup repositories

### Downloading the `fastq` files
To download the `fastq` files from the [European Nucleotide Archive (ENA)](https://www.ebi.ac.uk/ena/browser/), we used the [nf-core](https://nf-co.re/) pipeline [nf-core/fetchngs](https://nf-co.re/fetchngs/1.10.0).

### Cocos Island Transect data
```zsh
# Go to a directory where you have enough storage available
cd ~/analysis # this is an example, replace it with the path to the directory you want to place the files and analyse them
conda activate nextflow # if you haven't already  

# generate the required id.csv file with the ENA ID as per nf-core/fetchngs
echo "ERP149130" > ids.csv
nextflow run nf-core/fetchngs -r 1.10.0 -profile docker --input ids.csv --outdir ./cocos

# The data downloaded has both 12S and 16S samples. We just want the 16S samples
cat ./cocos/samplesheet/samplesheet.csv | grep 16S > ./cocos/samplesheet/samplesheet_16S.csv
```

### Rowley Shoals Islands
```zsh
# Go back to the analysis directory, where yo uhave downloaded the cocos data already
cd ~/analysis # this is an example, replace it with the path to the directory you want to place the files and analyse them
conda activate nextflow # if you haven't already  

# remove the existing 'ids.csv'
rm ids.csv

# generate the required id.csv file with the ENA ID as per nf-core/fetchngs
echo "SRP420753" > ids.csv

# Run the downlaod command again
nextflow run nf-core/fetchngs -r 1.10.0 -profile docker --input ids.csv --outdir ./rowley_shoals

# The data downloaded has both 12S and 16S samples. We just want the 16S samples
cat ./rowley_shoals/samplesheet/samplesheet.csv | grep 16S > ./rowley_shoals/samplesheet/samplesheet_16S.csv
```

### Nort West Western Australia  
This is a bit different, as both the metadata and the sample fastq files are located on Dryad rather than ENA or NCBI. So follow the below custom script to fetch the data.  

```zsh
cd ~/analysis
mkdir nwwa
cd nwwa
mkdir samplesheet fastq
cd fastq
wget -O 16S_Fish_data.zip https://datadryad.org/stash/downloads/file_stream/1140656 
unzip 16S_Fish_data.zip
rm 16S_Fish_data.zip
mv  16S_Fish_data/* .
rm -rf 16S_Fish_data
gzip *.fastq # files need to be gzipped for pipeline
```

#### Samplesheet

The metadata samplesheet will need to be created based on the sample names and information provided in the paper and supporting material. The below R code can be used to generate it, but the metadata.csv is also part of this repository (`data/metadata/NWWA_metadata.csv`)

```r
#-----------------------------------------------------------------------------
# Create NWWA metadata for nextflow pipeline to generate phyloseq object

library(phyloseq)
libary(tidyverse)

samples <- list()
samples$nw_false <- readRDS('data/phyloseq_objects/NWWA_16S_phyloseq_nt_FALSE.rds')


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


as_tibble(sample_data(samples$nw_false)) %>%
  mutate(sample() = rownames(sample_data(samples$nw_false))) %>%
  mutate(fastq_1 = paste0(sample, "~Fish16S.fastq.gz")) %>% 
  mutate(fastq_2 = "") %>%
  select(sample, fastq_1, fastq_2, assay, site, Control, Bioregion, Subregion, Habitat) %>%
  mutate(assay = as.character("16S")) %>%
  write_csv("data/metadata/NWWA_metadata.csv", na = "")

```


## Directory structure
Now that we have all the data downloaded, which is already demultiplexed on the sample level, we should have a directory structure as such:

```zsh
.
├── cocos
├── nw_wa
├── rowley_shoals
└── work
```

`work` was created by the `nf-core/fetchngs` run and can be ignored (contains information on executing the pipeline).
