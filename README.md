# [F.A.I.R.](https://www.go-fair.org/fair-principles/) and reproducible analysis for the paper: **"Exploring the data that explores the oceans: working towards robust eDNA workflows for ocean wildlife monitoring (submitted)"**
by Jessica R. Pearce 1, Philipp E. Bayer 1,2, Adam Bennett 1, Eric J. Raes 1,2, Marcelle E. Ayad 1, Shannon Corrigan 1,2, Matthew W. Fraser 1,2, Denise Anderson 3, Priscila Goncalves 1,2, Benjamin Callahan 4, Michael Bunce 5, Stephen Burnell 1,2, Sebastian Rauschert 1,2,*  

1 Minderoo Foundation, Perth 6000, WA  
2 The UWA Oceans Institute, The University of Western Australia, Crawley 6009, WA   
3 INSiGENe Pty Ltd.  
4 North Carolina State University, Raleigh, 27606, USA  
5 Department of Conservation, Wellington, New Zealand  

*Corresponding author  

Launch analysis: [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/MinderooFoundation/OceanOmics-amplicon-paper-analysis/HEAD?urlpath=rstudio)  
  
  
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/) [![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

## Analysis
This repository contains all data and code to generate the figures and statistics in the paper. Simply click on the above `binder` button to launch a Rstudio session in the browser, with access to all code and data in this GitHub repository. There, the code can interactively be changed and different plots and statistics can be (re-)created.

### What is binder?
For an overview of what binder is, please check out [this link](https://mybinder.org/).  

## Where does the data in this repo come from?
This repository contains the [`phyloseq`](https://joey711.github.io/phyloseq/index.html) objects for all three data sets analysed in the paper. The objects were generated with [Minderoo OceanOmics amplicon nextflow pipeline](https://github.com/MinderooFoundation/OceanOmics-amplicon-nf). Below is a detailed description, including code, to recreate the `phyloseq` objects.
The three data sets can be found here:  

- [West et al. 2021](https://doi.org/10.1111/ddi.13228): [North West Western Australia](https://datadryad.org/stash/dataset/doi:10.5061/dryad.8kprr4xmm)   
- [Minderoo OceanOmics](https://www.minderoo.org/oceanomics): [Cocos Keeling Island transect](https://www.ebi.ac.uk/ena/browser/view/PRJEB63982)  
- [Minderoo OceanOmics](https://www.minderoo.org/oceanomics): [Rowley Shoals Islands](https://www.ebi.ac.uk/ena/browser/view/PRJNA930913)  

Additionally, this repository includes a list of Australian marine fish species, named [Aust_fish_species_list.csv](https://github.com/MinderooFoundation/OceanOmics-amplicon-paper-analysis/blob/master/data/Aust_fish_species_list.csv). This was manually curated by domain experts, with data drawn from [Atlas of Living Australia](https://www.ala.org.au/) and the [Global Biodiversity Information Facility](https://www.gbif.org/).

The files in [metadata](https://github.com/MinderooFoundation/OceanOmics-amplicon-paper-analysis/tree/master/data/metadata) were generated as part of the data collection and sequencing, and are downloaded is part of the [downloading the data](https://github.com/MinderooFoundation/OceanOmics-amplicon-paper-analysis/tree/master/docs/get_data.md) description.  
Lastly, The [read_qc](https://github.com/MinderooFoundation/OceanOmics-amplicon-paper-analysis/tree/master/data/read_qc) folder contains QC output from the [`seqkit`](https://bioinf.shenwei.me/seqkit/) part of the nextflow pipeline and contains read QC statistics.

## Documentation: Generating the `phyloseq` objects

Everything in this section is optional and not required for re-analyzing the results of the paper. It is documented for full transparency and reproducibility of all results, should anyone desire to want to do so. Information on [setting up the compute environment](https://github.com/MinderooFoundation/OceanOmics-amplicon-paper-analysis/tree/master/docs/setup.md), [downloading the data](https://github.com/MinderooFoundation/OceanOmics-amplicon-paper-analysis/tree/master/docs/get_data.md) and [creating the `phyloseq` object](https://github.com/MinderooFoundation/OceanOmics-amplicon-paper-analysis/tree/master/docs/create_phyloseq.md) can be found in the `docs` folder or via the clickable links.

