# Analysis for the paper: "Exploring the data that explores the oceans: working towards robust eDNA workflows for ocean wildlife monitoring"
  

Launch analysis: [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/MinderooFoundation/OceanOmics-amplicon-paper-analysis/HEAD?urlpath=rstudio)

## Analysis
This repository contains all data and code to generate the figures and statistics in the paper. Simply click on the above binder button to launch a Rstudio session in the browser, with access to all code and data in this GitHub repository. There, the code can interactively be changed and different plots can be created.

### What is binder?
For an overview of what binder is, please check out [this link](https://mybinder.org/).  

## Where does the data in this repo come from?
This repository contains the [`phyloseq`](https://joey711.github.io/phyloseq/index.html) objects for all three data sets analysed in the paper. The objects were generated with [Minderoo OceanOmics amplicon nextflow pipeline](https://github.com/MinderooFoundation/OceanOmics-amplicon-nf). Below is a detailed description, including code, to recreate the `phyloseq` objects.
The three data sets can be found here:  

- [West et al. 2021](https://doi.org/10.1111/ddi.13228): [North West Western Australia](https://datadryad.org/stash/dataset/doi:10.5061/dryad.8kprr4xmm)   
- [Minderoo OceanOmics: Cocos Keeling Island transect:](https://www.ebi.ac.uk/ena/browser/view/PRJEB63982)  
- [Minderoo OceanOmics: Rowley Shoals Islands:](https://www.ebi.ac.uk/ena/browser/view/PRJNA930913)  

## Generating the `phyloseq` objects

Everything in this section is optional and not required for re-analyzing the results of the paper. It is documented for full transparency and reproducibility of all results, should anyone desire to want to so so.

### Setup compute environment
Unlike the easy to use setup with binder to analyse the `phyloseq` objects, the generation of the `phyloseq` objects themselves requires a bit more setup and likely access to an Apple or Linux machine, or a windows machine with [windows subsystem for linux](https://learn.microsoft.com/en-us/windows/wsl/install) installed. Also, free data storage of ~500Gb - 1TB is recommended to handle the `fastq` files. Below we will describe the setup for a Linux machine with a `Ubuntu 20.04` operating system. Apart from minoconda and nextflow, however, no other dependences are required, and our nextflow pipeline, in a fully reproducible manner, leverages [`docker`](https://www.docker.com/resources/what-container/#:~:text=A%20Docker%20container%20image%20is,tools%2C%20system%20libraries%20and%20settings.) container.

#### Install `miniconda`  
[`Miniconda`](https://docs.conda.io/en/latest/miniconda.html) is an application that allows to install software in a reproducible manner from an online repository. All this without the need for systems administration access or `root` access.  

Please follow the instructions for installing `miniconda` below, or on the [miniconda install instructions website](https://docs.conda.io/en/latest/miniconda-install.html).  

```zsh
# Download miniconda
cd ~ # go to home directory
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh # download miniconda installer 
bash Miniconda3-latest-Linux-x86_64.sh # execute the installer script and follow all instructions, including accepting the terms
```

The default install location will always be in the home directory, which can simply be accepted.  
Alternatively, you can also install the faster version of `conda`, [`mamba`](https://mamba.readthedocs.io/en/latest/mamba-installation.html#mamba-install)

#### Install nextflow  
`Nextflow` is very easy to install via the freshly installed `conda`/`mamba`. Alternatively, you can also install it via the [`nextflow` website instructions](https://www.nextflow.io/). The `conda` repository link for `nextflow` can be found [here](https://anaconda.org/bioconda/nextflow).

```zsh
# create a new conda environment for nextflow
conda create -n nextflow
conda activate nextflow  # we need to activate teh newly created environment, as per screen instruction
conda install -c bioconda nextflow # here we install nextflow into that environment, as per the instructions here https://anaconda.org/bioconda/nextflow
```

#### Install Docker  

This is highly recommended, even though the `nextflow` pipeline can be executed without it. `Docker` guarantees reproducibility through version control. 
Please follow the instruction [here to install `docker`](https://docs.docker.com/engine/install/ubuntu/) for Ubuntu, or any operating system you work off. This step will require sudo/root/admin access, so if you do not have that, consult with your IT team.  

In case you do have sudo/root/admin access:  
```zsh
# Update the apt package index and install packages to allow apt to use a repository over HTTPS
sudo apt-get update
sudo apt-get install ca-certificates curl gnupg

# Add Docker's official GPG key
sudo install -m 0755 -d /etc/apt/keyrings
curl -fsSL https://download.docker.com/linux/ubuntu/gpg | sudo gpg --dearmor -o /etc/apt/keyrings/docker.gpg
sudo chmod a+r /etc/apt/keyrings/docker.gpg

# Use the following command to set up the repository
echo \
  "deb [arch="$(dpkg --print-architecture)" signed-by=/etc/apt/keyrings/docker.gpg] https://download.docker.com/linux/ubuntu \
  "$(. /etc/os-release && echo "$VERSION_CODENAME")" stable" | \
  sudo tee /etc/apt/sources.list.d/docker.list > /dev/null
  
# Update the apt package index
sudo apt-get update

# Install the latest version of docker
sudo apt-get install docker-ce docker-ce-cli containerd.io docker-buildx-plugin docker-compose-plugin

# Verify it was succesfully installed
sudo docker run hello-world
```

And this is all the setup done.  

### Setup repositories

### Downloading the `fastq` files
To download the `fastq` files from the [European Nucleotide Archive (ENA)](https://www.ebi.ac.uk/ena/browser/), we used the [nf-core](https://nf-co.re/) pipeline [nf-core/fetchngs](https://nf-co.re/fetchngs/1.10.0).

### Cocos Island Transect data
```zsh
# Go to a directory where you have enough storage available
cd ~/analysis # this is an example, replace it with the path to the directory you want to place the files and analyse them
conda activate nextflow # if you haven't already  

# generate the required id.csv file with the ENA ID as per nf-core/fetchngs
echo "ERP149130" > ids.csv
nextflow run nf-core/fetchngs -profile docker --input ids.csv --outdir ./cocos
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
nextflow run nf-core/fetchngs -profile docker --input ids.csv --outdir ./rowley_shoals
```

### Nort West Western Australia  
This is a bit different, as both the metadata and the sample fastq files are located on DRyad rather than ENA or NCBI. So follow the below custom script to fetch the data.
```zsh
cd ~/analysis
mkdir nwwa
cd nwwa
mkdir metadata fastq
cd fastq
wget -O 16S_Fish_data.zip https://datadryad.org/stash/downloads/file_stream/1140656 
unzip 16S_Fish_data.zip
rm 16S_Fish_data.zip
mv  16S_Fish_data/* .
rm -rf 16S_Fish_data
### Setup & run amplicon pipeline
