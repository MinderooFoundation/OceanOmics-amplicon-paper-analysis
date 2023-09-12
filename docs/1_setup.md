# Setup

## Setup compute environment
Unlike the easy to use setup with binder to analyse the `phyloseq` objects, the generation of the `phyloseq` objects themselves requires a bit more setup and likely access to an Apple or Linux machine, or a windows machine with [windows subsystem for linux](https://learn.microsoft.com/en-us/windows/wsl/install) installed. Also, free data storage of ~500Gb - 1TB is recommended to handle the `fastq` files. Below we will describe the setup for a Linux machine with a `Ubuntu 20.04` operating system. Apart from minoconda and nextflow, however, no other dependencies are required, and our nextflow pipeline, in a fully reproducible manner, leverages [`docker`](https://www.docker.com/resources/what-container/#:~:text=A%20Docker%20container%20image%20is,tools%2C%20system%20libraries%20and%20settings.) container.

## Install `miniconda`  
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

## Install nextflow  
`Nextflow` is very easy to install via the freshly installed `conda`/`mamba`. Alternatively, you can also install it via the [`nextflow` website instructions](https://www.nextflow.io/). The `conda` repository link for `nextflow` can be found [here](https://anaconda.org/bioconda/nextflow).

```zsh
# create a new conda environment for nextflow
conda create -n nextflow
conda activate nextflow  # we need to activate teh newly created environment, as per screen instruction
conda install -c bioconda nextflow # here we install nextflow into that environment, as per the instructions here https://anaconda.org/bioconda/nextflow
```

## Install Docker (alternatively singularity)  

This is highly recommended, even though the `nextflow` pipeline can be executed without it. `Docker` guarantees reproducibility through version control. 
Please follow the instruction [here to install `docker`](https://docs.docker.com/engine/install/ubuntu/) for Ubuntu, or any operating system you work off. This step will require sudo/root/admin access, so if you do not have that, consult with your IT team.  

Alternatively, you can also use [singularity, with detailed installation isntructions here](https://docs.sylabs.io/guides/3.5/user-guide/quick_start.html)

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

## Set up blast nt database

For taxonomic annotation via the NCBI nt database, we need to set up the database on our machine first
```zsh
conda create -n blast
conda activate blast
conda install -c bioconda blast
cd ~/analysis
mkdir ncbi-nt
cd ncbi-nt
update_blastdb.pl
```

And this is all the setup done.  
