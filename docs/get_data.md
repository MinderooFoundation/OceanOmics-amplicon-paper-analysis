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
This is a bit different, as both the metadata and the sample fastq files are located on Dryad rather than ENA or NCBI. So follow the below custom script to fetch the data.  

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
