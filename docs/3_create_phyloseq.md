# Create Phyloseq wirh `OceanOmics-amplicon-nf`


## Execute the pipeline

Note: the paper contains three phyloseq objects for every data set. Those differ in the `pool` option of the `dada` function in `DADA2`. Run the nextflow pipeline three times, and change the `-pooled` parameter from `-pooled true` to `-pooled false` and then `-pooled pseudo`

Also: the analysis needs to be started in the directory with the fastq files, so follow the `cd` command in the code snippets. This code assumes the same directory setup as based on the `2_get_data` setup steps.

### Cocos Keeling Transect   

```zsh
conda activate nextflow
cd ~/analysis/cocos/fastq

nextflow run MinderooFoundation/OceanOmics-amplicon-nf main.nf \
                --input ../samplesheet/samplesheet.csv \
                --outdir ../amplicon_analysed \
                --dbfiles "/data/tools/databases/ncbi-nt/*" \ # If you want to blast against NCBI nt database, you need to first download it to your machine
                --filter_table https://raw.githubusercontent.com/a4000/test_data/main/other_files/filter.csv
                --bind_dir /data \
                -profile docker \
                --skip_demux true \
                -pooled true  # this is the option that needs to be changed to false and pseudo and re-run
```

### Rowley Shoals Islands  

```zsh
conda activate nextflow
cd ~/analysis/rowley_shoals/fastq

nextflow run MinderooFoundation/OceanOmics-amplicon-nf main.nf \
                --input ../samplesheet/samplesheet.csv \
                --outdir ../amplicon_analysed \
                --dbfiles "/data/tools/databases/ncbi-nt/*" \ # If you want to blast against NCBI nt database, you need to first download it to your machine
                --filter_table https://raw.githubusercontent.com/a4000/test_data/main/other_files/filter.csv
                --bind_dir /data \
                -profile docker \
                --skip_demux true \
                -pooled true  # this is the option that needs to be changed to false and pseudo and re-run
```

### North West Western Australia  

```zsh
conda activate nextflow
cd ~/analysis/nwwa

nextflow run MinderooFoundation/OceanOmics-amplicon-nf main.nf \
                --input ../samplesheet/NWWA_metadata.csv \
                --outdir ../amplicon_analysed \
                --dbfiles "/data/tools/databases/ncbi-nt/*" \ # If you want to blast against NCBI nt database, you need to first download it to your machine
                --filter_table https://raw.githubusercontent.com/a4000/test_data/main/other_files/filter.csv
                --bind_dir /data \
                -profile docker \
                --skip_demux true \
                -pooled true  # this is the option that needs to be changed to false and pseudo and re-run
```

## Optionally: Download the `MinderooFoundation/OceanOmics-amplicon-nf`

There is no need to download the pipeline from GitHub, as the execution command will automatically retrieve all necessary files; in case you are 
working offline, you can download the repository up front. Just make sure you reference the right path to the `main.nf` file in this repository and point to the files you need to analyse.

```zsh
git clone git@github.com:MinderooFoundation/OceanOmics-amplicon-nf.git
```
