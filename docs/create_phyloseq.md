# Create Phyloseq wirh `OceanOmics-amplicon-nf`

 
## Download the `MinderooFoundation/OceanOmics-amplicon-nf`

```zsh
git clone git@github.com:MinderooFoundation/OceanOmics-amplicon-nf.git
```

### Execute the pipeline

#### Cocos Keeling Transect   

```zsh
conda activate nextflow
cd ~/data

nextflow run OceanOmics-amplicon-nf main.nf \
                --input ./cocos/samplesheet/samplesheet.csv \
                --outdir ./cocos/amplicon_analysed \
                --dbfiles /data/tools/databases/ncbi-nt/ \ # If you want to blast against NCBI nt database, you need to first download it to your machine
                --bind_dir /data \
                -profile docker \
                --skip_demux true
```

#### Rowley Shoals Islands  

```zsh
conda activate nextflow
cd ~/data

nextflow run OceanOmics-amplicon-nf main.nf \
                --input ./rowley_shoals/samplesheet/samplesheet.csv \
                --outdir ./rowley_shoals/amplicon_analysed \
                --dbfiles /data/tools/databases/ncbi-nt/ \ # If you want to blast against NCBI nt database, you need to first download it to your machine
                --bind_dir /data \
                -profile docker \
                --skip_demux true
```