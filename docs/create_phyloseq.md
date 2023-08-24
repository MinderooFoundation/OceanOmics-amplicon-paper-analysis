# Create Phyloseq wirh `OceanOmics-amplicon-nf`

 
## Download the `MinderooFoundation/OceanOmics-amplicon-nf`

```zsh
git clone git@github.com:MinderooFoundation/OceanOmics-amplicon-nf.git
```

### Execute the pipeline

Note: the paper contains three phyloseq objects for every data set. Those differ in the `pool` option of the `dada` function in `DADA2`. Run the nextflow pipeline three times, and change the `-pooled` parameter from `-pooled true` to `-pooled false` and then `-pooled pseudo`

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
                --skip_demux true \
                -pooled true \
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
                --skip_demux true \
                -pooled true \
```

#### North West Western Australia  

```zsh
conda activate nextflow
cd ~/data

nextflow run OceanOmics-amplicon-nf main.nf \
                --input ./nw_wa/samplesheet/samplesheet.csv \
                --outdir ./rowleynw_wa_shoals/amplicon_analysed \
                --dbfiles /data/tools/databases/ncbi-nt/ \ # If you want to blast against NCBI nt database, you need to first download it to your machine
                --bind_dir /data \
                -profile docker \
                --skip_demux true \
                -pooled true \
```