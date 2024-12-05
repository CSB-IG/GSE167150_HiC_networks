# hic_networks_analysis

Analysis of TNBC Hi-C Networks

## HiEdge install

Create conda environment for HiEdge

```
conda env create -f workflow/envs/hiedge.yaml
```
Obtain HiEdge

```
mkdir -p resources/software

git clone git@github.com:Gabrielstav/HiEdge.git resources/software/HiEdge

```

## HiEdge prepare 

### ENCODE blacklist 

Obtain resources to run HiEdge. Specifically the [ENCODE blacklisted regions](https://github.com/Boyle-Lab/Blacklist?tab=readme-ov-file) of the genome. These regions are problematic because they may draw a lot of reads during alignment, artificially skew the counts and lead to inaccurate subsequent analyses (e.g. thresholding, normalizing, calling peaks, etc). The existence of these regions is related to the assembly of reference genome sequences and their treatment of repetitive regions.

This is [the link](https://www.encodeproject.org/annotations/ENCSR636HFF) to the blacklist in ENCODE.

Download blacklist bed file and set it up in resources

```
mkdir -p resources/blacklist

wget https://www.encodeproject.org/files/ENCFF356LFX/@@download/ENCFF356LFX.bed.gz -P resources/blacklist/

# decompress bed file
gunzip resources/blacklist/ENCFF356LFX.bed.gz

```

### Cytoband regions

The cytoband coordinates can be provided to HiEdge to filter out interactions involving centromeric regions. Centromeres are typically filtered out from Hi-C data because their repetitive sequence makes it difficult to map reads to those regions, so the interactions can be unreliable. People usually remove them rather than correcting for the potential technical biases.

Telomeres are probably filtered out too, but I don't recall Gabriel mentioning it in his thesis or the github page for HiEdge. Anywa, because we start from raw matrices to create the hi-c networks, we must filter out the "problematic" regions first (just as one would do before matrix balancing).

The cytoband positions are obtained as a bed file from the UCSC genome annotation database ([here](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/) is the dump for GRCh38).

```
mkdir -p resources/cytobands

wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cytoBand.txt.gz -P resources/cytobands/

# decompress the file
gunzip resources/cytobands/cytoBand.txt.gz

# rename for clarity
mv resources/cytobands/cytoBand.txt resources/cytobands/cytoBand_hg38.txt

```

## HiEdge config file

Set up HiEdge config file. We should do a jinja template and have a decent file structure to do this, but we'll just use this for now

```
#mkdir -p resources/hiedge

#cp resources/software/HiEdge/config/config.yaml resources/hiedge/

# edit the config file to set it up for this data

```

## HiEdge run


```
conda activate hiedge

python resources/HiEdge/main.py -c path/to/config_file.yaml -r run_name
```




