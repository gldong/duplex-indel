# Duplex-Indel

Duplex-Indel is a pipeline for somatic insertions and deletions (indel) calling from Tn5 transposase-based duplex sequencing data. It is inspired by and expands on the core methods for SNV calling in META-CS, [pre-pe](https://github.com/lh3/pre-pe/) and [lianti](https://github.com/lh3/lianti/). We thank Dr. Heng Li and other contributors to those repositories. Given the single-molecule resolution of duplex sequencing, Duplex-Indel can accomodate both single-cell and pooled-cell input, and is applicable to polyploid tissues. 

## Getting started

Duplex-Indel performs somatic indel calling by filtering out germline variants using matched bulk. Before running the pipeline, please make sure you have the following input data:

- `FASTQ` files of duplex sequencing reads (with strand barcodes in the reads)
- `BAM` file of matched bulk WGS data

## Installation

### 1. Create conda environment
If you don't have `conda` (or `miniconda`) already installed, see instructions [here](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html).

Install `mamba` in a `python3` base environment
```bash
conda create -n mybase python=3.7  # skip if you have a python3 environment
conda install -n mybase -c conda-forge mamba
conda activate mybase
```
Create a new environment for Duplex-Indel and install `snakemake` and `pysam`
```bash
mamba create -c conda-forge -c bioconda -n duplex-indel snakemake pysam
conda activate duplex-indel
```

### 2. Install Duplex-Indel
```bash
git clone https://github.com/gldong/duplex-indel.git
cd duplex-indel && make
```

### 3. Install other required tools
Duplex-Indel requires `seqtk`. 
```bash
# Install seqtk (under /path/to/duplex-indel)
git clone https://github.com/lh3/seqtk.git;
cd seqtk && make
cd ..
```

## Prepare reference files

Duplex-Indel supports both hg19 and hg38 alignments. Please make sure you are using the consistent version.

Some reference files are directly available in the `references` folder. For others, follow the instructions to download them. 

- **Reference genome regions**

	The genome is separated into smaller regions for parallel computing. We have provided region files for both [hg19](https://github.com/gldong/duplex-indel/blob/main/references/genome_regions_hg19.txt) and [hg38](https://github.com/gldong/duplex-indel/blob/main/references/genome_regions_hg38.txt).

- **Genome mask**

	High-quality regions of the genome for variant calling. We recommend using [UniMask](http://bit.ly/unimask) developed by Dr. Heng Li for [hg19](https://github.com/gldong/duplex-indel/blob/main/references/um75-hs37d5-comp.auto.bed) and easy regions (i.e. not in difficult regions) from Genome In A Bottle (GIAB) for hg38.
    ```bash
    # Download GIAB easy regions
    cd /path/to/duplex-indel/references
    wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.5/GRCh38@all/Union/GRCh38_notinalldifficultregions.bed.gz
    gunzip GRCh38_notinalldifficultregions.bed.gz
    ```

- **Reference genome**

	Make sure to use the same version as the bulk BAM. Generate necessary index files using the following command.
    ```bash
    cd /path/to/duplex-indel/references
    ```
    For hg19
    ```bash
    # Download hg19 with decoy
    wget ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz
    # Decompress
    gunzip hs37d5.fa.gz
    # Generate index files
    bwa index hs37d5.fa
    ```
    For hg38
    ```bash
    # Download hg38 without ALT contigs
    wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
    # Decompress
    gunzip GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
    # Change file extension to .fa
    mv GCA_000001405.15_GRCh38_no_alt_analysis_set.fna GCA_000001405.15_GRCh38_no_alt_analysis_set.fa
    # Generate index files
    bwa index GCA_000001405.15_GRCh38_no_alt_analysis_set.fa
    ```

- **Common indels**

	Download [gnomAD common indels](https://zenodo.org/records/15161320) (allele frequency >= 1%) to the `references` folder.

## Configure the run

A YAML file is required as the configuration for each run (i.e. each sample). We have provided example YAML files in the `scripts` folder for [hg19](https://github.com/gldong/duplex-indel/blob/main/scripts/config_hg19.yaml) and [hg38](https://github.com/gldong/duplex-indel/blob/main/scripts/config_hg38.yaml), which you can modify with your input file paths. 

## Run pipeline

```bash
conda activate duplex-indel
export PATH=/path/to/duplex-indel:/path/to/duplex-indel/seqtk:$PATH

# For single-cell input
bash run_pipeline.sh single SAMPLE_ID > LOG_DIR/log.SAMPLE_ID

# For pooled-cell input
bash run_pipeline.sh pooled SAMPLE_ID > LOG_DIR/log.SAMPLE_ID
```

## Cluster computing

In practice, processing duplex sequencing data requires cluster computing. An example using the SLURM system is provided here. 

### Cluster configuration

You can find example cluster configurations in the `scripts` folder for [single-cell](https://github.com/gldong/duplex-indel/tree/main/scripts/slurm.single_cell) and [pooled-cell](https://github.com/gldong/duplex-indel/tree/main/scripts/slurm.pooled_cell) input. 

### Run pipeline (sbatch)
```bash
conda activate duplex-indel
export PATH=/path/to/duplex-indel:/path/to/duplex-indel/seqtk:$PATH

# For single-cell input
sbatch -o LOG_DIR/SAMPLE_ID.%j run_pipeline.sh single SAMPLE_ID

# For pooled-cell input
sbatch -o LOG_DIR/SAMPLE_ID.%j run_pipeline.sh pooled SAMPLE_ID

# Note that "%j" is the job ID assigned on the cluster.
```

### [OPTIONAL] Apply post-calling filters

We also provide a few optional filters which can be customized for specific use cases, including:

1) Additional filtering of variants confounded by read merging using unmerged BAM [-m]

2) Filtering by unique Tn5 sites present at each variant [-t]

3) Filtering by unique barcode pairs covering each variant [-b] 

4) Additional filtering of variants based on their position towards the end of the read (default = 10) [-e <int>]

To apply one or multiple of these filters, use the following command:
```bash
bash apply_filters.sh -s SAMPLE_ID -c CONFIG_FILE -d SAMPLE_DIR [-m] [-t] [-b] [-e <int>]
```

Example:
```bash
# Apply filters for reading merging and variant position
bash apply_filters.sh -s SAMPLE_ID -c CONFIG_FILE -d SAMPLE_DIR -m -e 15
```




