# HmnFusion

[![Github Version](https://img.shields.io/github/v/release/guillaume-gricourt/HmnFusion?display_name=tag&sort=semver)](version) [![Conda Release](https://img.shields.io/conda/vn/bioconda/hmnfusion.svg)](https://anaconda.org/bioconda/hmnfusion)  
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
[![GitHub Super-Linter](https://github.com/guillaume-gricourt/HmnFusion/workflows/Tests/badge.svg)](https://github.com/marketplace/actions/super-linter)
[![Coverage](https://img.shields.io/endpoint?url=https://gist.githubusercontent.com/guillaume-gricourt/5b62753442bc7c44ae2995299575af0a/raw/coverage.json)](code_coverage)  
[![DOI](https://zenodo.org/badge/259869577.svg)](https://zenodo.org/badge/latestdoi/259869577)  

## Getting Started

### Installing

The software can be installed with `conda` (best way), `pip` and `docker`.  
The commands `worfklow-hmnfusion` and `workflow-fusion` must be launched with the docker image `hmnfusion`.  
To construct the BAM files, the command `workflow-align` must be launched with the docker image `hmnfusion-align`.  

Install with `conda`
```bash
conda install -c bioconda hmnfusion
```

Install with `pip`
```bash
# Replace the <version> for the version considered
# Download from github
wget https://github.com/guillaume-gricourt/HmnFusion/releases/download/<version>/pip.zip
unzip pip.zip
# Install
pip install hmnfusion-<version>-py3-none-any.whl
# Clean up
rm pip.zip hmnfusion-<version>-py3-none-any.whl hmnfusion-<version>.tar.gz
```

Install with `docker`  
Specify `<version>` with digits, not using `latest`
```bash
# Pull the image hmnfusion
docker pull ghcr.io/guillaume-gricourt/hmnfusion:<version>
or
docker pull ggricourt/hmnfusion:<version>
# Pull the image hmnfusion-align
docker pull ggricourt/hmnfusion-align:<version>
```

### Running
Software is available by
<pre>hmnfusion <i>command</i> <i>options</i></pre>

## Commands

### Extract Fusion

Aggregate results from Genefuse and Lumpy to produce a Json file.

<pre>
hmnfusion extractfusion \
    --input-genefuse-json <i>file</i> | --input-genefuse-html <i>file</i>\
    --input-lumpy-vcf <i>file</i> \
    --output-hmnfusion-json <i>file</i>
</pre>

### Quantification

Calculate allelic frequency given postion or Json file produced by `extractfusion` command.
A fusion is defined by two breakpoints. Only one must be in bed intervals, allelic depth is computed only on this side.
Name sample is used in vcf file.

<pre>
hmnfusion quantification \
    --input-hmnfusion-json <i>file</i> | --region <i>chromosome:position</i> \
    --input-sample-bam <i>file</i> | --input-sam <i>file</i> \
    --input-hmnfusion-bed <i>file</i> \
    --name <i>sample_name</i> \
    --output-hmnfusion-vcf <i>file</i>
</pre>

### MMEJ - Deletion

Extract MMEJ information from VCF file coming from classic variant caller (GATK, Varscan, ...)

<pre>
hmnfusion mmej-deletion \
    --input-sample-vcf <i>file</i> <i>file</i> ... \
    --input-reference-fasta <i>file</i> \
    --output-hmnfusion-xlsx <i>file</i>
</pre>

### MMEJ - Fusion

Extract MMEJ information from fusion breakpoints.
<pre>
hmnfusion mmej-fusion \
    --input-hmnfusion-json <i>file</i> \
    --input-sample-bam <i>file</i> \
    --input-reference-fasta <i>file</i> \
    --name <i>sample_name</i> \
    --output-hmnfusion-xlsx <i>file</i>
</pre>

## Docker

### Run ExtractFusion
Run *extractfusion* with docker like:  
<pre>
docker run -it \
    --rm \
    hmnfusion:latest \
    extractfusion \
    --input-genefuse-json <i>file</i> | --input-genefuse-html <i>file</i> \
    --input-lumpy-vcf <i>file</i> \
    --output-hmnfusion-json <i>file</i>
</pre>

### Run Workflow HmnFusion
Running combined *extractfusion* and *quantification* with one command-line:  
<pre>
docker run -it \
    --rm \
    hmnfusion:latest \
    workflow-hmnfusion \
    --input-genefuse-html <i>file</i> \
    --input-lumpy-vcf <i>file</i> \
    --input-sample-bam <i>file</i> \
    --input-hmnfusion-bed <i>file</i> \
    --name <i>sample_name</i> \
    --output-hmnfusion-vcf <i>file</i>
</pre>

### Run Workflow Fusion
Run with one command-line GeneFuse, Lumpy and HmnFusion:  
<pre>
docker run -it \
    --rm \
    hmnfusion:latest \
    workflow-fusion \
    --input-forward-fastq <i>file</i> \
    --input-reverse-fastq <i>file</i> \
    --input-sample-bam <i>file</i> \
    --output-hmnfusion-vcf <i>file</i> \
    --output-genefuse-html <i>file</i> \
    --output-lumpy-vcf <i>file</i> \
    --input-reference-fasta <i>file</i> \
    --name <i>sample_name</i> \
    --threads 4
</pre>

### Build BAM files
Only the docker image hmnfusion-align could be use for this feature.  
Be aware, the size of the image is near to 15Gb.  
The reference files use to build BAM files could be cited with this DOI:  
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6619597.svg)](https://doi.org/10.5281/zenodo.6619597)  
<pre>
docker run -it \
    --rm \
    hmnfusion-align:latest \
    workflow-align \
    --input-forward-fastq <i>file</i> \
    --input-reverse-fastq <i>file</i> \
    --output-directory <i>file</i> \
    --threads 4
</pre>


## Versioning

[SemVer](http://semver.org/) is used for versioning.

## Authors

* **Guillaume Gricourt**  
* **Dr. Ivan Sloma**  
