# HmnFusion

[![Release](https://img.shields.io/endpoint?url=https://gist.githubusercontent.com/guillaume-gricourt/5b62753442bc7c44ae2995299575af0a/raw/version.json)](version)  
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
[![GitHub Super-Linter](https://github.com/guillaume-gricourt/HmnFusion/workflows/Tests/badge.svg)](https://github.com/marketplace/actions/super-linter)
[![Coverage](https://img.shields.io/endpoint?url=https://gist.githubusercontent.com/guillaume-gricourt/5b62753442bc7c44ae2995299575af0a/raw/coverage.json)](code_coverage)  
[![DOI](https://zenodo.org/badge/259869577.svg)](https://zenodo.org/badge/latestdoi/259869577)  

## Getting Started

### Installing

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
```bash
# Install HmnFusion
docker pull ghcr.io/guillaume-gricourt/hmnfusion:<version>
or
docker pull ggricourt/hmnfusion:<version>
# Install HmnFusion-Align
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
The reference files use to build BAM files could be cite with this DOI:  
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


## Built with these main libraries

* [beautifulsoup4](https://pypi.org/project/beautifulsoup4) - Parsing efficiently HTML file
* [pysam](https://github.com/pysam-developers/pysam) - Essential library to work with BAM and VCF files
* [Pandas](https://github.com/pandas-dev/pandas) - Essential dataframe object
* [Snakemake](https://snakemake.readthedocs.io/en/stable/) - Use workflow

## Versioning

[SemVer](http://semver.org/) is used for versioning.

## Authors

* **Guillaume Gricourt**  
* **Dr. Ivan Sloma**  

## License

See the [LICENSE](LICENSE) file for details
