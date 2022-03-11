# HmnFusion

[![Release](https://img.shields.io/endpoint?url=https://gist.githubusercontent.com/guillaume-gricourt/5b62753442bc7c44ae2995299575af0a/raw/version.json)](version)  
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
[![GitHub Super-Linter](https://github.com/brsynth/rpFbaAnalysis/workflows/Tests/badge.svg)](https://github.com/marketplace/actions/super-linter)
[![Coverage](https://img.shields.io/endpoint?url=https://gist.githubusercontent.com/guillaume-gricourt/5b62753442bc7c44ae2995299575af0a/raw/coverage.json)](code_coverage)  

## Introduction
A tool to aggregate results of fusion produced by Genefuse and Lumpy and calculate allelic frequency.  

## Getting Started

### Installing
Installation is done with `pip`
```bash
pip install hmnfusion
```

### Running
Software is available by
<pre>HmnFusion <i>command</i> <i>options</i></pre>

## Commands

### Extract Fusion

Aggregate results from Genefuse and Lumpy to produce a Json file.

<pre>
HmnFusion extractfusion \
    --genefuse-json <i>file</i> | --genefuse-html <i>file</i>\
    --lumpy-vcf <i>file</i> \
    --output-json <i>file</i>
</pre>

### Quantification

Calculate allelic frequency given postion or Json file produced by `extractfusion` command.
A fusion is defined by two breakpoints. Only one must be in bed intervals, allelic depth is computed only on this side.
Name sample is used in vcf file.

<pre>
HmnFusion quantification \
    --hmnfusion-file <i>file</i> | --region <i>chromosome:position</i> \
    --input-bam <i>file</i> | --input-sam <i>file</i> \
    --input-bed <i>file</i> \
    --name <i>sample_name</i> \
    --output-vcf <i>file</i>
</pre>

## Docker

### Run ExtractFusion
Run *extractfusion* with docker like:  
<pre>
docker run -it \
    --rm \
    hmnfusion:latest \
    extractfusion \
    --genefuse-json <i>file</i> | --genefuse-html <i>file</i> \
    --lumpy-vcf <i>file</i> \
    --output-json <i>file</i>
</pre>

### Run Workflow HmnFusion
Running combined *extractfusion* and *quantification* with one commandi line:   
<pre>
docker run -it \
    --rm \
    hmnfusion:latest \
    workflow-hmnfusion \
    --genefuse-html <i>file</i> \
    --lumpy-vcf <i>file</i> \
    --input-bam <i>file</i> \
    --input-bed <i>file</i> \
    --name <i>sample_name</i> \
    --output-vcf <i>file</i>
</pre>

### Run Workflow Fusion
Run with one command line GeneFuse, Lumpy and HmnFusion:  
<pre>
docker run -it \
    --rm \
    hmnfusion:latest \
    workflow-fusion \
    --input-fastq-forward <i>file</i> \
    --input-fastq-reverse <i>file</i> \
    --input-bam <i>file</i> \
    --output-vcf <i>file</i> \
    --input-reference-fasta <i>file</i> \
    --input-bed-hmnfusion <i>file</i> \
    --name <i>sample_name</i> \
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

* **Guillaume Gricourt** - *Initial work*

## License

See the [LICENSE](LICENSE) file for details
