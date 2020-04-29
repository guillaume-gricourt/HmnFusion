<h3 align=center>HmnFusion</h3>

# Introduction
A command-line toolkit to compute differents metrics about quality, identity-vigilance and coverage from high-throughput sequencing provided by targeted DNA.

# Getting Started

## Installing
Installation is done with `pip`
```bash
pip install hmngenomics
```

## Running
Software is available by
<pre>HmnGenomics <i>command</i> <i>options</i><pre>

# Quality Metrics

Quality metrics is computed with two steps : extract statistics to produce JSON file and rendering

## Raw statistics
Compute raw statistics from FASTQ, BAM or VCF files
<pre>
HmnGenomics quality \
    --name <i>NAME</i> \
    --output <i>JSON</i> \
    --fastq-forward-before <i>FASTQ</i> \
    --fastq-reverse-before <i>FASTQ</i> \
    --bam <i>BAM</i> \
    --bed <i>BED</i> \
    --vcf <i>VCF</i>
</pre>

## Rendering
Render one or more JSON files produced above to build PDF or HTML file
<pre>HmnGenomics renderer -i <i>NAME</i> -o <i>PDF</i></pre>

# Coverage Metrics

## Position not covered
Extract position not covered under customizable cut off
<pre>HmnGenomics depthmin -i <i>BAM</i> -b <i>BED</i> --cut-off 30 -o <i>XLXS</i></pre>

## Coverage of bed file
Compute statistics of coverage from a bed file
<pre>HmnGenomics depthtarget -i <i>BAM</i> -b <i>BED</i> -m target -o <i>XLSX</i></pre>

# Identity vigilance

## Infer sexe of samples
Infer sexe from BAM files and BED file to produce XLSX file.
<pre>HmnGenomics infersexe -i <i>BAM</i> -b <i>BED</i> -o <i>XLSX</i></pre>

## Extract SNPs
Extract SNPs in VCF file from BAM files.
<pre>HmnGenomics extractvcf -i <i>BAM</i> --vcf-reference <i>VCF</i> -o <i>XLSX</i></pre>

# Built with these main libraries

* [CNVkit](https://github.com/etal/cnvkit) - Powerful library
* [pysam](https://github.com/pysam-developers/pysam) - Essential library to work with BAM and VCF files
* [biopython](https://github.com/biopython/biopython) - Essential library to work with FASTQ files
* [Pandas](https://github.com/pandas-dev/pandas) - Essential dataframe object
* [Django](https://github.com/django/django) - Build html/pdf reports with template

# Versioning

[SemVer](http://semver.org/) is used for versioning.

# Authors

* **Guillaume Gricourt** - *Initial work*

# License

See the [LICENSE.md](LICENSE.md) file for details
