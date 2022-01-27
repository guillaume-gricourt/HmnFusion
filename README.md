# HmnFusion

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


## Built with these main libraries

* [beautifulsoup4](https://pypi.org/project/beautifulsoup4) - Parsing efficiently HTML file
* [pysam](https://github.com/pysam-developers/pysam) - Essential library to work with BAM and VCF files
* [Pandas](https://github.com/pandas-dev/pandas) - Essential dataframe object

## Versioning

[SemVer](http://semver.org/) is used for versioning.

## Authors

* **Guillaume Gricourt** - *Initial work*

## License

See the [LICENSE](LICENSE) file for details
