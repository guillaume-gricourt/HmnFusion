Usage
=====

HmnFusion is designed to analyze some panel of genes sequenced with the aCAP-Seq library and performs with a high sensitivity and specificity `10.1016/j.jmoldx.2022.07.004 <https://www.sciencedirect.com/science/article/pii/S1525157822002185?via%3Dihub>`_

Steps
-----

1. BAM files are analyzed by `Lumpy <https://github.com/arq5x/lumpy-sv>`_ and `Genefuse <https://github.com/OpenGene/GeneFuse>`_
2. The results are aggregated by HmnFusion
3. Now you can calculate the fusion frequency and you can define the MMEJ sequence

These steps involved one or more commands.
A workflow is available to grouped these commands into one command.


How to start with HmnFusion
---------------------------

It depends what you have:

* If you have FASTQ files, you need first to create BAM files, please go to `Create BAM files`_
* If you have FASTQ, BAM files, please go to `Fusion callers`_
* If you have FASTQ, BAM files, HTML of JSON from Genefuse and VCF from Lumpy, please go to `Fusion frequency`_ or `MMEJ sequences`_

Questions about extra files: hg19 reference, BED files ? please read the :doc:`FAQ <faq>`.

Create BAM files
----------------

It's only available through ``docker``

.. code-block:: console

    $ docker run -it \
        --rm \
        hmnfusion-align:latest \
        workflow-align \
        --input-forward-fastq <FASTQ forward, file> \
        --input-reverse-fastq <FASTQ reverse, file> \
        --output-directory <Output directory> \
        --threads 4

Fusion callers
--------------

It will run Genefuse, Lumpy and HmnFusion, to detect and quantify fusions.

.. code-block:: console

    $ hmnfusion workflow-fusion \
        # Sample
        --name <Name of sample> \
        --input-forward-fastq <Fastq file forward> \
        --input-reverse-fastq <Fastq file reverse> \
        --input-sample-bam <Bam file> \
        # Bed
        --input-genefuse-bed <Genefuse bed file> \
        --input-lumpy-bed <Lumpy bed file> \
        --input-hmnfusion-bed <HmnFusion bed file> \
        # Reference
        --input-reference-fasta <Reference fasta file (hg19)> \
        # Output
        --output-hmnfusion-vcf <Vcf file output> \
        --output-genefuse-html <Genefuse html file output> \
        --output-lumpy-vcf <Lumpy vcf file output> \
        --threads [1-6]


Fusion frequency
----------------

It will extract fusions from Genefuse and Lumpy to quantify them.

.. code-block:: console

    $ hmnfusion workflow-fusion \
        # Sample
        --input-lumpy-vcf <Lumpy Vcf file> \
        --input-genefuse-json <Genefuse, json file> OR --input-genefuse-html <Genefuse, html file> \
        --input-sample-bam <Bam file> \
        --name <Name of sample> \
        # Bed
        --input-hmnfusion-bed <HmnFusion bed file> \
        # Output
        --output-hmnfusion-vcf <Vcf file output>

MMEJ sequences
--------------

Define fusions of interest.

.. code-block:: console

    $ hmnfusion extractfusion \
        # Sample
        --input-genefuse-json <Genefuse, json file> \
        --input-genefuse-html <Genefuse, html file> \
        --input-lumpy-vcf <Lumpy vcf file> \
        # Bed
        --input-hmnfusion-bed <Bed file> \
        # Output
        --output-hmnfusion-json <Json file output>

Extract MMEJ sequences.

.. code-block:: console

    $ hmnfusion mmej-fusion \
        # Sample
        --input-hmnfusion-json <HmnFusion, json file> \
        --input-sample-bam <Bam file> \
        --name <Name of sample> \
        # References
        --input-reference-fasta <Reference, fasta file> \
        # Output
        --output-hmnfusion-xlsx <Excel file output> \
        --output-hmnfusion-json <Json file output>
