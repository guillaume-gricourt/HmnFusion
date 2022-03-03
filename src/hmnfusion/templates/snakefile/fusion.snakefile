import tempfile

# Args.
FASTQ_FWD = config["fastq_fwd"]
FASTQ_REV = config["fastq_rev"]
BAM = config["bam"]
BED_GENEFUSE = config["bed_genefuse"]
BED_LUMPY = config["bed_lumpy"]
BED_HMNFUSION = config["bed_hmnfusion"]
NAME = config["name"]
OUTPUT = config["output_file"]

# Rules.
rule all:
	input:
		rules.quantification.output.vcf

rule extract:
	input:
		genefuse = GENEFUSE_FILE,
		lumpy = LUMPY
    params:
        genefuse = GENEFUSE_FMT
	output:
		json = temp(tempfile.NamedTemporaryFile(delete=False).name)
	shell:
		'hmnfusion extractfusion \
			--genefuse-"{params.genefuse}" "{input.genefuse}" \
			--lumpy-vcf "{input.lumpy}" \
			--output-json "{output.json}"'

rule quantification:
	input:
		json = rules.extract.output.json,
		bam = BAM,
		bed = BED_HMNFUSION
	params:
        name = NAME
	output:
		vcf = OUTPUT
	shell:
        'hmnfusion quantification \
			--hmnfusion-json "{input.json}" \
			--input-bam "{input.bam}" \
			--input-bed "{params.bed}" \
			--name "{params.name}" \
			--output-vcf "{output.vcf}"'
