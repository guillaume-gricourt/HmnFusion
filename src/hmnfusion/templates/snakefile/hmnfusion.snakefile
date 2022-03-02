import tempfile

print(config)
BAM = config["bam"]
LUMPY = config["lumpy"]
GENEFUSE = config["genefuse"]
HMNFUSION_QUANTIF = config["hmnfusion_quantification"]
BED = config["bed"]
OUTPUT = config["output"]

EXTRACT = tempfile.NamedTemporaryFile()
'''
rule all:
	input:
		"a"

rule extract:
	input:
		genefuse = GENEFUSE,
		lumpy = LUMPY
	output:
		json = temp(EXTRACT.name)
	shell:
		'hmnfusion extractfusion \
			--genefuse-html "{input.genefuse}" \
			--lumpy-vcf "{input.lumpy}" \
			--output-json "{output.json}"'

rule quantification:
	input:
		json = rules.extract.output.json,
		bam = BAM
	params:
		bed = BED
	output:
		vcf = OUTPUT
	shell:
        'hmnfusion quantification \
			--hmnfusion-json "{input.json}" \
			--input-bam "{input.bam}" \
			--input-bed {params.bed} \
			--name {wildcards.sfusion} \
			--output-vcf "{output.vcf}"'
'''
