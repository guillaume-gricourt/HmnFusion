import tempfile

# Args.
LUMPY_FILE = config["lumpy_file"]
GENEFUSE_FILE = config["genefuse_file"]
GENEFUSE_FMT = config["genefuse_fmt"]
BAM_FILE = config["bam_file"]
BED_FILE = config["bed_file"]
NAME = config["name"]
OUTPUT_FILE = config["output_file"]

# Rules.
rule all:
	input:
		rules.quantification.output.vcf

rule genefuse:
    input:
        fwd = 
        rev =
        bed = 
    output:
        
rule lumpy:
    input:
        bam = 
    out
rule extract:
	input:
		genefuse = GENEFUSE_FILE,
		lumpy = LUMPY_FILE
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
		bam = BAM_FILE,
		bed = BED_FILE
	params:
        name = NAME
	output:
		vcf = OUTPUT_FILE
	shell:
        'hmnfusion quantification \
			--hmnfusion-json "{input.json}" \
			--input-bam "{input.bam}" \
			--input-bed "{params.bed}" \
			--name "{params.name}" \
			--output-vcf "{output.vcf}"'
