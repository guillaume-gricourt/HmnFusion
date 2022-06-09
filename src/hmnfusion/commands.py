import argparse
import logging
import math
import os
import tempfile

from hmnfusion import _version
from hmnfusion import bed as ibed
from hmnfusion import (
    download_reference,
    extractfusion,
    fusion,
    fusion_flag,
    graph,
    install_software,
    mmej_deletion,
    mmej_fusion,
    quantification,
    region,
    utils,
    workflow,
)

AP = argparse.ArgumentParser(
    description="Extraction sequence and quantification of fusions \
        from DNA genomic analysis",
    epilog="See online documentation: ",
)
AP_subparsers = AP.add_subparsers(help="Sub-commnands (use with -h for more info)")


# Extract fusion.
def _cmd_extract_fusion(args):
    """Extract fusion from Genefuse and Lumpy analysis"""
    logging.info("Start analysis")
    # Grep args.
    finputs = {}
    finputs["genefuse"] = {}
    if args.input_genefuse_json:
        finputs["genefuse"]["path"] = os.path.abspath(args.input_genefuse_json)
        finputs["genefuse"]["format"] = "json"
    elif args.input_genefuse_html:
        finputs["genefuse"]["path"] = os.path.abspath(args.input_genefuse_html)
        finputs["genefuse"]["format"] = "html"
    finputs["lumpy"] = {}
    finputs["lumpy"]["path"] = os.path.abspath(args.input_lumpy_vcf)
    finputs["lumpy"]["format"] = "vcf"
    foutput = args.output_hmnfusion_json
    # Check if all exists.
    if not os.path.isfile(finputs["genefuse"]["path"]):
        utils.abort(
            AP, "File Genefuse doesn't exist : %s" % (finputs["genefuse"]["path"],)
        )
    if not os.path.isfile(finputs["lumpy"]["path"]):
        utils.abort(AP, "File Lumpy doesn't exist : %s" % (finputs["lumpy"]["path"],))
    if not os.path.isdir(os.path.dirname(os.path.abspath(foutput))):
        utils.abort(AP, "Outdir doesn't exist : %s" % (foutput,))
    bed_hmnfusion = args.input_hmnfusion_bed
    if bed_hmnfusion is None or not os.path.isfile(bed_hmnfusion):
        logging.warning("Use default hmnfusion bed")
        bed_hmnfusion = os.path.join(
            os.path.dirname(os.path.realpath(__file__)),
            "templates",
            "bed",
            "hmnfusion.bed",
        )

    # Parsing bed file.
    logging.info("Parsing bed file")
    bed = ibed.Bed.from_bed(bed_hmnfusion)

    # Run.
    g = graph.Graph(args.consensus_interval)
    # Genefuse.
    logging.info("Parse Genefuse")
    extractfusion.read_genefuse(
        g, finputs["genefuse"]["path"], finputs["genefuse"]["format"]
    )

    # Lumpy.
    logging.info("Parse Lumpy")
    extractfusion.read_lumpy(g, finputs["lumpy"]["path"], finputs["lumpy"]["format"])

    # Consensus.
    logging.info(
        "Build consensus with interval of %s pb"
        % (g.graph.graph["consensus_interval"],)
    )
    g.consensus_single()
    g.consensus_genefuse_lumpy()

    # Define nodes of interest.
    logging.info("Define nodes of interest")
    g.define_node_interest()

    # Filter same chrom.
    logging.info("Filter nodes")
    g.trim_node(bed)

    # Write output.
    if foutput:
        logging.info("Write output")
        g.to_json(path=foutput, metadata=dict(finputs=finputs, step="extractfusion"))
    logging.info("Analysis is finished")


P_extract_fusion = AP_subparsers.add_parser(
    "extractfusion", help=_cmd_extract_fusion.__doc__
)
P_efgg = P_extract_fusion.add_mutually_exclusive_group(required=True)
P_efgg.add_argument("--input-genefuse-json", help="Genefuse, json file")
P_efgg.add_argument("--input-genefuse-html", help="Genefuse, html file")
P_extract_fusion.add_argument("--input-lumpy-vcf", required=True, help="Lumpy vcf file")
P_extract_fusion.add_argument("--input-hmnfusion-bed", help="Bed file")
P_extract_fusion.add_argument(
    "--consensus-interval",
    type=int,
    default=500,
    help="Interval, pb, for which Fusion are considered equal if their chrom are",
)
P_extract_fusion.add_argument("--output-hmnfusion-json", help="Json file output")
P_extract_fusion.set_defaults(func=_cmd_extract_fusion)


# Quantification of fusions
def _cmd_quantification(args):
    """Quantify fusion candidates with a BAM file"""
    logging.info("Start analysis")

    # Check if all exists.
    logging.info("Check args")
    if args.input_hmnfusion_json and not os.path.isfile(args.input_hmnfusion_json):
        utils.abort(
            AP,
            "HmnFusion Json file doesn't exist : %s" % (args.input_hmnfusion_json,),
        )
    if args.region and not region.Region.check_region(args.region):
        utils.abort(
            AP,
            "Region format is not well formated. Required <chrom>:<position>",
        )
    if not os.path.isfile(args.input_sample_bam):
        utils.abort(
            AP,
            "Input alignment file doesn't exist : %s" % (args.input_sample_bam,),
        )
    if not utils.check_bam_index(args.input_sample_bam):
        utils.abort(
            AP, "Input alignment file must be in BAM format, index could not be build"
        )
    bed_hmnfusion = args.input_hmnfusion_bed
    if bed_hmnfusion is None or not os.path.isfile(bed_hmnfusion):
        logging.warning("Use default hmnfusion bed")
        bed_hmnfusion = os.path.join(
            os.path.dirname(os.path.realpath(__file__)),
            "templates",
            "bed",
            "hmnfusion.bed",
        )
    if args.output_hmnfusion_vcf:
        if not os.path.isdir(
            os.path.dirname(os.path.abspath(args.output_hmnfusion_vcf))
        ):
            utils.abort(AP, "Outdir doesn't exists: %s" % (args.output_hmnfusion_vcf,))
    if args.output_hmnfusion_json:
        if not os.path.isdir(
            os.path.dirname(os.path.abspath(args.output_hmnfusion_json))
        ):
            utils.abort(AP, "Outdir doesn't exists: %s" % (args.output_hmnfusion_json,))

    # Init.
    params = dict(
        falignment=dict(path=args.input_sample_bam, mode="rb"),
        clipped=dict(count=args.baseclipped_count, interval=args.baseclipped_interval),
    )

    if params["clipped"]["count"] < 1:
        logging.warning("Parameter base clipped count is too low, is set to 1 pb")
        params["clipped"]["count"] = 1
    elif params["clipped"]["count"] > 20:
        logging.warning("Parameter base clipped count is too high, is set to 20 pb")
        params["clipped"]["count"] = 20
    if params["clipped"]["interval"] < 1:
        logging.warning("Parameter base clipped interval is too low, is set to 1 pb")
        params["clipped"]["interval"] = 1
    elif params["clipped"]["interval"] > 20:
        logging.warning("Parameter base clipped interval is too high, is set to 20 pb")
        params["clipped"]["interval"] = 20
    if params["clipped"]["count"] > params["clipped"]["interval"]:
        logging.warning(
            "Parameter base clipped count is higher than parameter base \
            clipped interval -> count is set equal than interval : %s pb \
            instead of %s pb"
            % (params["clipped"]["interval"], params["clipped"]["count"])
        )
        params["clipped"]["count"] = params["clipped"]["interval"]

    if not utils.validate_name_sample(args.name):
        logging.warning("Name sample is not valid: %s" % (args.name,))

    # Parsing bed file.
    logging.info("Parsing bed file")
    bed = ibed.Bed.from_bed(bed_hmnfusion)

    # Parsing fusions.
    logging.info("Get region")
    g = graph.Graph()
    if args.region:
        fus = fusion.Fusion()
        r = region.Region.from_str(args.region)
        fus.set_region(r)
        fus.evidence.raw = 0
        g.add_node(fus, 0, False, True)
    elif args.input_hmnfusion_json:
        g = graph.Graph.from_json(args.input_hmnfusion_json)

    # Process
    logging.info("Calcul VAF fusion")
    quantification.run(params, bed, g)

    # Write output.
    if args.output_hmnfusion_vcf:
        logging.info("Write output vcf")
        quantification.write(args.output_hmnfusion_vcf, args.name, g)
    if args.output_hmnfusion_json:
        logging.info("Write output json")
        g.to_json(
            path=args.output_hmnfusion_json,
            metadata=dict(name=args.name, step="quantification"),
        )

    logging.info("Analysis is finished")


P_quantification = AP_subparsers.add_parser(
    "quantification", help=_cmd_quantification.__doc__
)
P_qpg = P_quantification.add_mutually_exclusive_group(required=True)
P_qpg.add_argument("--region", help="Region format <chrom>:<postion>")
P_qpg.add_argument(
    "--input-hmnfusion-json", help='Output Json produced by command "extractfusion"'
)
P_quantification.add_argument("--input-sample-bam", required=True, help="Bam file")
P_quantification.add_argument("--input-hmnfusion-bed", help="Bed file")
P_quantification.add_argument("--name", required=True, help="Name of sample")
P_quantification.add_argument(
    "--baseclipped-interval",
    type=int,
    default=6,
    help="Interval to count hard/soft-clipped bases from fusion point (pb)",
)
P_quantification.add_argument(
    "--baseclipped-count",
    type=int,
    default=4,
    help="Number of base hard/soft-clipped bases to count in interval (pb)",
)
P_quantification.add_argument("--output-hmnfusion-vcf", help="Vcf file output")
P_quantification.add_argument("--output-hmnfusion-json", help="Json file output")
P_quantification.set_defaults(func=_cmd_quantification)


# mmej signatures.
def _cmd_mmej_deletion(args):
    """Detect MMEJ signatures from deletion, fusion"""
    logging.info("Start analysis")
    # Grep args.
    for finput in args.input_sample_vcf:
        if not os.path.isfile(finput):
            utils.abort(AP, 'Vcf file doesn"t exist: %s' % (finput,))
    if not os.path.isfile(args.input_reference_fasta):
        utils.abort(
            AP, 'Reference file doesn"t exist: %s' % (args.input_reference_fasta,)
        )
    if not os.path.isdir(os.path.dirname(os.path.abspath(args.output_hmnfusion_xlsx))):
        utils.abort(AP, 'Outdir doesn"t exist: %s' % (args.output_hmnfusion_xlsx,))

    # Run.
    logging.info("Check index fasta reference")
    if not utils.check_fasta_index(args.input_reference_fasta):
        utils.abort(
            'Index of fasta file doesn"t exist: %s' % (args.input_reference_fasta,)
        )

    logging.info("Extract events from files")
    mmej_deletions = []
    for vcf_file in args.input_sample_vcf:
        mmej_deletions.extend(mmej_deletion.MmejDeletion.from_vcf(path=vcf_file))

    logging.info("Caracterize events with reference file")
    for mmej_del in mmej_deletions:
        mmej_del.set_value_sequence(path=args.input_reference_fasta)

    logging.info("Write output")
    mmej_deletion.MmejDeletion.to_excel(
        path=args.output_hmnfusion_xlsx,
        mmej_deletions=mmej_deletions,
    )
    logging.info("Analysis is finished")


P_mmej_deletion = AP_subparsers.add_parser(
    "mmej-deletion", help=_cmd_mmej_deletion.__doc__
)
P_mmej_deletion.add_argument("--input-sample-vcf", nargs="+", help="Vcf file")
P_mmej_deletion.add_argument(
    "--input-reference-fasta", required=True, help="Genome of reference"
)
P_mmej_deletion.add_argument("--output-hmnfusion-xlsx", help="Output file")
P_mmej_deletion.set_defaults(func=_cmd_mmej_deletion)


def _cmd_mmej_fusion(args):
    """Identify MMEJ from fusion"""
    logging.info("Start analysis - MMEJ Fusion")
    # Check if all exists.
    if not os.path.isfile(args.input_hmnfusion_json):
        utils.abort(AP, "File input doesn't exist: %s" % (args.input_hmnfusion_json,))
    if not os.path.isfile(args.input_reference_fasta):
        utils.abort(
            AP, "File reference doesn't exist: %s" % (args.input_reference_fasta,)
        )
    if not os.path.isfile(args.input_sample_bam):
        utils.abort(AP, "File bam doesn't exist: %s" % (args.input_sample_bam,))
    if not os.path.isdir(os.path.dirname(os.path.abspath(args.output_hmnfusion_xlsx))):
        utils.abort(AP, "Outdir doesn't exist: %s" % (args.output_hmnfusion_xlsx,))
    if args.output_hmnfusion_json:
        if not os.path.isdir(
            os.path.dirname(os.path.abspath(args.output_hmnfusion_json))
        ):
            utils.abort(AP, "Outdir doesn't exist: %s" % (args.output_hmnfusion_json,))

    if args.size_to_extract % 2 != 0:
        utils.abort(
            AP,
            "Size_to_extract argument should be an even number: %s"
            % (args.size_to_extract,),
        )
    if not utils.validate_name_sample(args.name):
        logging.warning("Name sample is not valid: %s" % (args.name,))

    # Init.
    utils.check_fasta_index(args.input_reference_fasta)
    logging.info("Check index fasta reference")
    if not utils.check_fasta_index(args.input_reference_fasta):
        utils.abort(
            'Index of fasta file doesn"t exist: %s' % (args.input_reference_fasta,)
        )

    # Run.
    logging.info("Load input file")
    g = graph.Graph.from_json(args.input_hmnfusion_json)

    # Subset.
    logging.info("Select fusion")
    fusions = g.subset_graph(
        include=args.fusion_include_flag, exclude=args.fusion_exclude_flag
    )

    # Process fusion.
    logging.info("Get mmej motif from %s fusions" % (len(fusions),))
    dfs = []
    for f in sorted(fusions):
        f.set_mmej(
            path_reference=args.input_reference_fasta,
            path_bam=args.input_sample_bam,
            interval=args.size_to_extract,
        )
        dfs.append(f.mmej_dataframe())

    logging.info("Write fusions to output file")
    mmej_fusion.MmejFusion.to_excel(path=args.output_hmnfusion_xlsx, dfs=dfs)
    if args.output_hmnfusion_json:
        g.to_json(
            path=args.output_hmnfusion_json,
            metadata=dict(name=args.name, step="mmej-fusion"),
        )
    logging.info("End analysis - MMEJ Fusion")


P_mmej_fusion = AP_subparsers.add_parser("mmej-fusion", help=_cmd_mmej_fusion.__doc__)
P_mmej_fusion.add_argument(
    "--input-hmnfusion-json", required=True, help="HmnFusion, json file"
)
P_mmej_fusion.add_argument("--input-sample-bam", required=True, help="Bam file")
P_mmej_fusion.add_argument(
    "--input-reference-fasta", required=True, help="Reference, fasta file"
)
P_mmej_fusion.add_argument(
    "--fusion-include-flag",
    default=8,
    type=int,
    help="Select fusions with fusion-flag",
)
P_mmej_fusion.add_argument(
    "--fusion-exclude-flag",
    default=7,
    type=int,
    help="Exclude fusions with fusion-flag",
)
P_mmej_fusion.add_argument(
    "--size-to-extract",
    type=int,
    default=60,
    help="Size of sequence to extract before and after the genomic coordinate (even number)",
)
P_mmej_fusion.add_argument("--name", required=True, help="Name of sample")
P_mmej_fusion.add_argument(
    "--output-hmnfusion-xlsx", required=True, help="Excel file output"
)
P_mmej_fusion.add_argument("--output-hmnfusion-json", help="Json file output")
P_mmej_fusion.set_defaults(func=_cmd_mmej_fusion)


def _cmd_wkf_align(args):
    """Worflow to build BAM files from FASTQ files"""
    logging.info("Start - Worfklow Align")
    # Args.
    if not os.path.isfile(args.input_forward_fastq):
        utils.abort(
            AP, "File Fastq Forward doesn't exist : %s" % (args.input_forward_fastq,)
        )
    if not os.path.isfile(args.input_reverse_fastq):
        utils.abort(
            AP, "File Fastq Reverse doesn't exist : %s" % (args.input_reverse_fastq,)
        )
    input_config_json = args.input_config_json
    if input_config_json is None or not os.path.isfile(input_config_json):
        logging.warning("Use default config base")
        input_config_json = os.path.join(
            os.path.dirname(os.path.realpath(__file__)),
            "templates",
            "snakefile",
            "config.align.json",
        )
    input_design_bed = args.input_design_bed
    if input_design_bed is None or not os.path.isfile(input_design_bed):
        logging.warning("Use default design bed")
        input_design_bed = os.path.join(
            os.path.dirname(os.path.realpath(__file__)),
            "templates",
            "bed",
            "design.bed",
        )

    # Run.
    logging.info("Load config file")
    config = utils.read_json(path=input_config_json)
    config_run = {
        "run": {
            "input_forward_fastq": args.input_forward_fastq,
            "input_reverse_fastq": args.input_reverse_fastq,
            "name": args.name,
            "input_design_bed": input_design_bed,
            "output_directory": args.output_directory,
            "tmpdir": args.tmpdir,
            "platform": args.platform,
            "threads": args.threads,
            "cores": args.threads * 2,
        }
    }
    config.update(config_run)

    logging.info("Run Worflow Align")
    workflow.run(
        snakefile=os.path.join(
            os.path.dirname(os.path.realpath(__file__)),
            "templates",
            "snakefile",
            "Snakefile.align",
        ),
        config=config,
        cores=args.threads,
    )
    logging.info("End - Worfklow Align")


P_wkf_align = AP_subparsers.add_parser("workflow-align", help=_cmd_wkf_align.__doc__)
P_wkf_align.add_argument("--input-config-json", help="Input config file")
P_wkf_align.add_argument(
    "--input-forward-fastq", required=True, help="Fastq file forward"
)
P_wkf_align.add_argument(
    "--input-reverse-fastq", required=True, help="Fastq file reverse"
)
P_wkf_align.add_argument("--name", required=True, help="Name of sample")
P_wkf_align.add_argument("--input-design-bed", help="Design bed file")
P_wkf_align.add_argument(
    "--output-directory", required=True, help="Directory to output"
)
P_wkf_align.add_argument(
    "--tmpdir",
    default=tempfile.gettempdir(),
    help="Directory used for temporary results",
)
P_wkf_align.add_argument(
    "--platform",
    type=str,
    default="ILLUMINA",
    help="Platform label to indicate into RGLINE of the BAM file",
)
P_wkf_align.add_argument(
    "--threads",
    type=int,
    default=1,
    choices=range(1, 7),
    metavar="[1-6]",
    help="Threads used",
)
P_wkf_align.set_defaults(func=_cmd_wkf_align)


def _cmd_wkf_hmnfusion(args):
    """Worflow to run HmnFusion ExtractFusion & Quantification"""
    logging.info("Start - Worfklow HmnFusion")
    # Args.
    genefuse_file = ""
    genefuse_fmt = ""
    if args.input_genefuse_json is None:
        genefuse_file = args.input_genefuse_html
        genefuse_fmt = "html"
    else:
        genefuse_file = args.input_genefuse_json
        genefuse_fmt = "json"
    if not os.path.isfile(genefuse_file):
        utils.abort(AP, "File Genefuse doesn't exist : %s" % (genefuse_file,))
    if not os.path.isfile(args.input_lumpy_vcf):
        utils.abort(AP, "File Lumpy doesn't exist : %s" % (args.input_lumpy_vcf,))
    if not os.path.isfile(args.input_sample_bam):
        utils.abort(AP, "File bam doesn't exist : %s" % (args.input_sample_bam,))
    input_hmnfusion_bed = args.input_hmnfusion_bed
    if input_hmnfusion_bed is None or not os.path.isfile(input_hmnfusion_bed):
        logging.warning("Use default hmnfusion bed")
        input_hmnfusion_bed = os.path.join(
            os.path.dirname(os.path.realpath(__file__)),
            "templates",
            "bed",
            "hmnfusion.bed",
        )
    if not os.path.isdir(os.path.dirname(os.path.abspath(args.output_hmnfusion_vcf))):
        utils.abort(AP, "Outdir doesn't exist : %s" % (args.output_hmnfusion_vcf,))

    if not utils.validate_name_sample(args.name):
        logging.warning("Name sample is not valid: %s" % (args.name,))

    # Run.
    logging.info("Run Worflow HmnFusion")
    config = {
        "input_sample_bam": args.input_sample_bam,
        "name": args.name,
        "input_hmnfusion_bed": input_hmnfusion_bed,
        "output_hmnfusion_vcf": args.output_hmnfusion_vcf,
        "lumpy": args.input_lumpy_vcf,
        "genefuse": genefuse_file,
        "genefuse_fmt": genefuse_fmt,
    }
    workflow.run(
        snakefile=os.path.join(
            os.path.dirname(os.path.realpath(__file__)),
            "templates",
            "snakefile",
            "Snakefile.hmnfusion",
        ),
        config=config,
        cores=1,
    )
    logging.info("End - Worfklow HmnFusion")


P_wkf_hmnfusion = AP_subparsers.add_parser(
    "workflow-hmnfusion", help=_cmd_wkf_hmnfusion.__doc__
)
P_wkf_hmnfusion.add_argument("--input-lumpy-vcf", required=True, help="Lumpy Vcf file")
P_wkf_hf_g = P_wkf_hmnfusion.add_mutually_exclusive_group(required=True)
P_wkf_hf_g.add_argument("--input-genefuse-json", help="Genefuse, json file")
P_wkf_hf_g.add_argument("--input-genefuse-html", help="Genefuse, html file")
P_wkf_hmnfusion.add_argument("--input-sample-bam", required=True, help="Bam file")
P_wkf_hmnfusion.add_argument("--input-hmnfusion-bed", help="HmnFusion bed file")
P_wkf_hmnfusion.add_argument("--name", required=True, help="Name of sample")
P_wkf_hmnfusion.add_argument(
    "--output-hmnfusion-vcf", required=True, help="Vcf file output"
)
P_wkf_hmnfusion.set_defaults(func=_cmd_wkf_hmnfusion)


def _cmd_wkf_fusion(args):
    """Worflow to detect & quantify fusions with Genefuse, Lumpy and HmnFusion"""
    logging.info("Start - Worfklow Fusion")
    # Args.
    if not os.path.isfile(args.input_sample_bam):
        utils.abort(AP, "File bam doesn't exist : %s" % (args.input_sample_bam,))
    if not utils.check_bam_index(args.input_sample_bam):
        utils.abort(
            AP, "Input alignment file must be in BAM format, index could not be build"
        )

    input_forward_fastq = args.input_forward_fastq
    input_reverse_fastq = args.input_reverse_fastq
    is_convert_bam = False
    if (
        input_forward_fastq is None
        or not os.path.isfile(input_forward_fastq)
        or input_reverse_fastq is None
        or not os.path.isfile(input_reverse_fastq)
    ):
        is_convert_bam = True

    input_hmnfusion_bed = args.input_hmnfusion_bed
    if input_hmnfusion_bed is None or not os.path.isfile(input_hmnfusion_bed):
        logging.warning("Use default hmnfusion bed")
        input_hmnfusion_bed = os.path.join(
            os.path.dirname(os.path.realpath(__file__)),
            "templates",
            "bed",
            "hmnfusion.bed",
        )
    input_genefuse_bed = args.input_genefuse_bed
    if input_genefuse_bed is None or not os.path.isfile(input_genefuse_bed):
        logging.warning("Use default genefuse bed")
        input_genefuse_bed = os.path.join(
            os.path.dirname(os.path.realpath(__file__)),
            "templates",
            "bed",
            "genefuse.bed",
        )
    input_lumpy_bed = args.input_lumpy_bed
    if input_lumpy_bed is None or not os.path.isfile(input_lumpy_bed):
        logging.warning("Use default lumpy bed")
        input_lumpy_bed = os.path.join(
            os.path.dirname(os.path.realpath(__file__)), "templates", "bed", "lumpy.bed"
        )
    if not os.path.isfile(args.input_reference_fasta):
        utils.abort(
            AP, "File reference doesn't exist : %s" % (args.input_reference_fasta,)
        )
    if not os.path.isdir(os.path.dirname(os.path.abspath(args.output_hmnfusion_vcf))):
        utils.abort(AP, "Outdir doesn't exist : %s" % (args.output_hmnfusion_vcf,))
    if not os.path.isdir(os.path.dirname(os.path.abspath(args.output_genefuse_html))):
        utils.abort(AP, "Outdir doesn't exist : %s" % (args.output_genefuse_html,))
    if not os.path.isdir(os.path.dirname(os.path.abspath(args.output_lumpy_vcf))):
        utils.abort(AP, "Outdir doesn't exist : %s" % (args.output_lumpy_vcf,))

    if not utils.validate_name_sample(args.name):
        logging.warning("Name sample is not valid: %s" % (args.name,))

    # Threads
    threads_genefuse = 1
    if args.threads > 1:
        threads_genefuse = math.ceil(args.threads * 0.8)

    # Run.
    if is_convert_bam:
        logging.info("Convert bam to fastq")
        input_forward_fastq, input_reverse_fastq = utils.bam_to_fastq(
            path=args.input_sample_bam, threads=args.threads
        )

    logging.info("Run Workflow Fusion")
    config = {
        "input_forward_fastq": input_forward_fastq,
        "input_reverse_fastq": input_reverse_fastq,
        "input_sample_bam": args.input_sample_bam,
        "name": args.name,
        "input_genefuse_bed": input_genefuse_bed,
        "input_hmnfusion_bed": input_hmnfusion_bed,
        "input_lumpy_bed": input_lumpy_bed,
        "input_reference_fasta": args.input_reference_fasta,
        "output_hmnfusion_vcf": args.output_hmnfusion_vcf,
        "genefuse": args.output_genefuse_html,
        "lumpy": args.output_lumpy_vcf,
        "threads_genefuse": threads_genefuse,
    }
    workflow.run(
        snakefile=os.path.join(
            os.path.dirname(os.path.realpath(__file__)),
            "templates",
            "snakefile",
            "Snakefile.fusion",
        ),
        config=config,
        cores=args.threads,
    )
    logging.info("End - Worfklow Fusion")


P_wkf_fusion = AP_subparsers.add_parser("workflow-fusion", help=_cmd_wkf_fusion.__doc__)
P_wkf_fusion.add_argument("--input-forward-fastq", help="Fastq file forward")
P_wkf_fusion.add_argument("--input-reverse-fastq", help="Fastq file reverse")
P_wkf_fusion.add_argument("--input-sample-bam", required=True, help="Bam file")
P_wkf_fusion.add_argument("--input-genefuse-bed", help="Genefuse bed file")
P_wkf_fusion.add_argument("--input-lumpy-bed", help="Lumpy bed file")
P_wkf_fusion.add_argument("--input-hmnfusion-bed", help="HmnFusion bed file")
P_wkf_fusion.add_argument(
    "--input-reference-fasta", required=True, help="Reference fasta file (hg19)"
)
P_wkf_fusion.add_argument("--name", required=True, help="Name of sample")
P_wkf_fusion.add_argument(
    "--output-hmnfusion-vcf", required=True, help="Vcf file output"
)
P_wkf_fusion.add_argument(
    "--output-genefuse-html", required=True, help="Genefuse html file output"
)
P_wkf_fusion.add_argument(
    "--output-lumpy-vcf", required=True, help="Lumpy vcf file output"
)
P_wkf_fusion.add_argument(
    "--threads",
    type=int,
    default=1,
    choices=range(1, 7),
    metavar="[1-6]",
    help="Threads used",
)
P_wkf_fusion.set_defaults(func=_cmd_wkf_fusion)


# Installation software.
def _cmd_install_software(args):
    """Install all softwares"""
    logging.info("Start - install-software")

    logging.info("Check if software required are installed")
    if not install_software.InstallSoftware.required():
        utils.abort("Failed")
    if args.uninstall:
        if not install_software.InstallSoftware.uninstall_genefuse():
            logging.warning("GeneFuse could not be uninstall")
        if not install_software.InstallSoftware.uninstall_lumpy():
            logging.warning("Lumpy could not be uninstall")
    else:
        path = install_software.InstallSoftware.get_path_folder()
        if path is None:
            utils.abort(
                "No folder in the path could not be written, change user to more level privileges"
            )
        if not install_software.InstallSoftware.install_genefuse(path=path):
            utils.abort("Failed")
        if not install_software.InstallSoftware.install_lumpy():
            utils.abort("Failed")

    logging.info("End - install-software")


P_install_software = AP_subparsers.add_parser(
    "install-software", help=_cmd_install_software.__doc__
)
P_install_software.add_argument(
    "--uninstall",
    action="store_true",
    help="Remove GeneFuse & lumpy-sv conda environment",
)
P_install_software.set_defaults(func=_cmd_install_software)


# fusion flag.
def _cmd_fusion_flag(args):
    """Show all fusion flag"""
    # Reset logging
    for handler in logging.root.handlers:
        logging.root.removeHandler(handler)
    logging.basicConfig(
        level=logging.INFO,
        format="%(message)s",
    )
    # Grep all attributes
    values = []
    for attr in dir(fusion_flag.FusionFlag):
        if not attr.startswith("__"):
            values.append((attr, fusion_flag.FusionFlag[attr].value))
    # Print values
    for value in sorted(values, key=lambda x: x[1]):
        logging.info("%s\t%s" % (value[0], value[1]))


P_fusion_flag = AP_subparsers.add_parser("fusion-flag", help=_cmd_fusion_flag.__doc__)
P_fusion_flag.set_defaults(func=_cmd_fusion_flag)


def _cmd_download_zenodo(args):
    """Download reference files from Zenodo"""
    logging.info("Start - download-zenodo")
    # Parse args.
    output_directory = os.path.abspath(args.output_directory)
    if not os.path.isdir(output_directory):
        logging.warning("Create output-directory: %s" % (output_directory,))
        os.makedirs(output_directory, exist_ok=True)
    zenodo_id = args.input_zenodo_str
    access_token = None
    if args.input_token_str:
        access_token = args.input_token_str
    elif args.input_token_txt:
        if not os.path.isfile(args.input_token_txt):
            utils.abort(
                AP, "File TokenAccess doesn't exist : %s" % (args.input_token_txt,)
            )
        access_token = download_reference.DownloadZenodo.load_token(
            path=args.input_token_txt
        )
    params = dict()
    if access_token:
        params["access_token"] = access_token
    # Init
    dwl = download_reference.DownloadZenodo(id=zenodo_id, params=params)
    # Grep data.
    items = dwl.show_files()
    for item in items:
        logging.info("Process: %s" % (item["filename"],))
        foutput = os.path.join(output_directory, item["filename"])
        dwl.download_file(url=item["links"]["download"], path=foutput)
        download_reference.DownloadZenodo.check_checksum(
            checksum=item["checksum"], path=foutput
        )

    logging.info("End - download-zenodo")


P_download_zenodo = AP_subparsers.add_parser(
    "download-zenodo", help=_cmd_download_zenodo.__doc__
)
P_download_zenodo.add_argument(
    "--input-zenodo-str", required=True, help="Id of Zenodo repository"
)
P_token = P_download_zenodo.add_mutually_exclusive_group(required=False)
P_token.add_argument(
    "--input-token-str",
    help="Token string to access into repository with restricted access",
)
P_token.add_argument(
    "--input-token-txt",
    help="File with token to access into repository with restricted access",
)
P_download_zenodo.add_argument(
    "--output-directory", required=True, help="Directory to write files"
)
P_download_zenodo.set_defaults(func=_cmd_download_zenodo)


# Version.
def print_version(_args):
    """Display this program"s version"""
    print(_version.__version__)


P_version = AP_subparsers.add_parser("version", help=print_version.__doc__)
P_version.set_defaults(func=print_version)


# Help.
def print_help():
    """Display this program"s help"""
    print(AP_subparsers.help)
    AP.exit()


# Main.
def parse_args(args=None):
    """Parse the command line"""
    return AP.parse_args(args=args)
