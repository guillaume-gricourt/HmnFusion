import argparse
import logging
import math
import os

from hmnfusion import _version
from hmnfusion import bed as ibed
from hmnfusion import (
    extractfusion,
    fusion,
    graph,
    mmej,
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
        extractfusion.write_hmnfusion_json(foutput, finputs, g)
    logging.info("Analysis is finished")


P_extract_fusion = AP_subparsers.add_parser(
    "extractfusion", help=_cmd_extract_fusion.__doc__
)
P_efgg = P_extract_fusion.add_mutually_exclusive_group(required=True)
P_efgg.add_argument("--input-genefuse-json", help="Genefuse, json file")
P_efgg.add_argument("--input-genefuse-html", help="Genefuse, html file")
P_extract_fusion.add_argument("--input-lumpy-vcf", required=True, help="Lumpy vcf file")
P_extract_fusion.add_argument("--input-hmnfusion-bed", required=True, help="Bed file")
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
    finputs = {}
    foutput = args.output_hmnfusion_vcf
    finputs["output"] = foutput
    if args.input_hmnfusion_json:
        finputs["hmnfusion_json"] = args.input_hmnfusion_json
        if not os.path.isfile(finputs["hmnfusion_json"]):
            utils.abort(
                AP,
                "HmnFusion Json file doesn't exist : %s" % (finputs["hmnfusion_json"],),
            )
    if args.region:
        if not region.Region.check_region(args.region):
            utils.abort(
                AP,
                "Region format is not well formated. \
                Required <chrom>:<position>",
            )
        finputs["region"] = args.region

    finputs["alignment"] = {}
    finputs["alignment"]["path"] = args.input_sample_bam
    finputs["alignment"]["mode"] = "rb"
    if not os.path.isfile(finputs["alignment"]["path"]):
        utils.abort(
            AP,
            "Input alignment file doesn't exist : %s" % (finputs["alignment"]["path"],),
        )
    if not utils.check_bam_index(finputs["alignment"]["path"]):
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
    finputs["bed"] = bed_hmnfusion
    if not os.path.isdir(os.path.dirname(os.path.abspath(finputs["output"]))):
        utils.abort(AP, "Outdir doesn't exist : %s" % (finputs["output"],))

    params = dict(
        falignment=dict(
            path=finputs["alignment"]["path"], mode=finputs["alignment"]["mode"]
        ),
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

    # Parsing bed file.
    logging.info("Parsing bed file")
    bed = ibed.Bed.from_bed(args.input_hmnfusion_bed)

    # Parsing fusions.
    logging.info("Get region")
    g = graph.Graph()
    if args.region:
        fus = fusion.Fusion()
        reg = region.Region(
            finputs["region"].split(":")[0],
            finputs["region"].split(":")[1].split("-")[0],
        )
        fus.set_region(reg)
        fus.evidence.raw = 0
        g.add_node(fus, 0, False, True)
    elif args.input_hmnfusion_json:
        g = extractfusion.read_hmnfusion_json(finputs["hmnfusion_json"])

    # Process
    logging.info("Calcul VAF fusion")
    quantification.run(params, bed, g)

    # Write output.
    if foutput:
        logging.info("Write output")
        quantification.write(foutput, args.name, g)
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
P_quantification.set_defaults(func=_cmd_quantification)


# mmej signatures.
def _cmd_mmej_deletion(args):
    """Detect MMEJ signatures from deletion, fusion"""
    logging.info("Start analysis")
    # Grep args.
    for finput in args.input_sample_vcf:
        if not os.path.isfile(finput):
            utils.abort(AP, 'Vcf file doesn"t exist : %s' % (finput,))
    if not os.path.isfile(args.input_reference_fasta):
        utils.abort(
            AP, 'Reference file doesn"t exist : %s' % (args.input_reference_fasta,)
        )
    if not os.path.isdir(os.path.dirname(os.path.abspath(args.output_hmnfusion_xlsx))):
        utils.abort(AP, 'Outdir doesn"t exist : %s' % (args.output_hmnfusion_xlsx,))

    # Run.
    logging.info("Extract events from files")
    df = mmej.extract(args.input_sample_vcf)

    logging.info("Caracterize events with reference file")
    df = mmej.signatures(args.input_reference_fasta, df)

    logging.info("MMEJ signatures significance")
    df = mmej.conclude(df)

    logging.info("Write output")
    mmej.write(args.output_hmnfusion_xlsx, df)

    logging.info("Analysis is finished")


P_mmej = AP_subparsers.add_parser("mmej-deletion", help=_cmd_mmej_deletion.__doc__)
P_mmej.add_argument("--input-sample-vcf", nargs="+", help="Vcf file")
P_mmej.add_argument(
    "--input-reference-fasta", required=True, help="Genome of reference"
)
P_mmej.add_argument("--output-hmnfusion-xlsx", help="Output file")
P_mmej.set_defaults(func=_cmd_mmej_deletion)


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
    if not os.path.isfile(args.input_forward_fastq):
        utils.abort(
            AP, "File Fastq Forward doesn't exist : %s" % (args.input_forward_fastq,)
        )
    if not os.path.isfile(args.input_reverse_fastq):
        utils.abort(
            AP, "File Fastq Reverse doesn't exist : %s" % (args.input_reverse_fastq,)
        )
    if not os.path.isfile(args.input_sample_bam):
        utils.abort(AP, "File bam doesn't exist : %s" % (args.input_sample_bam,))
    if not utils.check_bam_index(args.input_sample_bam):
        utils.abort(
            AP, "Input alignment file must be in BAM format, index could not be build"
        )
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

    # Threads
    threads_genefuse = 1
    if args.threads > 1:
        threads_genefuse = math.ceil(args.threads * 0.8)

    # Run.
    logging.info("Run Workflow Fusion")
    config = {
        "input_forward_fastq": args.input_forward_fastq,
        "input_reverse_fastq": args.input_reverse_fastq,
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
P_wkf_fusion.add_argument(
    "--input-forward-fastq", required=True, help="Fastq file forward"
)
P_wkf_fusion.add_argument(
    "--input-reverse-fastq", required=True, help="Fastq file reverse"
)
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


def _cmd_mmej_fusion(args):
    """Identify MMEJ from fusion"""
    logging.info("Start analysis - MMEJ Fusion")
    # Grep args.
    fmt = ""
    input_path = ""
    if args.genefuse_json:
        fmt = "genefuse_json"
        input_path = args.genefuse_json
    elif args.genefuse_html:
        fmt = "genefuse_html"
        input_path = args.genefuse_html
    elif args.lumpy_vcf:
        fmt = "lumpy_vcf"
        input_path = args.lumpy_vcf
    elif args.hmnfusion_json:
        fmt = "hmnfusion_json"
        input_path = args.hmnfusion_json

    # Check if all exists.
    if not os.path.isfile(input_path):
        utils.abort(
            AP, "File input doesn't exist : %s" % (input_path,)
        )
    if not os.path.isdir(os.path.dirname(os.path.abspath(args.output_xlsx))):
        utils.abort(AP, "Outdir doesn't exist : %s" % (args.output_xlsx,))

    # Run.
    logging.info("Load input file")
    g = graph.Graph()
    g = mmej_fusion.load_file(input_path, fmt, g)

    # Subset.
    fusions = mmej_fusion.subset_graph(g)

    # Process fusion.
    print('nb of fusion', len(fusions))
    for fusion in fusions[:1]:
        print(fusion)
        path_bam_filter = mmej_fusion.filter_sequence(path=args.input_bam, fus=fusion)
        print(path_bam_filter)
        # path_bam_consensus = mmej_fusion.create_consensus(path_reference=args.input_reference, path_bam=path_bam_filter)
    logging.info("End analysis - MMEJ Fusion")


P_mmej_fus = AP_subparsers.add_parser(
    "mmej-fusion", help=_cmd_mmej_fusion.__doc__
)
P_mmej_fus_input = P_mmej_fus.add_mutually_exclusive_group(required=True)
P_mmej_fus_input.add_argument("--genefuse-json", help="Genefuse, json file")
P_mmej_fus_input.add_argument("--genefuse-html", help="Genefuse, html file")
P_mmej_fus_input.add_argument("--lumpy-vcf", help="Genefuse, html file")
P_mmej_fus_input.add_argument("--hmnfusion-json", help="HmnFusion, json file")
P_mmej_fus.add_argument("--input-bam", help="Bam file")
P_mmej_fus.add_argument("--input-reference", help="Reference, fasta file")
P_mmej_fus.add_argument(
    "--size-to-extract",
    type=int,
    default=30,
    help="Size of sequence to extract before and after the genomic coordinate",
)
P_mmej_fus.add_argument("--output-xlsx", help="Excel file output")
P_mmej_fus.set_defaults(func=_cmd_mmej_fusion)


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
