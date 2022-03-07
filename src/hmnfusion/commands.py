import argparse
import logging
import os

from hmnfusion import _version
from hmnfusion import bed as ibed
from hmnfusion import (
    extractfusion,
    fusion,
    graph,
    mmej,
    quantification,
    region,
    utils,
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
    if args.genefuse_json:
        finputs["genefuse"]["path"] = os.path.abspath(args.genefuse_json)
        finputs["genefuse"]["format"] = "json"
    elif args.genefuse_html:
        finputs["genefuse"]["path"] = os.path.abspath(args.genefuse_html)
        finputs["genefuse"]["format"] = "html"
    finputs["lumpy"] = {}
    finputs["lumpy"]["path"] = os.path.abspath(args.lumpy_vcf)
    finputs["lumpy"]["format"] = "vcf"
    foutput = args.output_json
    # Check if all exists.
    if not os.path.isfile(finputs["genefuse"]["path"]):
        utils.abort(
            AP, "File Genefuse doesn't exist : %s" % (finputs["genefuse"]["path"],)
        )
    if not os.path.isfile(finputs["lumpy"]["path"]):
        utils.abort(AP, "File Lumpy doesn't exist : %s" % (finputs["lumpy"]["path"],))
    if not os.path.isdir(os.path.dirname(os.path.abspath(foutput))):
        utils.abort(AP, "Outdir doesn't exist : %s" % (foutput,))

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
    g.trim_node()

    # Write output.
    if foutput:
        logging.info("Write output")
        extractfusion.write_hmnfusion_json(foutput, finputs, g)
    logging.info("Analysis is finished")


P_extract_fusion = AP_subparsers.add_parser(
    "extractfusion", help=_cmd_extract_fusion.__doc__
)
P_efgg = P_extract_fusion.add_mutually_exclusive_group(required=True)
P_efgg.add_argument("--genefuse-json", help="Genefuse, json file")
P_efgg.add_argument("--genefuse-html", help="Genefuse, html file")
P_extract_fusion.add_argument("--lumpy-vcf", required=True, help="Lumpy vcf file")
P_extract_fusion.add_argument(
    "--consensus-interval",
    type=int,
    default=500,
    help="Interval, pb, for which Fusion are considered equal if their chrom are",
)
P_extract_fusion.add_argument("--output-json", help="Json file output")
P_extract_fusion.set_defaults(func=_cmd_extract_fusion)


# Quantification of fusions
def _cmd_quantification(args):
    """Quantify fusion candidates with a BAM file"""
    logging.info("Start analysis")

    # Check if all exists.
    logging.info("Check args")
    finputs = {}
    foutput = args.output_vcf
    finputs["output"] = foutput
    if args.hmnfusion_json:
        finputs["hmnfusion_json"] = args.hmnfusion_json
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
    finputs["alignment"]["path"] = args.input_bam
    finputs["alignment"]["mode"] = "rb"
    if not os.path.isfile(finputs["alignment"]["path"]):
        utils.abort(
            AP,
            "Input alignment file doesn't exist : %s" % (finputs["alignment"]["path"],),
        )
    finputs["bed"] = args.input_bed
    if not os.path.isfile(args.input_bed):
        utils.abort(AP, "Bed file doesn't exist : %s" % (args.input_bed,))
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
    bed = ibed.Bed.from_bed(args.input_bed)

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
    elif args.hmnfusion_json:
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
    "--hmnfusion-json", help='Output Json produced by command "extractfusion"'
)
P_quantification.add_argument("--input-bam", required=True, help="Bam file")
P_quantification.add_argument("--input-bed", help="Bed file")
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
P_quantification.add_argument("--output-vcf", help="Vcf file output")
P_quantification.set_defaults(func=_cmd_quantification)


# mmej signatures.
def _cmd_mmej(args):
    """Detect MMEJ signatures from deletion, fusion"""
    logging.info("Start analysis")
    # Grep args.
    for finput in args.input_vcf:
        if not os.path.isfile(finput):
            utils.abort(AP, 'Vcf file doesn"t exist : %s' % (finput,))
    if not os.path.isfile(args.reference):
        utils.abort(AP, 'Reference file doesn"t exist : %s' % (args.reference,))
    if not os.path.isdir(os.path.dirname(os.path.abspath(args.output_xlsx))):
        utils.abort(AP, 'Outdir doesn"t exist : %s' % (args.output_xlsx,))

    # Run.
    logging.info("Extract events from files")
    df = mmej.extract(args.input_vcf)

    logging.info("Caracterize events with reference file")
    df = mmej.signatures(args.reference, df)

    logging.info("MMEJ signatures significance")
    df = mmej.conclude(df)

    logging.info("Write output")
    mmej.write(args.output_xlsx, df)

    logging.info("Analysis is finished")


P_mmej = AP_subparsers.add_parser("mmej", help=_cmd_mmej.__doc__)
P_mmej.add_argument("-i", "--input-vcf", nargs="+", help="Vcf file")
P_mmej.add_argument("-r", "--reference", required=True, help="Genome of reference")
P_mmej.add_argument("--output-xlsx", help="Output file")
P_mmej.set_defaults(func=_cmd_mmej)


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
    if not os.path.isdir(os.path.dirname(os.path.abspath(foutput))):
        utils.abort(AP, "Outdir doesn't exist : %s" % (foutput,))

    # Run.
    logging.info("Load input file")
    g = graph.Graph()
    g = mmej_fusion.load_file(input_path, fmt, g)

    # Subset.
    fusions = mmej_fusion.subset_graph(g)

    # Process fusion.
    for fusion in fusions:
        path_bam_filter = mmej_fusion.filter_sequence(path=args.input_bam, fus=fusion)
        path_bam_consensus = mmej_fusion.create_consensus(path_reference=args.input_reference, path_bam=path_bam_filter)
    logging.info("End analysis - MMEJ Fusion")


P_mmej_fus = AP_subparsers.add_parser(
    "mmej-fusion", help=_cmd_mmej_fusion.__doc__
)
P_mmej_fus_genefuse = P_mmej_fus.add_mutually_exclusive_group(required=True)
P_mmej_fus_genefuse.add_argument("--genefuse-json", help="Genefuse, json file")
P_mmej_fus_genefuse.add_argument("--genefuse-html", help="Genefuse, html file")
P_mmej_fus.add_argument("--lumpy-vcf", help="Genefuse, html file")
P_mmej_fus.add_argument("--hmnfusion-json", help="HmnFusion, json file")
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
