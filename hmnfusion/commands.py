import argparse
import logging
import os
import sys

from .utils import abort
from . import (extractfusion, quantification, mmej, graph)
from ._version import __version__

AP = argparse.ArgumentParser(
    description="Extraction sequence and quantification of fusions from DNA genomic analysis",
    epilog="See online documentation: ")
AP_subparsers = AP.add_subparsers(
    help="Sub-commnands (use with -h for more info)")

# Extract fusion.
def _cmd_extract_fusion(args):
    """Extract fusion from Genefuse and Lumpy analysis"""
    logging.info('Start analysis')
    # Grep args.
    finputs = {}
    finputs['genefuse'] = {}
    if args.genefuse_json:
        finputs['genefuse']['path'] = os.path.abspath(args.genefuse_json)
        finputs['genefuse']['format'] = 'json'
    elif args.genefuse_html:
        finputs['genefuse']['path'] = os.path.abspath(args.genefuse_html)
        finputs['genefuse']['format'] = 'html'
    finputs['lumpy'] = {}
    finputs['lumpy']['path'] = os.path.abspath(args.lumpy_vcf)
    finputs['lumpy']['format'] = 'vcf'
    foutput = args.output_json
    # Check if all exists.
    if not os.path.isfile(finputs['genefuse']['path']):
        abort(AP, "File Genefuse doesn't exist : %s"%(finputs['genefuse']['path'],))
    if not os.path.isfile(finputs['lumpy']['path']):
        abort(AP, "File Lumpy doesn't exist : %s"%(finputs['lumpy']['path'],))
    if not os.path.isdir(os.path.dirname(os.path.abspath(foutput))):
        abort(AP, "Outdir doesn't exist : %s"%(foutput,))

    # Run.
    g = graph.Graph(args.consensus_interval)
    # Genefuse.
    logging.info('Parse Genefuse')
    extractfusion.read_genefuse(g, finputs['genefuse']['path'], finputs['genefuse']['format'])

    # Lumpy.    
    logging.info('Parse Lumpy')
    extractfusion.read_lumpy(g, finputs['lumpy']['path'], finputs['lumpy']['format'])

    # Consensus.
    logging.info('Build consensus with interval of %s pb'%(g.graph.graph['consensus_interval'],))
    g.consensus_single()
    g.consensus_genefuse_lumpy()
   
    # Define nodes of interest.
    logging.info('Define nodes of interest')
    g.define_node_interest()

    # Filter same chrom.
    logging.info('Filter nodes')
    g.trim_node()

    # Write output.
    if foutput:
        logging.info('Write output')
        extractfusion.write_hmnfusion_json(foutput, finputs, g)
    logging.info("Analysis is finished")

P_extract_fusion = AP_subparsers.add_parser('extractfusion', help=_cmd_extract_fusion.__doc__)
P_extract_fusion_genefuse_group = P_extract_fusion.add_mutually_exclusive_group(required=True)
P_extract_fusion_genefuse_group.add_argument('--genefuse-json', 
    help='Genefuse, json file')
P_extract_fusion_genefuse_group.add_argument('--genefuse-html',
    help='Genefuse, html file')
P_extract_fusion.add_argument('--lumpy-vcf', required=True,
    help='Lumpy vcf file')
P_extract_fusion.add_argument('--consensus-interval', type=int, default=500,
    help='Interval, pb, for which Fusion are considered equal if their chrom are')
P_extract_fusion.add_argument('--output-json',
    help='Json file output')
P_extract_fusion.set_defaults(func=_cmd_extract_fusion)

# Quantification of fusions
def _cmd_quantification(args):
    """Quantify fusion candidates with a BAM file"""
    logging.info('Start analysis')

    # Check if all exists.
    logging.info('Check args')
    finputs = {}
    foutput = args.output_vcf
    finputs['output'] = foutput
    if args.hmnfusion_json:
        finputs['hmnfusion_json'] = args.hmnfusion_json
        if not os.path.isfile(finputs['hmnfusion_json']):
            abort(AP, "HmnFusion Json file doesn't exist : %s"%(finputs['hmnfusion_json'],))
    if args.region:
        if not quantification.check_region(args.region):
            abort(AP, "Region format is not well formated. Required <chrom>:<position>")
        finputs['region'] = args.region

    falignment_path = ''
    falignment_mode = ''
    if args.input_bam:
        falignment_path = args.input_bam
        falignment_mode = 'rb'
    if args.input_sam:
        falignment_path = args.input_sam
        falignment_mode = 'r'
    finputs['alignment'] = {}
    finputs['alignment']['path'] = falignment_path
    finputs['alignment']['mode'] = falignment_mode
    if not os.path.isfile(finputs['alignment']['path']):
        abort(AP, "Input alignment file doesn't exist : %s"%(finputs['alignment']['path'],))
    finputs['bed'] = args.input_bed
    if not os.path.isfile(args.input_bed):
        abort(AP, "Bed file doesn't exist : %s"%(args.input_bed,))
    if not os.path.isdir(os.path.dirname(os.path.abspath(finputs['output']))):
        abort(AP, "Outdir doesn't exist : %s"%(finputs['output'],))

    params = dict(falignment=dict(path=finputs['alignment']['path'], mode=finputs['alignment']['mode']), clipped=dict(count=args.baseclipped_count, interval=args.baseclipped_interval))

    if params['clipped']['count'] < 1: 
        logging.warning('Parameter base clipped count is too low, is set to 1 pb')
        params['clipped']['count'] = 1
    elif params['clipped']['count'] > 20:
        logging.warning('Parameter base clipped count is too high, is set to 20 pb')
        params['clipped']['count'] = 20
    if params['clipped']['interval'] < 1: 
        logging.warning('Parameter base clipped interval is too low, is set to 1 pb')
        params['clipped']['interval'] = 1
    elif params['clipped']['interval'] > 20:
        logging.warning('Parameter base clipped interval is too high, is set to 20 pb')
        params['clipped']['interval'] = 20
    if params['clipped']['count'] > params['clipped']['interval']:
        logging.warning('Parameter base clipped count is higher than parameter base clipped interval -> count is set equal than interval : %s pb instead of %s pb'%(params['clipped']['interval'], params['clipped']['count']))
        params['clipped']['count'] = params['clipped']['interval']
    
    # Parsing bed file.
    logging.info('Parsing bed file')
    bed = quantification.read_bed(args.input_bed)
    
    # Parsing fusions.
    logging.info('Get region')
    g = graph.Graph()
    if args.region:
        fusion = Fusion()
        region = Region(finputs['region'].split(':')[0], finputs['region'].split(':')[1].split('-')[0])
        fusion.set_region(region)
        fusion.evidence.raw = 0
        g.add_node(fusion, 0, False, True)
    elif args.hmnfusion_json:
        g = extractfusion.read_hmnfusion_json(finputs['hmnfusion_json'])

    # Process
    logging.info('Calcul VAF fusion')
    fusions = quantification.run(params, bed, g)

    for ix, fusion in enumerate(fusions):
        logging.info('%s - %s' % (ix,fusion))
    # Write output.
    if foutput:
        logging.info('Write output')
        quantification.write(foutput, args.name, g)
    logging.info("Analysis is finished")

P_quantification = AP_subparsers.add_parser('quantification', help=_cmd_quantification.__doc__)
P_quantification_position_group = P_quantification.add_mutually_exclusive_group(required=True)
P_quantification_position_group.add_argument('--region', 
    help='Region format <chrom>:<postion>')
P_quantification_position_group.add_argument('--hmnfusion-json',
    help='Output Json produced by command "extractfusion"')
P_quantification_input_alignment_group = P_quantification.add_mutually_exclusive_group(required=True)
P_quantification_input_alignment_group.add_argument('--input-bam', 
    help='Bam file')
P_quantification_input_alignment_group.add_argument('--input-sam', 
    help='Sam file')
P_quantification.add_argument('--input-bed', 
    help='Bed file')
P_quantification.add_argument('--name', required=True,
    help='Name of sample')
P_quantification.add_argument('--baseclipped-interval', type=int, default=6, 
    help='Interval to count hard/soft-clipped bases from fusion point (pb)')
P_quantification.add_argument('--baseclipped-count', type=int, default=4, 
    help='Number of base hard/soft-clipped bases to count in interval (pb)')
P_quantification.add_argument('--output-vcf',
    help='Vcf file output')
P_quantification.set_defaults(func=_cmd_quantification)

# mmej signatures.
def _cmd_mmej(args):
    """Detect MMEJ signatures from deletion, fusion"""
    logging.info('Start analysis')
    # Grep args.
    for finput in args.input_vcf:
        if not os.path.isfile(finput):
            abort(AP, "Vcf file doesn't exist : %s"%(finput,))
    if not os.path.isfile(args.reference):
        abort(AP, "Reference file doesn't exist : %s"%(args.reference,))
    if not os.path.isdir(os.path.dirname(os.path.abspath(args.output_xlsx))):
        abort(AP, "Outdir doesn't exist : %s"%(args.output_xlsx,))

    # Run.
    logging.info('Extract events from files')
    df = mmej.extract(args.input_vcf, event=args.event)

    logging.info('Caracterize events with reference file')
    df = mmej.signatures(args.reference, df)

    logging.info('MMEJ signatures significance')
    df = mmej.conclude(df)

    logging.info('Write output')
    mmej.write(args.output_xlsx, df)

    logging.info("Analysis is finished")

P_mmej = AP_subparsers.add_parser('mmej', help=_cmd_mmej.__doc__)
P_mmej.add_argument('-i', '--input-vcf', nargs='+',
    help='Vcf file')
P_mmej.add_argument('-r', '--reference', required=True,
    help='Genome of reference')
P_mmej.add_argument('-e', '--event', default='deletion', choices=['deletion', 'fusion', 'all'], 
    help='Event from which to detect MMEJ')
P_mmej.add_argument('--output-xlsx',
    help='Output file')
P_mmej.set_defaults(func=_cmd_mmej)


# Version.
def print_version(_args):
    """Display this program's version"""
    print(__version__)

P_version = AP_subparsers.add_parser('version', help=print_version.__doc__)
P_version.set_defaults(func=print_version)

# Help.
def print_help():
    """Display this program's help"""
    print(AP_subparsers.help)
    AP.exit()

# Main.
def parse_args(args=None):
    """Parse the command line"""
    return AP.parse_args(args=args)
