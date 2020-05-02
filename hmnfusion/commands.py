import argparse
import logging
import os
import re
import sys

import numpy as np
import pandas as pd

from collections import namedtuple
from multiprocessing import Pool
from .utils import abort
from . import (extractfusion, quantification)
from ._version import __version__

from concurrent import futures

AP = argparse.ArgumentParser(
	description="Extraction sequence and quantification of fusions from DNA genomic analysis",
	epilog="See online documentation: ")
AP_subparsers = AP.add_subparsers(
	help="Sub-commnands (use with -h for more info)")

# Extract fusion.
def _cmd_extract_fusion(args):
	'''Extract fusion from Genefuse and Lumpy analysis'''
	logging.info('Start analysis')
	# Grep args.
	finputs = {}
	finputs['genefuse'] = {}
	if args.genefuse_json:
		finputs['genefuse']['path'] = args.genefuse_json
		finputs['genefuse']['format'] = 'json'
	elif args.genefuse_html:
		finputs['genefuse']['path'] = args.genefuse_html
		finputs['genefuse']['format'] = 'html'
	finputs['lumpy'] = {}
	finputs['lumpy']['path'] = args.lumpy_vcf
	finputs['lumpy']['format'] = 'vcf'
	foutput = args.output
	# Check if all exists.
	if not os.path.isfile(finputs['genefuse']['path']):
		abort("File Genefuse doesn't exist : %s"%(fgenefuse['path'],))
	if not os.path.isfile(finputs['lumpy']['path']):
		abort("File Lumpy doesn't exist : %s"%(flumpy['path'],))
	if not os.path.isdir(os.path.dirname(foutput)):
		abort("Outdir doesn't exist : %s"%(foutput,))

	# Run.
	fusions = {}
	# Genefuse.
	logging.info('Genefuse')
	fgenefuse = finputs['genefuse']
	if fgenefuse['format'] == 'json':
		logging.info('\tExtract fusions from Json')
		fusions['genefuse'] = extractfusion.read_genefuse_json(fgenefuse['path'])
	elif fgenefuse['format'] == 'html':
		logging.info('\tExtract fusions from Html')
		fusions['genefuse'] = extractfusion.read_genefuse_html(fgenefuse['path'])
	# Lumpy.	
	logging.info('Lumpy')
	flumpy = finputs['lumpy']
	if flumpy['format'] == 'vcf':
		fusions['lumpy'] = {}
		logging.info('\tExtract fusions')
		fusions['lumpy']['raw'] = extractfusion.read_lumpy(flumpy['path'])
		logging.info('\tBuild consensus')
		fusions['lumpy']['consensus'] = extractfusion.consensus_lumpy(fusions['lumpy']['raw'], args.consensus_interval)
	# Consensus.
	logging.info('Build consensus with interval of %s pb'%(args.consensus_interval,))
	fusions['consensus'] = extractfusion.consensus_genefuse_lumpy(fusions, args.consensus_interval)

	logging.info('Find %s fusion(s)'%(len(fusions['consensus']),))
	for ix, fusion in enumerate(fusions['consensus']):
		logging.info('%s - %s' % (ix,fusion))
	# Write output.
	if foutput:
		logging.info('Write output')
		extractfusion.write(foutput, finputs, fusions)
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
P_extract_fusion.add_argument('-o','--output',
	help='Json file output')
P_extract_fusion.set_defaults(func=_cmd_extract_fusion)

# Quantification of fusions
def _cmd_quantification(args):
	logging.info('Start analysis')

	# Check if all exists.
	foutput = args.output
	if args.hmnfusion_json and not os.path.isfile(args.hmnfusion_json):
		abort("HmnFusion Json file doesn't exist : %s"%(args.hmnfusion_json,))
	if args.region:
		chrom, pos = 0, 0
		m = re.search(r'\D*(\d+):(\d+)', args.region)
		if m:
			chrom = int(m.group(1))
			pos = int(m.group(2))
		if chrom == 0 or pos == 0:
			abort("Region format is not well formated. Required <chrom>:<position>")
	if not os.path.isfile(args.input_bam):
		abort("Bam file doesn't exist : %s"%(args.input_bam,))
	if not os.path.isfile(args.input_bed):
		abort("Bed file doesn't exist : %s"%(args.input_bed,))
	if not os.path.isdir(os.path.dirname(foutput)):
		abort("Outdir doesn't exist : %s"%(foutput,))

	# Parsing bed file.

	fusions = []	
	elif args.hmnfusion_json:
		fusions = quantification.parse_hmnfusion_json(args.hmnfusion_json)

P_quantification = AP_subparsers.add_parser('quantification', help=_cmd_quantification.__doc__)
P_quantification_position_group = P_quantification.add_mutually_exclusive_group(required=True)
P_quantification_position_group.add_argument('--region', 
	help='Region format <chrom>:<postion>')
P_extract_fusion_genefuse_group.add_argument('--hmnfusion-json',
	help='Output Json produced by command "extractfusion"')
P_quantification.add_argument('--input-bam', 
	help='Bam file')
P_quantification.add_argument('--input-bed', 
	help='Bed file')
P_quantification.add_argument('-o','--output',
	help='Json file output')
P_quantification.set_defaults(func=_cmd_quantification)

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
