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
	'''Extract fusion from Genfuse and Lumpy analysis'''
	# Checking args.
	fgenefuse = dict(path='',form='')
	if args.genefuse_json:
		fgenefuse['path'] = args.genefuse_json
		fgenefuse['form'] = 'json'
	else if args.genefuse_html:
		fgenefuse['path'] = args.genefuse_html
		fgenefuse['form'] = 'html'
	flumpy = dict(path=args.lumpy_vcf, form='vcf')
	foutput = args.output
	
	if ! os.path.isfile(fgenefuse['path']):
		abort("File Genefuse doesn't exist : %s"%(fgenefuse['path'],))
	if ! os.path.isfile(flumpy['path']):
		abort("File Lumpy doesn't exist : %s"%(flumpy['path'],))
	if ! os.path.isdir(os.path.dirname(foutput)):
		abort("Outdir doesn't exist : %s"%(foutput,))
		
	# Parsing genefuse.
	extractfusion.parse_genefuse(fgenefuse)

	# Parsing lumpy.
	extractfusion.parse_lumpy(flumpy)
	
	# Build.
	
	# Write output.
	
	
	logging.info("Analysis is finished")

P_extract_fusion = AP_subparsers.add_parser('extractfusion', help=_cmd_extract_fusion.__doc__)
P_extract_fusion_genefuse_group = parser.add_mutually_exclusive_group(required=True)
P_extract_fusion_genefuse_group.add_argument('--genefuse-json', 
	help='Genefuse, json file')
P_extract_fusion_genefuse_group.add_argument('--genefuse-html',
	help='Genefuse, html file')
P_extract_fusion.add_argument('--lumpy-vcf', required=True,
	help='Lumpy vcf file')
P_extract_fusion.add_argument('-o','--output', required=True, 
	help='Txt output')

P_extract_fusion.set_defaults(func=_cmd_extract_fusion)

# Quantification of fusions
def _cmd_quantification(args):
	logging.info("Start analysis")
	pass

P_quantification = AP_subparsers.add_parser('quantification', help=_cmd_quantification.__doc__)
P_quantification.add_argument('-i', '--input', nargs='+',
	help='Vcf files')
P_quantification.add_argument('-m', '--mode', choices=extractvcf.MODES, default=extractvcf.MODES[0],
	help='Which information to extract')
  
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
