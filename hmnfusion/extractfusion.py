import bs4
import json
import logging
import os
import pysam
import re
import sys

from pysam import VariantFile

from .utils import read_json, update_list, write_json, Region, Fusion

RE_GENEFUSE_LABEL = re.compile(r'[+-](\w+):(\d+)')
RE_GENEFUSE_TOTAL = re.compile(r'total:\s{1}(\d+),')
RE_LUMPY_CHROM = re.compile(r'\D*(\d+)')
RE_LUMPY_ALT = re.compile(r'(\w+):(\d+)')

# Genefuse.
def parse_genefuse_label(label):
	fusion = Fusion('genefuse')
	# Fusion breakpoint.
	sleft, sright = label.split('___')
	m = re.search(RE_GENEFUSE_LABEL, sleft)
	if m:
		region = Region(m.group(1), m.group(2), 'left')
		fusion.setRegion(region)
	m = re.search(RE_GENEFUSE_LABEL, sright)
	if m:
		region = Region(m.group(1), m.group(2), 'right')
		fusion.setRegion(region)
	# Evidence.
	mev = re.search(RE_GENEFUSE_TOTAL, label)
	if mev:
		evidence = mev.group(1)
	fusion.evidence = evidence
	return fusion

def read_genefuse_json(filename):
	data = read_json(filename)
	fusions = []
	for label in data.get('fusions', []).keys():
		if not label.lower().startswith('fusion'):
			continue
		fusions.append(parse_genefuse_label(label))
	return fusions
	
def read_genefuse_html(filename):
	fusions = []
	# Read.
	bhtml = ""
	with open(filename) as fid:
		bhtml = fid.read()
	soup = bs4.BeautifulSoup(bhtml, features="lxml")

	# Parsing.
	div_menu = soup.find("div", attrs={"id":"menu"})
	for record in div_menu.findAll("a"):
		label = record.text
		if not "fusion" in label.lower():
			continue
		fusions.append(parse_genefuse_label(label))
	return fusions

# Lumpy.
def read_lumpy(flumpy):
	fusions = []
	treats = set()

	vcf_in = VariantFile(flumpy)
	for record in vcf_in.fetch():
		# Check if variant is already seen.
		if 'SVTYPE' in record.info.keys() and record.info.get('SVTYPE') != 'BND':
			continue
		ident_number, ident_paired = record.id.split('_')
		if ident_number in treats:	
			continue
		treats.add(ident_number)

		# Build fusion.
		fusion = Fusion('lumpy')

		region = Region(record.chrom, int(record.pos)-1)
		fusion.setRegion(region)

		alt = record.alts[0]
		m = re.search(RE_LUMPY_ALT, alt)
		alt_chrom, alt_pos = 0, 0
		if m:
			alt_chrom = m.group(1)
			alt_pos = int(m.group(2))-1
		region = Region(alt_chrom, alt_pos)
		fusion.setRegion(region)

		evidence = 0
		if 'SU' in record.info.keys():
			evidence = record.info.get('SU')[0]
		fusion.setEvidence(evidence)

		fusions.append(fusion)
	return fusions

# Consensus.
def consensus_lumpy(records, consensus_interval):
	if len(records) == 0:
		return records

	fusions = []
		
	# Sort by evidence.
	records.sort(reverse=True)
	fusions.append(records[0])
	del records[0]

	# Merge records.
	for fusion in records:		
		if fusion.is_near(fusions[0], consensus_interval):
			fusions[0].evidence += fusion.evidence
		else:
			fusions.append(fusion)
	return fusions

def consensus_genefuse_lumpy(data, consensus_interval):
	lumpy = data['lumpy']['consensus']
	genefuse = data['genefuse']

	if len(lumpy) == 0:
		return genefuse
	if len(genefuse) == 0:
		return lumpy

	genefuse.sort(reverse=True)
	to_delete = []
	for gfusion in genefuse:
		for ix in range(0, len(lumpy)):
			if gfusion.is_near(lumpy[ix], consensus_interval):
				gfusion.evidence = (gfusion.evidence + lumpy[ix].evidence) / 2
				gfusion.software = 'consensus'
				to_delete.append(ix)
		# Update.
		lumpy = update_list(lumpy, to_delete)

	if len(lumpy) > 0:
		genefuse += lumpy

	# Filter if fusion is on same chrom.
	to_delete = []
	for ix in range(len(genefuse)):
		if genefuse[ix].is_same_chrom():
			to_delete.append(ix)
	genefuse = update_list(genefuse, to_delete)

	genefuse.sort(reverse=True)	
	return genefuse

# Write.
def write(filename, finputs, fusions):
	data = {}
	data['inputs'] = finputs
	data['fusions'] = {}
	for ix, fusion in enumerate(fusions['consensus']):
		data['fusions'][str(ix)] = fusion.to_dict()
	write_json(filename, data)
