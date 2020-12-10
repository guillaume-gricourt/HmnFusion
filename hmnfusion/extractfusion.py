import bs4
import copy
import json
import logging
import os
import pysam
import re
import sys

from pysam import VariantFile

from .utils import read_json, update_list, write_json, Region, Fusion
from ._version import __app_name__

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
	for ix, fusion in enumerate(fusions):
		fusion.ident = ix+1
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
	for ix, fusion in enumerate(fusions):
		fusion.ident = ix+1
	return fusions

# Lumpy.
def read_lumpy(flumpy):
	fusions = []
	treats = set()

	vcf_in = VariantFile(flumpy)
	ident = 1
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
		fusion.evidence = evidence

		fusion.ident = ident
		ident += 1
		fusions.append(fusion)
	return fusions

# Consensus.
def consensus_single(records, consensus_interval):
	if len(records) == 0:
		return records

	fusions = []
	# Sort by evidence.
	records.sort(reverse=True)
	merging = {}
	# Merge records.
	for ia in range(len(records)):
		for ib in range(ia+1, len(records)):
			if records[ia].is_near(records[ib], consensus_interval):
				if not ia in merging.keys():
					merging[ia] = set() 
				merging[ia].add(ib)

	# Merge non adjency fusions.
	isNotConvergence = True
	maxi = max(list(merging.keys()) + [0]) + 1
	while isNotConvergence:
		todel = set()
		isFound = False
		for ia in range(maxi):
			for ib in range(maxi):
				if ia == ib or not ia in merging.keys() or not ib in merging.keys():
					continue
				ix = set()
				ix.add(ib)
				for j in merging.get(ib, set()):
					ix.add(j)
				for i in ix:
					if ib in merging.get(ia, set()):
						for j in merging.get(ib, set()):
							merging[ia].add(j)
						del merging[ib]
						isFound = True
						break
				
		isNotConvergence = isFound
	# Build consensus
	for k, v in merging.items():
		ones = []
		for i in [k] + list(v):
			ones.append(copy.deepcopy(records[i]))
		ones.sort(reverse=True)

		one = ones.pop(0)
		one.buildFrom = one.ident
		for fusion in ones:
			if ('genefuse' in one.software and 'genefuse' in fusion.software) or ('lumpy' in one.software and 'lumpy' in fusion.software):
				one.evidence += fusion.evidence
			one.buildFrom = fusion.ident
			one.software = fusion.software
		one.isConsensus = True
		one.software = __app_name__
		one.removeBuildNameCons()
		fusions.append(one)
	fusions.sort(reverse=True)
	for ix, fusion in enumerate(fusions):
		fusion.ident = ix + 1
	return fusions

def consensus_genefuse_lumpy(genefuse_raw, lumpy_raw, genefuse_consensus, lumpy_consensus, consensus_interval):

	genefuse_not_shown = []
	for fusion_consensus in genefuse_consensus:
		for fusion_raw in genefuse_raw:
			if not fusion_raw.ident in fusion_consensus.buildFrom:
				genefuse_not_shown.append(fusion_raw)

	lumpy_not_shown = []
	for fusion_consensus in lumpy_consensus:
		for fusion_raw in genefuse_raw:
			if not fusion_raw.ident in fusion_consensus.buildFrom:
				lumpy_not_shown.append(fusion_raw)

	#consensus_raw = consensus_single(genefuse_not_shown + lumpy_not_shown, consensus_interval, 'mean')
	consensus = consensus_single(genefuse_not_shown + lumpy_not_shown + genefuse_consensus + lumpy_consensus, consensus_interval)
	consensus.sort(reverse=True)
	return consensus

def filter_same_chrom(fusions):
	# Filter if fusion is on same chrom.
	to_delete = []
	for ix in range(len(fusions)):
		if fusions[ix].is_same_chrom():
			to_delete.append(ix)
	fusions = update_list(fusions, to_delete)
	return fusions

# Write.
def write(filename, finputs, fusions):
	data = {}
	data['inputs'] = finputs
	data['fusions'] = {}
	for ix, fusion in enumerate(fusions['genefuse']['raw'] + fusions['lumpy']['raw'] + fusions['consensus']):
		data['fusions'][str(ix)] = fusion.to_dict()
	write_json(filename, data)
