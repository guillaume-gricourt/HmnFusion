import json
import logging
import pysam

import numpy as np 
import pandas as pd

from .utils import read_json, update_list, write_json, Region, Fusion

# Bed.
def read_bed(filename):
	isHeader = 0
	with open(filename) as fid:
		if 'track' in fid.readline():
			isHeader=1
	bed = pd.read_csv(filename, sep='\t', usecols=[0,1,2], names=["chrom", "start", "end"], dtype={"chrom":str, "start":int, "end":int}, skiprows=isHeader)
	return bed

# Parsing regions.
def _parse_region(region):
	m = region.split(':')
	chrom, pos = '', 0
	try:
		chrom = m[0]
		pos = int(m[1])
	except:
		pass
	return (chrom, pos)

def check_region(sregion):
	if not ':' in sregion:
		return False
	if len(sregion.split(':')) != 2:
		return False 
	region = _parse_region(sregion)
	if region[1] == 0:
		return False
	return True

def build_region(sregion):
	fusion = Fusion()
	fusion.first.chrom, fusion.first.position = _parse_region(sregion)
	return fusion

def parse_hmnfusion_json(filename):
	data = read_json(filename)
	fusions = []
	for key in sorted(data["fusions"].keys()):
		fusions.append(Fusion.from_dict(data["fusions"][key]))
	return fusions

def _cigar2position(cigars, start):
	data = {}
	if cigars[0][0] in [4, 5]:
		for i in range(cigars[0][1]):
			i += 1
			data[start-i] = cigars[0][0]
		del cigars[0]
	for cigar in cigars:
		op, nb = cigar[0], cigar[1]
		if op in [4, 5]:
			for i in range(nb):
				data[start+i] = op
				i += 1
		elif op == 0:
			for i in range(nb):
				data[start] = op
				start+=1
		elif op == 1:
			pass
		else:
			start += 1
	return data


def _select_bed(x, region):
	if x['chrom'] == region.chrom:
		if x['start'] <= region.position and x['end'] >= region.position:
			return True
	return False
		

def run(params, bed, fusions):

	alignment = pysam.AlignmentFile(params['falignment']['path'], params['falignment']['mode'])

	to_delete = []
	isSkip = False
	for ix, fusion in enumerate(fusions):

		# Check fusion against bed.
		sub_first, sub_second = pd.DataFrame(columns=bed.columns), pd.DataFrame(columns=bed.columns)
		if fusion.first.is_init():
			sel = bed.apply(_select_bed, axis=1, args=(fusion.first,))
			sub_first = bed[sel]
		if fusion.second.is_init():
			sel = bed.apply(_select_bed, axis=1, args=(fusion.second,))
			sub_second = bed[sel]
		if len(sub_first) > 1 or len(sub_second) > 1:
			logging.warning('Fusion %s is found multiple times in bed -> skipping'%(fusion,))
			isSkip = True
		if len(sub_first) + len(sub_second) == 2 :
			logging.warning('Fusion %s is found on left and right of breakpoint in the bed -> skipping'%(fusion,))
			isSkip = True
		if len(sub_first) + len(sub_second) == 0 :
			logging.warning("Fusion %s isn't found on left or right of breakpoint in the bed -> skipping"%(fusion,))
			isSkip = True

		if isSkip:
			to_delete.append(ix)
			continue
		# Init.
		bed_sel = pd.DataFrame(columns=bed.columns)
		region = Region()
		if len(sub_first) == 1:
			bed_sel = sub_first
			region = fusion.first
		elif len(sub_second) == 1:
			bed_sel = sub_second
			region = fusion.second
		else:
			logging.warning("Fusion %s, something bad happened -> skipping"%(fusion,))
		count = dict(coverage=0, split=0, mate=0, clipped=0)
		# Run. 
		for aligned_segment in alignment.fetch(bed_sel.iloc[0, 0], bed_sel.iloc[0, 1], bed_sel.iloc[0, 2]):
			# Filtering.
			if aligned_segment.is_unmapped or aligned_segment.is_duplicate or aligned_segment.is_supplementary:
				continue

			cigar2pos = _cigar2position(aligned_segment.cigartuples, aligned_segment.reference_start)				
			if not region.position in cigar2pos.keys(): 
				continue

			count['coverage'] += 1
			# Count split reads.
			if aligned_segment.has_tag("SA"):
				count['split'] += 1
				continue
			
			# Count other Chrom.
			if aligned_segment.is_paired:	
				if not aligned_segment.mate_is_unmapped and not aligned_segment.is_unmapped:
					if aligned_segment.next_reference_id != aligned_segment.reference_id:
						count['mate'] += 1
						continue
		
			# Count reads clipped.
			count_clipped = np.zeros((2, params['clipped']['interval']))
			for i in range(params['clipped']['interval']):
				if cigar2pos.get(region.position-i-1,0) in [4, 5]:
					count_clipped[0][i] = 1
				if cigar2pos.get(region.position+i+1,0) in [4, 5]:
					count_clipped[1][i] = 1
			
			if np.max(np.sum(count_clipped, axis=1)) >= params['clipped']['count']:
				count['clipped'] += 1
		fusion.depth = count['coverage']
		fusion.evidence = count['split'] + count['mate'] + count['clipped']
	fusions = update_list(fusions, to_delete)
	return fusions

# Write.
def write(filename, finputs, params, fusions):
	data = {}
	data['inputs'] = finputs
	data['parameters'] = params
	data['fusions'] = {}
	for ix, fusion in enumerate(fusions):
		data['fusions'][str(ix)] = fusion.to_dict()
	write_json(filename, data)
