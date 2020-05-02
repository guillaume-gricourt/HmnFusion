import bs4
import json
import logging
import os
import pysam
import re
import sys

from pysam import VariantFile

from .utils import readJson

class Region():

	def __init__(self, chrom=0, position=0, orientation='undefined'):
		self.setChrom(chrom)
		self.setPosition(position)
		self.setOrientation(orientation)

	def getChrom(self):
		return self._chrom
		
	def setChrom(self, chrom):
		self._chrom = int(chrom)
		
	def getPosition(self):
		return self._position
		
	def setPosition(self, position):
		self._position = int(position)

	def getOrientation(self):
		return self._orientation
	
	def setOrientation(self, orientation):
		if orientation in ['left', 'right']:
			self._orientation = orientation
		else:
			self._orientation = 'undefined'

	def isDefined(self):
		if self._chrom == 0:
			return False
		return True

	def to_dict(self):
		return dict(orientation=self._orientation,
			chrom=str(self._chrom),
			position=str(self._position))

	def __repr__(self):
		return '%s %s:%s'%(self._orientation, self._chrom, self._position)

	def __eq__(self, other):
		return self._chrom == other.chrom and self._position == other.position and self._orientation == other.orientation

	chrom = property(getChrom, setChrom)
	position = property(getPosition, setPosition)
	orientation = property(getOrientation, setOrientation)
	
class Fusion():
	
	def __init__(self, software='undefined'):
		self._first = Region()
		self._second = Region()
		self._evidence = 0
		self.setSoftware(software)
	
	def getFirst(self):
		return self._first
		
	def setFirst(self, first):
		self._first = first
		
	def getSecond(self):
		return self._second
		
	def setSecond(self, second):
		self._second = second
	
	def getEvidence(self):
		return self._evidence

	def setEvidence(self, evidence):
		self._evidence = int(evidence)

	def getSoftware(self):
		return self._software

	def setSoftware(self, software):
		if software in ['genefuse', 'lumpy', 'consensus']:
			self._software = software
		else:
			self._software = 'undefined'

	def setRegion(self, region):
		if not self._first.isDefined():
			self._first = region
		else:
			self._second = region

	def isNear(self, other, consensus_interval):
		gaps = []
		if self._first.chrom == other.first.chrom:
			gaps.append(abs(self._first.position-other.first.position))
			if self._second.chrom == other.second.chrom:
				gaps.append(abs(self._second.position-other.second.position))
		elif self._first.chrom == other.second.chrom:
			gaps.append(abs(self._first.position-other.second.position))
			if self._second.chrom == other.first.chrom:
				gaps.append(abs(self._second.position-other.first.position))
		if len(gaps) == 2 and min(gaps) <= consensus_interval:
			return True
		return False

	def to_dict(self):
		return dict(software=self._software,
			first=self._first.to_dict(),
			second=self._second.to_dict(),
			evidence=str(self._evidence))

	def __repr__(self):
		return 'From %s (%s) %s %s'%(self._software, self._evidence, self._first, self._second)
 
	def __eq__(self, other):
		return self._software == other.software and self._first == other.first and self._second == other.second and self._evidence == other.evidence

	def __lt__(self, other):
		return self._evidence < other.evidence

	def __gt__(self, other):
		return self._evidence > other.evidence

	first = property(getFirst, setFirst)
	second = property(getSecond, setSecond)
	evidence = property(getEvidence, setEvidence)
	software = property(getSoftware, setSoftware)

RE_GENEFUSE_LABEL = re.compile(r'[+-]\D*(\d+):(\d+)')
RE_GENEFUSE_TOTAL = re.compile(r'total:\s{1}(\d+),')
RE_LUMPY_CHROM = re.compile(r'\D*(\d+)')
RE_LUMPY_ALT = re.compile(r'\D+(\d+):(\d+)')

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
	data = readJson(filename)
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

		ref_chrom = re.search(RE_LUMPY_CHROM, record.chrom).group(1)
		ref_pos = record.pos
		region = Region(ref_chrom, ref_pos)
		fusion.setRegion(region)

		alt = record.alts[0]
		m = re.search(RE_LUMPY_ALT, alt)
		alt_chrom, alt_pos = 0, 0
		if m:
			alt_chrom = m.group(1)
			alt_pos = m.group(2)
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
		if fusion.isNear(fusions[0], consensus_interval):
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
	for gfusion in genefuse:
		for ix in range(0, len(lumpy)):
			if gfusion.isNear(lumpy[ix], consensus_interval):
				gfusion.evidence = (gfusion.evidence + lumpy[ix].evidence) / 2
				gfusion.software = 'consensus'
				del lumpy[ix]
	if len(lumpy) > 0:
		genefuse += lumpy
	genefuse.sort(reverse=True)	
	return genefuse

# Write
def write(filename, finputs, fusions):
	data = {}
	data['inputs'] = finputs
	data['fusions'] = {}
	for ix, fusion in enumerate(fusions['consensus']):
		data['fusions'][str(ix)] = fusion.to_dict()
	with open(filename, 'w') as fod:
		json.dump(data, fod)
