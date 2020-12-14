import json
import logging
import os
import sys

from pkg_resources import parse_version
from platform import python_version

class Region():

	def __init__(self, chrom='', position=0, orientation='undefined'):
		self.setChrom(chrom)
		self.setPosition(position)
		self.setOrientation(orientation)

	def getChrom(self):
		return self._chrom
		
	def setChrom(self, chrom):
		self._chrom = chrom
		
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

	def is_init(self):
		if self._chrom == '':
			return False
		return True

	def to_dict(self):
		return dict(orientation=self._orientation,
			chrom=self._chrom,
			position=str(self._position))

	@classmethod
	def from_dict(cls, data):
		region = Region()
		region.orientation = data.get('orientation', '')
		region.chrom = data.get('chrom', '')
		region.position = data.get('position', 0)
		return region
		
	def __repr__(self):
		return '%s %s:%s'%(self._orientation, self._chrom, self._position)

	def __eq__(self, other):
		return self._chrom == other.chrom and self._position == other.position and self._orientation == other.orientation

	chrom = property(getChrom, setChrom)
	position = property(getPosition, setPosition)
	orientation = property(getOrientation, setOrientation)
	
class Fusion():
	
	name_consensus = 'CONS'

	def __init__(self, software='undefined'):
		self._first = Region()
		self._second = Region()
		self._evidence = 0
		self._evidence_details = {}
		self._depth = 0
		self._ident = ''
		self._buildFrom = set()
		self._isConsensus = False
		self._software = set()
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

	def getEvidenceDetails(self):
		return self._evidence_details

	def setEvidenceDetails(self, evidence_details):
		self._evidence_details = evidence_details.copy()

	def getDepth(self):
		return self._depth

	def setDepth(self, depth):
		self._depth = int(depth)

	def getIdent(self):
		return self._ident
	
	def setIdent(self, ident):
		if type(ident) == int:
			self.removeSoftware('undefined')
			name = ''
			if self._isConsensus:
				name = Fusion.name_consensus
			else:
				name = '-'.join([ x[:3].upper() for x in self._software])
			self._ident = name + str(ident)
		else:
			self._ident = ident

	def getFrom(self):
		return self._buildFrom

	def setFrom(self, ident):
		if type(ident) == list or type(ident) == set:
			self._buildFrom = list2set(ident)
		else:
			self._buildFrom.add(ident)

	def getIsConsensus(self):
		return self._isConsensus
	
	def setIsConsensus(self, value):
		self._isConsensus = value

	def removeFomNameCons(self):
		self._buildFrom = set([ x for x in self._buildFrom if not x.startswith(Fusion.name_consensus)])

	def getSoftware(self):
		return self._software

	def setSoftware(self, software):
		if type(software) is str:
			self._software.add(software)
		else:
			for x in software:
				self._software.add(x)

	def removeSoftware(self, software):
		if software in self._software:
			self._software.remove(software)

	def setRegion(self, region):
		if not self._first.is_init():
			self._first = region
		else:
			self._second = region

	def is_near(self, other, consensus_interval):
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

	def is_same_chrom(self):
		if self._first.is_init() and self._second.is_init():
			if self._first.chrom == self._second.chrom:
				return True
		return False

	def swap(self):
		tmp = self._first
		self._first = self._second
		self.second = tmp

	def vaf(self):
		return self._depth / self._evidence * 100

	def to_dict(self):
		return dict(software=list(self._software), 
				first=self._first.to_dict(), 
				second=self._second.to_dict(), 
				evidence=str(self._evidence),
				evidence_details=self._evidence_details,
				ident=self._ident, 
				buildFrom=list(self._buildFrom),
				isConsensus=self._isConsensus,
				depth=str(self._depth))

	@classmethod
	def from_dict(cls, data):
		fusion = Fusion()
		softwares = data.get('software', [])
		if len(softwares) > 0:		
			fusion.setSoftware(softwares)
			fusion.removeSoftware('undefined')
		fusion.first = Region.from_dict(data.get('first', {}))
		fusion.second = Region.from_dict(data.get('second', {}))
		fusion.evidence = data.get('evidence', 0)
		fusion.evidence_details = data.get('evidence_details', {})
		fusion.depth = data.get('depth', 0)
		fusion.ident = data.get('ident', '')
		fusion.isConsensus = data.get('isConsensus', True)
		fusion.buildFrom = list2set(data.get('buildFrom', []))
		return fusion


	def __repr__(self):
		return 'From %s %s %s (%s) %s %s Evidence %s Depth %s'%(','.join(self._software), self._ident, ','.join(self._buildFrom), self._evidence, self._first, self._second, self._evidence, self._depth)
 
	def __eq__(self, other):
		
		return self._software.difference(other.software) == 0 and other.software.difference(self._software) == 0 and self._ident == other.ident and self._first == other.first and self._second == other.second and self._evidence == other.evidence

	def __lt__(self, other):
		return self._evidence < other.evidence

	def __gt__(self, other):
		return self._evidence > other.evidence

	first = property(getFirst, setFirst)
	second = property(getSecond, setSecond)
	evidence = property(getEvidence, setEvidence)
	evidence_details = property(getEvidenceDetails, setEvidenceDetails)
	depth = property(getDepth, setDepth)
	ident = property(getIdent, setIdent)
	buildFrom = property(getFrom, setFrom)
	isConsensus = property(getIsConsensus, setIsConsensus)
	software = property(getSoftware, setSoftware)	

def abort(parser, msg=""):
	parser.error(msg)

def read_json(path):
	with open(path) as fid:
		return json.load(fid)

def write_json(path, data):
	with open(path, 'w') as fod:
		json.dump(data, fod)

def dict2list(dico):
	l = []
	for score in sorted(dico.keys()):
		l += [score] * dico[score]
	return l

def calculate_percentile(dico, target):
	total = 0
	ct = (target * (sum(dico.values())+1))/100 
	for key, value in sorted(dico.items(), key=lambda x:x[0]):
		total += value
		if total >= ct:
			return key
	return -1

def isPython3():
	return parse_version(python_version()) >= parse_version('3.0.0')

def isNaN(num):
	return num != num

def update_list(li, indexes):
	up = 0
	for ix in indexes:
		del li[ix-up]
		up += 1
	return li

def list2set(li):
	a = set()
	for i in li:
		a.add(i)
	return a
