import logging
import os
import pysam
import sys

from .utils import readJson

class Region():
	def __init__(self):
		self._chrom = None
		self._position = None
		
	def getChrom(self):
		return self._chrom
		
	def setChrom(self, chrom):
		self._chrom = chrom
		
	def getPosition(self):
		return self._position
		
	def setPosition(self, position):
		self._position = position
		
	chrom = property(getChrom, setChrom)
	position = property(getPosition, setPosition)
	
class Fusion():
	def __init__(self):
		self._left = Region()
		self._right = Region()
		
	def getLeft(self):
		return self._left
		
	def setLeft(self, left):
		self._left = left
		
	def getRight(self):
		return self._right
		
	def setRight(self, right):
		self._right = right
		
	left = property(getLeft, setLeft)
	right = property(getRight, setRight)
	
# Genefuse.
def read_genefuse_json(filename, ):
	data = readJson(filename)
	fusions = []
	for event in data.get('fusions', []).keys():
		if not event.lower().startswith('fusion'):
			continue
		fusion = Fusion()
		for side in ['left', 'right']:
			chrom, pos = None, None
			if side in data['fusions'][event].keys():
				if 'gene_chr' in data['fusions'][event][side].keys():
					mgene = re.search(r'\D*(\d+)\D*', data['fusions'][event][side]['gene_chr'])
					if mgene:
						chrom = mgene.group(1)
				if 'pos_str' in data['fusions'][event][side].keys():
					mpos = re.search(r'\w+:\w+:\w+\|[+-]\w+:(\d+)', data['fusions'][event][side]['pos_str'])
					if mpos:
						pos = mpos.group(1)
			if side == 'left':
				fusion.left.chrom = chrom
				fusion.left.position = pos
			if side == 'right':
				fusion.right.chrom = chrom
				fusion.right.position = pos
	return fusions
	
def read_genefuse_json(filename, ):



def read_genefuse(fgenefuse):
	fusions = []
	if fgenefuse['format'] == 'json':
		fusions = read_genefuse_json(fgenefuse['path'], )
	elif fgenefuse['format'] == 'html':
		fusions = read_genefuse_html(fgenefuse['path', )
	return fusions
# Lumpy.
