import logging
import os
import pysam
import sys

from .utils import readJson

# Genefuse.
def read_genefuse_json(filename, ):
	data = readJson(filename)
	for fusion in data.get('fusions', []).keys():
		df
		
def read_genefuse_json(filename, ):



def read_genefuse(fgenefuse):
	if fgenefuse['format'] == 'json':
		read_genefuse_json(fgenefuse['path'], )
	elif fgenefuse['format'] == 'html':
		read_genefuse_html(fgenefuse['path', )
	
# Lumpy.
