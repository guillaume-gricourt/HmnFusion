import logging
import os
import sys

from pkg_resources import parse_version
from platform import python_version

def abort(help, msg=""):
	logging.error(msg)
	logging.info(help, format='%(message)s')
	sys.exit(1)

def readJson(path):
	with open(path) as fid:
		return json.load(fid)

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
