#!/usr/bin/env python
# coding: utf8

import argparse
import copy
import io
import json
import os
import pdfkit
import pysam
import re
import sys
import matplotlib

import numpy as np 

##############
## Function ##
##############
def abort(help, msg=""):
	print(msg)
	print(help)
	sys.exit(1)

def checkFile(path):
	if os.path.exists(path):
		return True
	return False	
		
##################
## Parsing Args ##
##################

description = 'Create Pdf report metagenomic shotgun'
parser = argparse.ArgumentParser(description=description)

parser.add_argument('-i', required=True, dest='bam',
	action='store',help='Bam file')
parser.add_argument('-o', required=True, dest='output',
	action='store',help='File count txt')

if len(sys.argv) <= 1: 
	abort(parser.print_help() )
else: 
	args = parser.parse_args()

fbam = args.bam
foutput = args.output

BASESOFTCLIP = 4
BASESOFTCLIP_TOLERATE = 2

bam = pysam.AlignmentFile(fbam, "rb")

print(fbam)
data = {}


def cigar2position(cigars, start):
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
	
	

#bed_bcr = ("chr4", 54265897, 54280889)
#fusion_bcr = (54266930, 54266931)

bed_bcr = ("chr4", 55140698, 55141140)
fusion_bcr = (55141043, 55141044)


#bed_bcr = ("chr22", 23631700, 23634833)
#fusion_bcr = (23632549, 23632550)
#bed_abl = ("chr9", 133615450, 133615900)
#fusion_abl = (133615629, 133615630)

#cov = bam.count_coverage(contig=chrom, start=start, stop=stop)
#print("coverage %s" %(cov,))

coverage = 0
diffchr = 0
split_read = 0
softclip = 0

col = 0

dup = 0

already_count = set()

reads_soft_clips = []
data = {}
data["split"] = []
data["chrom"] = []
data["clip"] = []

reads = [ x for x in bam.fetch(bed_bcr[0], bed_bcr[1], bed_bcr[2]) if not x.is_duplicate and not x.is_unmapped]
#reads += [ x for x in bam.fetch(bed_abl[0], bed_abl[1], bed_abl[2]) if not x.is_duplicate and not x.is_unmapped]

for aligned_segment in reads: 
	
	reference_start = aligned_segment.reference_start
	cigars = aligned_segment.cigartuples
	cigar2pos = cigar2position(cigars, reference_start)
	
	sens = "2" if aligned_segment.is_read2 else "1"
	name_segment = aligned_segment.query_name + "-" + sens

	fusionSecondary = False
	start = 0
	'''
	if aligned_segment.reference_name == "chr9":
		start = fusion_abl[0]
		fusionSecondary = True
	elif aligned_segment.reference_name == "chr22":
		start = fusion_bcr[0]
	'''
	if aligned_segment.reference_name == "chr4":
		start = fusion_bcr[0]
	if not start in cigar2pos.keys(): 
		continue


	if aligned_segment.is_supplementary:
		continue
			
	coverage += 1


	isFound = False


	
	
	#Split
	if aligned_segment.has_tag("SA"):
		split_read += 1
		#print("split read " + name_segment)
		
		if name_segment in already_count:
			dup += 1
			print([ str(x) for x in reads if x.query_name == "M03869:187:000000000-CNMKC:1:1103:10744:20313"])
			print(aligned_segment)
			sys.exit(0)	
		#print(aligned_segment)
		data["split"] += [name_segment]
		already_count.add(name_segment)
		isFound = True

		continue
		
	#Other Chrom
	if not aligned_segment.is_secondary and not aligned_segment.is_supplementary and aligned_segment.is_paired:	
		if not aligned_segment.mate_is_unmapped and not aligned_segment.is_unmapped:
			if aligned_segment.next_reference_id != aligned_segment.reference_id:
				diffchr += 1
				#print("other " + name_segment)
				if name_segment in already_count:
					dup += 1
					print([ str(x) for x in reads if x.query_name == "M03869:187:000000000-CNMKC:1:1103:10744:20313"])
					print(aligned_segment)
					sys.exit(0)

				#print(name_segment)
				#print(aligned_segment)
				data["chrom"] += [name_segment]
				already_count.add(name_segment)
				isFound = True

				continue
	
	#Check clip		
	#reference_start = aligned_segment.reference_start
	reference_end = aligned_segment.reference_end
	reference_length = aligned_segment.reference_length
	
	#cigars = aligned_segment.cigartuples
	#cigar2pos = cigar2position(cigars, reference_start)
	
	#cigars = aligned_segment.cigartuples
	isSoftClip = []
	#Down
	softClipFound = False
	
	#if name_segment == "M03869:187:000000000-CNMKC:1:1113:3350:7281":
	#	print(cigar2pos)
		

	for i in range(BASESOFTCLIP + BASESOFTCLIP_TOLERATE):
		if cigar2pos.get(start-i-1,0) in [4, 5]:
			isSoftClip.append(True)
		else:
			isSoftClip.append(False)
	
	if sum(isSoftClip) >= BASESOFTCLIP:
		softclip += 1
		softClipFound = True
		#print("Down : " + aligned_segment.query_name)
	isSoftClip = []	
	if not softClipFound:
		for i in range(BASESOFTCLIP + BASESOFTCLIP_TOLERATE):
			if cigar2pos.get(start+i+2,0) in [4, 5]:
				isSoftClip.append(True)
			else:
				isSoftClip.append(False)
		if sum(isSoftClip) >= BASESOFTCLIP:
			softclip += 1
			#print("Up : " + aligned_segment.query_name)
			softClipFound = True
			
	
	if softClipFound:
		reads_soft_clips.append(name_segment)
		if name_segment in already_count:
			dup += 1
			print([ str(x) for x in reads if x.query_name == "M03869:187:000000000-CNMKC:1:1103:10744:20313"])
			print(aligned_segment)
			sys.exit(0)			
		already_count.add(name_segment)
		isFound = True
		

		
			
reads_soft_clips = sorted(reads_soft_clips)

#print('\n'.join(reads_soft_clips))
print("Coverage %s\n" % (coverage,))
print("SplitRead %s \nSoftClip %s \nDiffChr %s" % (split_read, softclip, diffchr))
print("Dup %s\n" % (dup,))
with open(foutput, "w") as fod:
	json.dump(data, fod, indent=4)

