#!/usr/bin/python3

import argparse
import sys
import gzip

parser = argparse.ArgumentParser(description="This script compares genetic variants in two individual samples from their VCF files. This comparison is executed only if a given variant is genotyped in both samples, and the genotype quality is above the 'minQ' threshold. The VCF files must contain only one sample each. They must be gzip compressed. The result is a table with physical coordinates for every variant and the type of observed difference between the samples: no difference, homozygous difference or heterozygous difference. This script writes to stdout.")
parser.add_argument("VCF1", type=str, help="VCF file for sample1.")
parser.add_argument("VCF2", type=str, help="VCF file for sample2.")
arguments = parser.parse_args()

# inputCmpVCF = arguments.VCF1
# outputFile = arguments.VCF2

def getLineInfo(sample,minQ=60):
	GQ = int(sample[-1].split(':')[sample[8].split(':').index('GQ')])
	if GQ <= minQ or sample[0].startswith('scaffold'):
		return None
	GT = sample[-1].split(':')[sample[8].split(':').index('GT')].split('/')
	REFALT = [i for i in sample[4].split(',')]
	REFALT.insert(0, sample[3])
	GT[0], GT[1] = REFALT[int(GT[0])], REFALT[int(GT[1])]
	sample[0] = sample[0].replace('Chr', '')
	return [sample[0],sample[1],GT,int(GQ)]

def addZygocity(line1, line2):
	line1 = getLineInfo(line1)
	line2 = getLineInfo(line2)
	if line1 and line2:
		mySet1 = set(line1[2])
		mySet2 = set(line2[2])
		line1[2] = line1[2][0] + '/' + line1[2][1]
		line2[2] = line2[2][0] + '/' + line2[2][1]
		line1[3] = line2[2]
		if mySet1 == mySet2:
			line1.append('noDiff')
			return line1
		elif len(mySet1) != len(mySet2):
			line1.append('HeteroDiff')
			return line1
		elif (len(mySet1) != 1) or (len(mySet2) != 1):
			line1.append('HeteroDiff')
			return line1
		else:
			line1.append('HomoDiff')
			return line1
	else:
		None

def getCommonSNPs(inputVCF1, inputVCF2):
	def getLine(string):
		return string.readline().decode().strip().split('\t')
	print('chromosome','position','sample1','sample2','diffType', sep='\t')
	with gzip.open(inputVCF1,'r') as inFile1, gzip.open(inputVCF2, 'r') as inFile2:
		theLine1 = getLine(inFile1)
		theLine2 = getLine(inFile2)
		while theLine1 or theLine2:
			if theLine1 == [''] or theLine2 == ['']:
				break
			elif (theLine1[0][0] == '#') or (theLine2[0][0] == '#'):
				theLine1 = getLine(inFile1)
				theLine2 = getLine(inFile2)
				continue
			elif theLine1[0] == theLine2[0]:
				currentChrom = theLine1[0]
				if int(theLine1[1]) == int(theLine2[1]):
					list2print = addZygocity(theLine1, theLine2)
					print(*list2print, sep='\t') if list2print else None
					theLine1 = getLine(inFile1)
					theLine2 = getLine(inFile2)
					continue
				elif int(theLine1[1]) > int(theLine2[1]):
					theLine2 = getLine(inFile2)
					continue
				elif int(theLine1[1]) < int(theLine2[1]):
					theLine1 = getLine(inFile1)
					continue
			elif theLine2[0] != currentChrom:
				theLine1 = getLine(inFile1)
				continue
			elif theLine1[0] != currentChrom:
				theLine2 = getLine(inFile2)
				continue
			else:
				break
		return None

try:
	getCommonSNPs(arguments.VCF1, arguments.VCF2)
except IOError:
	try:
		sys.stdout.close()
	except IOError:
		pass
	try:
		sys.stderr.close()
	except IOError:
		pass