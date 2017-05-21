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

def getLineInfo(line,minQ=40):
	GQ = int(line[-1].split(':')[line[8].split(':').index('GQ')])
	# if GQ <= minQ: #or line[0].startswith('scaffold'):
		# return None
	GT = line[-1].split(':')[line[8].split(':').index('GT')].split('/')
	REFALT = [i for i in line[4].split(',')]
	REFALT.insert(0, line[3])
	GT[0], GT[1] = REFALT[int(GT[0])], REFALT[int(GT[1])]
	# line[0] = line[0].replace('Chr', '')
	return [line[0] + ':' + line[1],GT,int(GQ)]

def addZygocity(line1, line2):
	if line1 and line2:
		mySet1 = set(line1[1])
		mySet2 = set(line2[1])
		line1[1] = line1[1][0] + '/' + line1[1][1]
		line2[1] = line2[1][0] + '/' + line2[1][1]
		line1[2] = line2[1]
		line1 = line1[0].split(':') + line1[1:]
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
	print('chromosome','position','sample1','sample2','diffType', sep='\t')
	with gzip.open(inputVCF1,'r') as inFile1, gzip.open(inputVCF2, 'r') as inFile2:
		myDict = dict()
		for line1 in inFile1:
			theLine1 = line1.decode().strip().split('\t')
			if theLine1[0][0] == '#':
				continue
			else:
				theLine1 = getLineInfo(theLine1)
				if theLine1 and not myDict.get(theLine1[0]):
					myDict[theLine1[0]] = [theLine1]
				elif theLine1 and myDict.get(theLine1[0]):
					myDict[theLine1[0]].append(theLine1)
		for line2 in inFile2:
			theLine2 = line2.decode().strip().split('\t')
			if theLine2[0][0] == '#':
				continue
			else:
				theLine2 = getLineInfo(theLine2)
				if theLine2 and myDict.get(theLine2[0]):
					list2print = addZygocity(myDict[theLine2[0]].pop(0), theLine2)
					print(*list2print, sep='\t') if list2print else None
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
