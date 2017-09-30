#!/usr/bin/python3

import argparse
import sys


parser = argparse.ArgumentParser(description="This script produces a gff3 file with the genes that are contained in the coordinates specified in coordinates_file.")
parser.add_argument("input_GFF3", type=str, help="Reference GFF3")
parser.add_argument("coordinates_file", type=str, help="File with one region per line, every line containing chromosome_name, start, end positions separated by tab (the whole file should have 3 columns)")
parser.add_argument("output_GFF3", type=str, help="Name of output GFF3 file with the genes that are contained in the coordinates specified in the coordinates_file")
arguments = parser.parse_args()

inputGFF3 = arguments.input_GFF3
inputTxt = arguments.coordinates_file
outputGFF3 = arguments.output_GFF3

def getGenes(inGFF3, inTxt, outGFF3):
	with open(inGFF3) as inputGFF3, open(inTxt) as inputTxt, open(outGFF3, "a") as outputGFF3:

		myDict = dict()

		for coord in inputTxt:
			theLine = coord.strip().split('\t')
			if myDict.get(theLine[0]):
				myDict[theLine[0]].append([int(theLine[1]), int(theLine[2])])
			else:
				myDict[theLine[0]] = [[int(theLine[1]), int(theLine[2])]]
				continue

		for line in inputGFF3:
			theLine = line.strip().split('\t')
			if theLine[0][0] == '#':
				outputGFF3.write(line)
				continue
			elif theLine[2] == 'gene' and myDict.get(theLine[0]):
				for element in myDict.get(theLine[0]):
						theRange = range(int(theLine[3]), int(theLine[4]) + 1)
						start = element[0]
						end = element[1]
						if (start in theRange) or (end in theRange):
#							outputGFF3.writelines(['\t'.join(i) + '\n' for i in [theLine]])
							outputGFF3.write(line)
						else:
							continue
			else:
				continue

		print ("\n Completed! The output file was named {} \n".format(outGFF3))

getGenes(inputGFF3,inputTxt,outputGFF3)
