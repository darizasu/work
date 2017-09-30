#!/usr/bin/python3

import argparse
import sys


parser = argparse.ArgumentParser(description="This script produces a simple tab-separated file with merged coordinates for repeats produced by NGSEP for different samples")
parser.add_argument("input_GFF3", type=str, help="GFF3 file with combined and sorted repeats (only repeats) produced by NGSEP.")
parser.add_argument("minQ", type=int, help="Minimum quality for the repeat that will be merged (from 0 to 40)")
parser.add_argument("output_file", type=str, help="Name of output GFF3 file with the genes that are contained in the coordinates specified in the coordinates_file")
arguments = parser.parse_args()

inputGFF3 = arguments.input_GFF3
minQual = arguments.minQ
output_file = arguments.output_file

def mergeReps(inGFF3, outFile, minQual):
	with open(inGFF3) as inputGFF3, open(outFile, "a") as outputFile:

		myDict = dict()
		print ("\nMerging regions in {} input file\n".format(inGFF3))

		for line in inputGFF3:
			theLine = line.strip().split('\t')

			if (int(theLine[5]) < minQual) or (theLine[2] != 'REPEAT'):
				continue

			else:

				if myDict.get(theLine[0]):
					start = myDict[theLine[0]][-1][0]
					end = myDict[theLine[0]][-1][1] + 1
					theRange = range(start, end)

					if (int(theLine[3]) in theRange) and (int(theLine[4]) in theRange):
						continue
					elif (int(theLine[3]) in theRange) and (int(theLine[4]) not in theRange):
						myDict[theLine[0]][-1][1] = int(theLine[4])
					elif (int(theLine[3]) not in theRange) and (int(theLine[4]) not in theRange):
						myDict[theLine[0]].append([int(theLine[3]), int(theLine[4])])
					else:
						continue

				else:
					myDict[theLine[0]] = [[int(theLine[3]), int(theLine[4])]]

		print ("\nWriting in {} output file now\n".format(outFile))
		for key, values in myDict.items():
			for value in values:
				outputFile.writelines([key, '\t', str(value[0]), '\t', str(value[1]), '\n'])

	print ("\nThe process seems to be completed\n")

mergeReps(inputGFF3, output_file, minQual)