#!/usr/bin/python3

import argparse
import sys


parser = argparse.ArgumentParser(description="This script converts the output file from NGSEP-CompareVCF into a comparison matrix among samples. This script produces a single matrix with the information of one single column of the output file of CompareVCF. That single column could be 'SNPs_Both', '%HomoDifferences', 'TotalDifferences', etc. Check NGSEP-CompareVCF for more details")
parser.add_argument("CompareVCF_input", type=str, help="The output file from NGSEP-CompareVCF produced after comparing two VCF files")
parser.add_argument("colNumber", type=int, help="The number of the column that contains the information that you want to include in the matrix. First column is number '1', second column is number '2' and so on. Please specify only ONE number")
parser.add_argument("outputFile", type=str, help="The filename of the tab-separated comparison matrix. Be careful with the filename, as this script appends its output to the given filename")
arguments = parser.parse_args()

inputCmpVCF = arguments.CompareVCF_input
colNum = arguments.colNumber - 1
outputFile = arguments.outputFile

def CompareVCF2matrix(inFile, outFile, colNum):
	with open(inFile) as inputCmpVCF, open(outFile, "a") as outputFile:

		myList = list()
		mySet = set()
		myDict = dict()
		tab = 1
		print ("\nMerging regions in {} input file\n".format(inFile))

		for line in inputCmpVCF:
			theLine = line.strip().split('\t')
			if theLine[0] == 'Sample1':
				continue
			else:
					myList.append([theLine[0], theLine[1], theLine[colNum]])
					mySet.add(theLine[0])
		myList = sorted(myList, key=lambda x: (x[0], x[1]))
		mySet_keys = sorted(mySet)

		for key in mySet_keys:
			myDict[key] = tab
			tab += 1

		outList = [['s' for y in range(tab)] for x in range(tab)]
		outList[0][1:] = mySet_keys
		
		for obs in myList:
			if colNum == 8:
				outList[myDict[obs[0]]][myDict[obs[1]]] = str(100 - float(obs[2]))
			else:
				outList[myDict[obs[0]]][myDict[obs[1]]] = obs[2]
			outList[myDict[obs[0]]][0] = obs[0]

		for newLine in outList:
			outputFile.writelines(['\t'.join(i) + '\n' for i in [newLine]])

		return outList

CompareVCF2matrix(inputCmpVCF, outputFile, colNum)