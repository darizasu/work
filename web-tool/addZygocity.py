#!/usr/bin/python3

import sys

try:
	for line in sys.stdin:
		theLine = line.strip().split('\t')
		name = theLine[4].split('=')[0] + '-' + theLine[5].split('=')[0]
		if len(theLine) < 2:
			name = theLine[0]
			continue
		elif theLine[4].split('=')[1] == theLine[5].split('=')[1]:
			print(line.strip(),'noDiff',name,sep='\t',end='\n')
			continue
		else:
			mySet1 = {i for i in theLine[4].split('=')[1].split('/')}
			mySet2 = {i for i in theLine[5].split('=')[1].split('/')}
			if len(mySet1) != len(mySet2):
				print(line.strip(),'HeteroDiff',name,sep='\t',end='\n')
				continue
			elif (len(mySet1) != 1) or (len(mySet2) != 1):
				print(line.strip(),'HeteroDiff',name,sep='\t',end='\n')
				continue
			else:
				print(line.strip(),'HomoDiff',name,sep='\t',end='\n')
				continue
except IOError:
	try:
		sys.stdout.close()
	except IOError:
		pass
	try:
		sys.stderr.close()
	except IOError:
		pass