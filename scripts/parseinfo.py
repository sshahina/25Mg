#!/bin/python
import csv
import pprint

with open('../info/17Oang_FixedRunSheet.dat','r') as csvfile:
	reader = csv.DictReader(csvfile,delimiter='\t')
	for row in reader:
		print(row['Run'],row['Ea (keV)'],row['Q(target)'])
	#	print(row)
	#readerList = list(reader)
	pprint.pprint(readerList)

	for entry in readerList:
	#	if entry['good?'] == '#':
			#print("\tskipping line")
	#	else:
			#print(entry['Run'])
	
	print(reader['Run'])#.index("10154"))

