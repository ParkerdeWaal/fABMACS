# restraintGen.py - Parker de Waal 2016
# Generates a restraint file for CA atoms 
#
# Usage: python restraintGen.py {system .GRO file} {ligand resname} {restraint radius in nm}
#        python restraintGen.py isob.gro LIG 1 > chainA.itp

import sys
import numpy as np
from __future__ import print_function


inFile = sys.argv[1]
ligandName = sys.argv[2]
restraintRadius = sys.argv[3]
groFile = open(inFile, "r")

# Define new array types and initilze arrays/list
raw = np.dtype([('ID', int),('pos','(3,)f8')])
ligPos = np.empty([0,2],dtype=raw)
caPos = np.empty([0,2],dtype=raw)
resList = []

for line in groFile:
	if ligandName in line:
		ligPos = np.append(ligPos, np.array([(line.split()[2],(line.split()[3],line.split()[4], line.split()[5]))], dtype=raw))
	if line[13:15] == 'CA':
		caPos = np.append(caPos, np.array([(line.split()[2],(line.split()[3],line.split()[4], line.split()[5]))], dtype=raw))

# calculate distanes
for d1 in np.nditer(ligPos):
	for d2 in np.nditer(caPos):
		tempDist1 = (d1['pos']-d2['pos'])**2
		tempDist1 = tempDist1.sum(axis=-1)
		tempDist1 = np.sqrt(tempDist1)
		if tempDist1 >= float(restraintRadius):
			resList.append(int(d2['ID']))
			b = set(resList)
# sort list
b = [int(x) for x in b]
b.sort()

print """
[ position_restraints ]
;  i funct	 fcx        fcy        fcz
"""

for item in b:
	print '%5i  1 1000 1000 1000' % item
			

