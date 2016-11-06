# scaleLine.py - Parker de Waal 2016
# This tool helps to extend cylinders drawn in VMD
#
# Usage: ipython -i scaleLine.py
# To increase the distance between point A and B by 1.5 times: extendB(1.5) or extendA(1.5)

from __future__ import print_function
import numpy as np

# Set user arrays here:
A = np.array([-2.478,18.258,28.034])
B = np.array([37.912,31.968,27.864])


def extendB(scalar):
	array=(B-A)*scalar+A
	print(*array, sep=' ')
	print(*array, sep=',')

def extendA(scalar):
	array=(A-B)*scalar+B
        print(*array, sep=' ')
	print(*array, sep=',')
