# Parse the coefficients from https://math.mit.edu/~drew/ClassicalModPolys.html 
# e.g. from 
# [3,0] 1
# [2,0] -162000
# [2,1] 1488
# [2,2] -1
# [1,0] 8748000000
# [1,1] 40773375
# [0,0] -157464000000000
# to X3 - X2Y2 + 1488X2Y - 162000X2 + 1488XY2+ 4077375XY + 8748000000X + Y3 - 162000Y2 + 8748000000Y -157464000000000. 

import math
import sys
import time
import numpy as np
#from scipy import linalg

NUM_ROW = 7

def main():

	infile = open("phi_j_2")
	coeff = []
	for i in xrange(0,NUM_ROW):
		it = infile.readline()
		print i, it
		s = it.replace("[","").replace("\n","").replace("]",",").split(',')
		v = []
		for e in s:
			v.append(int(e))
		coeff.append( v )
	print len(coeff), coeff
			
	infile.close()
	
main()
