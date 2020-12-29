
import random

q = 35
RRR = Integers(q)  # ring of Z/qZ

def AG_1dim(ES):
	""" Test Arora-Ge over composite molulus q. 
		1 is the LWE secret dimension, ES is the support of error, size = size of ES.  """
	Rx.<x> = PolynomialRing(RRR)
	s = Integers(q).random_element()
	size_ES = len(ES)
	m = 100*(2**size_ES)  # number of samples
	print("s = ", s, "m = ", m, "ES = ", ES)

	M_eq = []   # init the matrix of the linear equations
	y_eq = []   # init the matrix of the linear equations

	for i in range(m):
		# create an LWE sample
		a = Integers(q).random_element()
		e = ES[Integers(size_ES).random_element()]
		y = s*a+e  
		# print(a, e, y)
		
		# build a polynomial poly(x) = prod_{e\in ES}(y - ax - e)
		poly = Rx(1)
		for j in range(size_ES):
			poly *= Rx(y - a*x - ES[j])
		#print(poly)

		# linearization
		if poly != 0:
			L = [ poly[k] for k in range(1,size_ES+1) ]
			if L not in M_eq:
				M_eq.append(L)
				y_eq.append( [ -poly[0] ] )
	print(M_eq)
	print(y_eq)

	# build a square matrix of dim size_ES out of possibly many samples

	M = Matrix(RRR, [[ M_eq[j][i] for i in range(size_ES)] for j in range(size_ES)] )
	#print(M, M.det())
	row = 0  # starting point, will increase if the first choice is not invertible mod q
	while(gcd(M.det(), q) !=1 and row+size_ES<len(y_eq)):
		row+=1
		M = Matrix(RRR, [[ M_eq[j][i] for i in range(size_ES)] for j in range(row, row+size_ES)] )
	if row+size_ES<len(y_eq):
		Y = Matrix(RRR, [[ y_eq[j][0] ] for j in range(row, row+size_ES)] )
		print("Y = ", Y)
		print("[ s^1, s^2, ..., s^k] = ", M.inverse()*Y)
	else:
		print("not enough samples")


ES = [2,1,14]
AG_1dim(ES)

