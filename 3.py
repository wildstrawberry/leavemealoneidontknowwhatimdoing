# Brouver_boost_v_31
# Try boosting an abitrary (non-zero-sum) equilibria
# Previous: 
#   v_3, 2013.07.24, initialize
#   v_31, 2013.07.29, allowing non-zero-sum, simplify the column oracle
#   v_314, 2013.08.01, double-boosting
#   v_3141, 2013.08.06, plus random matrix generation
#   v_31415, 2013.08.07, giving up double boosting

import math
import psyco
import time
import sys
from pylab import *
#import numpy as np
import matplotlib.pyplot as plt
import random

psyco.full()

epsilon = 0.001

def main():
	
	# Notation: xAy: x wants to maximize his payoff; xBy: y wants to maximize via B.
	# when B=-A: zerosum

	M = 3	#number of rows
	N = 3	#number of cols	
	
#	print "payoff_A ="
#	payoff_A =	Random_matrix_generation(M, N, 0.0, 1.0)	
#	print "payoff_B ="
#	payoff_B =	Random_matrix_generation(M, N, 0.0, 1.0)	
	
#	payoff_A =	[[0.4	,0.6	, 1	],
#				[1		,0.2	,0.6	],
#				[0.6	,1		, 0	]]

#	payoff_B =	[[0.2,1, 0.6	],
#				[0.6	,0	,1	],
#				[1	,0.6	,0.4	]]


	payoff_A =	[[-2	,-4, 1	],
				[1		,-3	,-2	],
				[-3		,1	, -4	]]

	payoff_B =	[[-3,1	,-1	],
				[-1	,-4	,1	],
				[1	,0	,-2	]]

# do for zero sum game

#	for	i in xrange(0,M):				
#		for j in xrange(0,N):
#			payoff_B[i][j]=-payoff_A[i][j]

#	double_boost(M,N,payoff_A,payoff_B)
	oracle_boost(M,N,payoff_A,payoff_B)




def oracle_boost(M,N,payoff_A,payoff_B):	

	start = time.time()

#Uniform strategies initialization
	A = []
	B = []
	for i in xrange(0,M):
		A.append(1.0/M)
	for j in xrange(0,N):
		B.append(1.0/N)

	T = int(math.log(N)/epsilon/epsilon)
	percent = int(T/1000)
	
	print "Epsilon:", epsilon
	print "Number of rounds:", T
	print "Initialized x(1): ", A
	print "Initialized y(1): ", B
	
	# for taking average
	record_A = []
	record_B = []
	
	avg_payoff_A = 0
	avg_payoff_B = 0
	
	for i in xrange(0,M):
		record_A.append(0)
	for j in xrange(0,N):
		record_B.append(0)

	x_axis = []
	y0_axis = []
	y1_axis = []
	y2_axis = []

	x_value = []
	y_value = []

	for t in xrange(0,T):
		
		B = oracle_B(M,N,payoff_B,A,B)	

		A = update_A(M,N,payoff_A,A,B)	
		
		for i in xrange(0,M):
			record_A[i]+=A[i]
		for j in xrange(0,N):
			record_B[j]+=B[j]
			
		avg_payoff_A += Expected_payoff(M,N,payoff_A,A,B)
		avg_payoff_B += Expected_payoff(M,N,payoff_B,A,B)
		
		if (t%percent==1):				#show something intermediate
#		if (t<10000):				#show something intermediate
	#		print "round", t, A, B
			x_axis.append(t)
			y0_axis.append(A[0])
			y1_axis.append(A[1])
			y2_axis.append(A[2])
			A_t = normalize(M, record_A)
			B_t = normalize(N, record_B)
			x_value.append(Expected_payoff(M,N,payoff_A,A_t,B_t))
			y_value.append(Expected_payoff(M,N,payoff_B,A_t,B_t))

		
	print "x(T):", A
	print "y(T):", B
	print "Expected payoff A by x(T), y(T)", Expected_payoff(M,N,payoff_A,A,B)
	print "Expected payoff B by x(T), y(T)", Expected_payoff(M,N,payoff_B,A,B)

	A_avg = normalize(M, record_A)
	B_avg = normalize(N, record_B)

	print "sum A(x(T), y(T))/T:", avg_payoff_A/T
	print "sum B(x(T), y(T))/T:", avg_payoff_B/T
	print "x(avg):", A_avg
	print "y(avg):", B_avg
#	print "Expected payoff A by x(avg), y(avg)", Expected_payoff(M,N,payoff_A,A_avg,B_avg)
#	print "Expected payoff B by x(avg), y(avg)", Expected_payoff(M,N,payoff_B,A_avg,B_avg)
	
	check_equilibrium(M, N, payoff_A, payoff_B, A_avg, B_avg)

	print time.time() - start, "seconds went away"

	plt.figure()
	plt.plot(x_axis, y0_axis)
	plt.plot(x_axis, y1_axis)
	plt.plot(x_axis, y2_axis)

	plt.figure()
	plt.plot(x_axis, x_value)
	plt.plot(x_axis, y_value)
	grid(True)
	show()

	
def check_equilibrium(M, N, payoff_A, payoff_B, A_suspected, B_suspected):
	
	A_value = Expected_payoff(M, N, payoff_A, A_suspected, B_suspected)	
	print "Expected payoff A by x(avg), y(avg)", A_value
	B_value = Expected_payoff(M, N, payoff_B, A_suspected, B_suspected)	
	print "Expected payoff B by x(avg), y(avg)", B_value

	# check B_sus = B*?
	# if any pure strategy of A gets better payoff, then fail.
	for i in xrange(0,M):
		value = 0.0
		for j in xrange(0,N):
			value += B_suspected[j]*payoff_A[i][j]
		if (value > A_value + epsilon):
			print "warning: B_suspected may not be equilibrium", "i=", i, "value_A=", value

	# check A_sus = A*?
	# if any pure strategy of B gets better payoff, then fail.
	for j in xrange(0,N):
		value = 0.0
		for i in xrange(0,M):
			value += A_suspected[i]*payoff_B[i][j]
		if (value > B_value + epsilon):
			print "warning: A_suspected may not be equilibrium", "j=", j, "value_B=", value

def update_A(M,N,payoff,A,B):		

	for	i in xrange(0,M):				
		m_i = 0
		for j in xrange(0,N):
			m_i -= payoff[i][j]*B[j]
		
		A[i] = A[i]*(1.0-epsilon)**m_i
	A = normalize( M, A )	
	return A
	
def update_B(M,N,payoff,A,B):		

	for	j in xrange(0,N):				
		m_j = 0
		for i in xrange(0,M):
			m_j -= payoff[i][j]*A[i]
		
		B[j] = B[j]*(1.0-epsilon)**m_j
	B = normalize( N, B )	
	return B
	
	
def oracle_B(M,N,payoff,A,B):		
# generate a pure col strategy each round

	v_B = Expected_payoff(M,N,payoff,A,B)

	G = list(B)								# G as the temperal trial for B's strategy.
	Best = list(B)							# Best takes the best trial for B's strategy.

	for	t in xrange(0,N):				
		for j in xrange(0,N):
			if (j==t):
				G[j]=1.0
			else:
				G[j]=0.0
		v_G = Expected_payoff(M,N,payoff,A,G)
		if (v_G > v_B):
			Best = list(G)
			v_B = v_G
#			print Best, v_G
	return Best


def oracle_B_mixed(M,N,payoff,A,B):		
# generate a mixed strategy each round, pick one which would maximize y (B)'s payoff
	v_B = Expected_payoff(M,N,payoff,A,B)

	G = list(B)								# G as the temperal trial for B's strategy.
	Best = list(B)							# Best takes the best trial for B's strategy.

	for	t in xrange(0,1000):				# numbers of trials of finding the best strategy for B
		G = generate_random_B( N, G )
		v_G = Expected_payoff(M,N,payoff,A,G) #- distance (N,G,B)
		if (v_G > v_B):
			Best = list(G)
			v_B = v_G
			print t, Best, v_G
	return Best


def generate_random_B( N , G ):
	for k in xrange(0,N):
		G[k] = random.random()
	G = normalize (N, G)
	return G


def distance( N, G, B ):
	d = 0
	for k in xrange(0,N):
		d += (B[k]-G[k])**2
	return d


def normalize( N , V ):
	G = list(V)
	normalized_factor = 0.0
	for k in xrange(0,N):
		normalized_factor += G[k]
	for l in xrange(0,N):
		G[l] = G[l]/normalized_factor
	return G


def Expected_payoff(M,N,payoff,A,B):		
	value = 0.0
	for i in xrange(0,M):
		for j in xrange(0,N):
			value+= A[i]*B[j]*payoff[i][j]
	return value

	
def Random_matrix_generation(M,N, lower_bound, higher_bound):	
	G = []
	for i in xrange(0,M):
		row = []
		for j in xrange(0,N):
			row.append( (higher_bound-lower_bound)*random.random() + lower_bound )
		print row
		G.append(row)
	return G

	
def double_boost(M,N,payoff_A,payoff_B):	

	start = time.time()

#Uniform strategies initialization	

	A = []
	B = []
	for i in xrange(0,M):
		A.append(1.0/M)
	for j in xrange(0,N):
		B.append(1.0/N)
		
	T = int(math.log(N)/epsilon/epsilon)
	percent = int(T/1000)
	
	print "Epsilon:", epsilon
	print "Number of rounds:", T
	print "Initialized x(1): ", A
	print "Initialized y(1): ", B
	
	# for taking average
	record_A = []
	record_B = []
	
	avg_payoff_A = 0
	avg_payoff_B = 0
	
	for i in xrange(0,M):
		record_A.append(0)
	for j in xrange(0,N):
		record_B.append(0)

	x_axis = []
	y0_axis = []
	y1_axis = []
	y2_axis = []

	x_value = []
	y_value = []

	for t in xrange(0,T):
		
		B = update_B(M,N,payoff_B,A,B)	

		A = update_A(M,N,payoff_A,A,B)	
		
		for i in xrange(0,M):
			record_A[i]+=A[i]
		for j in xrange(0,N):
			record_B[j]+=B[j]
			
		avg_payoff_A += Expected_payoff(M,N,payoff_A,A,B)
		avg_payoff_B += Expected_payoff(M,N,payoff_B,A,B)
		
		if (t%percent==1):				#show something intermediate
#		if (t<10000):				#show something intermediate
	#		print "round", t, A, B
			x_axis.append(t)
			y0_axis.append(A[0])
			y1_axis.append(A[1])
			y2_axis.append(A[2])
			A_t = normalize(M, record_A)
			B_t = normalize(N, record_B)
			x_value.append(Expected_payoff(M,N,payoff_A,A_t,B_t))
			y_value.append(Expected_payoff(M,N,payoff_B,A_t,B_t))

		
	print "x(T):", A
	print "y(T):", B
	print "Expected payoff A by x(T), y(T)", Expected_payoff(M,N,payoff_A,A,B)
	print "Expected payoff B by x(T), y(T)", Expected_payoff(M,N,payoff_B,A,B)

	A_avg = normalize(M, record_A)
	B_avg = normalize(N, record_B)

	print "sum A(x(T), y(T))/T:", avg_payoff_A/T
	print "sum B(x(T), y(T))/T:", avg_payoff_B/T
	print "x(avg):", A_avg
	print "y(avg):", B_avg
#	print "Expected payoff A by x(avg), y(avg)", Expected_payoff(M,N,payoff_A,A_avg,B_avg)
#	print "Expected payoff B by x(avg), y(avg)", Expected_payoff(M,N,payoff_B,A_avg,B_avg)
	
	check_equilibrium(M, N, payoff_A, payoff_B, A_avg, B_avg)

	print time.time() - start, "seconds went away"

	plt.figure()
	plt.plot(x_axis, y0_axis)
	plt.plot(x_axis, y1_axis)
	plt.plot(x_axis, y2_axis)

	plt.figure()
	plt.plot(x_axis, x_value)
	plt.plot(x_axis, y_value)
	grid(True)
	show()

main()
