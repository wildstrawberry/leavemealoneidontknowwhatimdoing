# stable

import math
import cmath
import numpy as np
import scipy.special as sp
import matplotlib.pyplot as plt
plt.matplotlib.rc('text', usetex = True)
plt.matplotlib.rc('grid', linestyle = 'dotted')
plt.matplotlib.rc('figure', figsize = (6.4,4.8)) # (width,height) inches

# scipy.special.sph_harm(m, n, theta, phi, out=None)
#|m| <= n,  n >= 0, theta in [0, 2*pi], phi in [0, pi].

pi = math.pi

def SHsum(ell, theta, phi, coeff):
	# compute the sum of spherical harmonics of order ell
	# sph_harm(m, n, theta, phi, out=None)  theta \in [0, 2*pi], phi \in [0, pi]
	s = sp.sph_harm(0, ell, theta, phi)
	for m in range(1, ell+1):
		s += sp.sph_harm(m, ell, theta, phi) * coeff[m]   #* cmath.exp(2*pi *2*m* 1j /(ell+1))
		s += sp.sph_harm(-m, ell, theta, phi) * coeff[-m] #* 1j#cmath.exp(2*pi *(2*m)* 1j /(ell))
	return s

def plot_SH(ell):
	slices = 90
	range_phi = pi
	range_theta = 2*pi
	phi_list = np.linspace(0, pi, slices+1)
	theta_list = np.linspace(0, 2*pi, 2*slices+1)

	coeffreal_best = []
	coeff_best = []
	maxx_best = ell*100

	for trials in range(1):
		coeffreal = np.random.rand(2*ell+1)
		coeff = np.array( [ cmath.exp(2*pi * 1j * coeffreal[i]) for i in range(2*ell+1)  ] )
		#print coeffreal
		#print coeff

		y_list = np.array([[ math.sqrt(4*pi/(2*ell+1))*np.abs(SHsum(ell, theta, phi, coeff)) \
								for theta in theta_list] for phi in phi_list] )
		maxx = np.amax(y_list)
	
		y_sum = sum([ np.linalg.norm(y_list[i])**2 * np.sin(pi* i /slices) for i in range(slices+1)  ])
	
		print trials, ell, maxx, np.amin(y_list), y_sum/((slices+1)*(2*slices+1))
	
		if maxx<maxx_best:
			coeffreal_best = coeffreal
			coeff_best = coeff
			maxx_best = maxx

	print coeffreal_best
	print maxx_best
	y_list = np.array([[ math.sqrt(4*pi/(2*ell+1))*np.abs(SHsum(ell, theta, phi, coeff_best)) \
								for theta in theta_list] for phi in phi_list] )
	plt.imshow(y_list, interpolation='none', cmap='Blues')
	plt.colorbar()
	plt.show()
	plt.close()


plot_SH(3)

def sample_SH(ell):
	''' test SH on several points '''
	slices = 50
	range_phi = pi
	range_theta = 2*pi
	phi_list = np.linspace(0, pi, slices+1)
	theta_list = np.linspace(0, 2*pi, 2*slices+1)
	print phi_list
	print theta_list

	coeffreal_best = []
	coeff_best = []
	maxx_best = ell*100

	y_list_sum = np.array([[ 0.0 for theta in theta_list] for phi in phi_list] )

	for trials in range(100):
		coeffreal = np.random.rand(2*ell+1)
		coeff = np.array( [ cmath.exp(2*pi * 1j * coeffreal[i]) for i in range(2*ell+1)  ] )
		#print coeffreal
		#print coeff

		y_list = np.array([[ math.sqrt(4*pi/(2*ell+1))*np.abs(SHsum(ell, theta, phi, coeff)) \
								for theta in theta_list] for phi in phi_list] )
		maxx = np.amax(y_list)

		y_list_sum+=y_list
	
		y_sum = sum([ np.linalg.norm(y_list[i])**2 * np.sin(pi* i /slices) for i in range(slices+1)  ])
	
		print trials, ell, "max:", maxx, "min:", np.amin(y_list), "norm:", y_sum/((slices+1)*(2*slices+1))
	
		if maxx<maxx_best:
			coeffreal_best = coeffreal
			coeff_best = coeff
			maxx_best = maxx

	#print coeffreal_best
	print "best score:", maxx_best
	#y_list = np.array([[ math.sqrt(4*pi/(2*ell+1))*np.abs(SHsum(ell, theta, phi, coeff_best)) \
	#							for theta in theta_list] for phi in phi_list] )]
	print y_list_sum[0][slices/2], y_list_sum[1][slices/2], y_list_sum[slices/2][slices/2], y_list_sum[slices/2][slices], y_list_sum[slices/2][3*slices/2], y_list_sum[slices/2+1][3*slices/2], y_list_sum[slices/2-1][3*slices/2]
	plt.imshow(y_list_sum, interpolation='none', cmap='Blues')
	plt.colorbar()
	plt.show()
	plt.close()

#sample_SH(17)


