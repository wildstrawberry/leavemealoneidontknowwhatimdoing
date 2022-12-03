"""
Tested with python 2.7. 
"""

import math
import cmath
import random
import numpy as np
import scipy.special as sp
import pylab
import matplotlib.pyplot as plt

#plt.matplotlib.rc('text', usetex = True)
#plt.matplotlib.rc('grid', linestyle = 'dotted')
#plt.matplotlib.rc('figure', figsize = (6.4,4.8)) # (width,height) inches
# G[col][row]

pi2j = cmath.pi*2j

def convol(a, b, q):
    #convolution of two vectors of dimension q
    c = np.array( [ 0.0 for i in range(q)] )
    for i in range(q):
        for j in range(q):
            c[i]+= a[(i-j)%q]*b[j]
    return c

def FT(f, q):
    #input: f, output, hat f over Z_q
    FT_array = np.array( [ 0.0 for i in range(q)] )
    for z in range(q):
        Wz = 0
        for x in range(q):
            Wz += f[x] * cmath.exp(pi2j*float(x*z)/q).real / math.sqrt(q) 
            # since we only take the real part, our FT only works for symmetric functions
        FT_array[z] = Wz
    return FT_array

def Gauss(x2, sigma):
    """  exp(- |x|^2 ), it seems that we only need p = 2.0  """
    x2 = float(x2) / (sigma*sigma)
    return math.exp(-3.1415*x2)


def FTsquare(f, q):
    #input: f, output, (hat f)^2 over Z_q
    FT_array = np.array( [ 0.0 for i in range(q)] )
    for z in range(q):
        Wz = 0
        for x in range(q):
            Wz += f[x] * cmath.exp(pi2j*float(x*z)/q).real / math.sqrt(q) 
            # since we only take the real part, our FT only works for symmetric functions
        FT_array[z] = Wz**2
    return FT_array

def FT2d(f, q, square_bool):
    #input: f, output, hat f
    #if square_bool =2 then return the square of FT
    FT_array = np.array( [[ 0.0 for i in range(q)] for j in range(q) ] )
    for y in range(q):
      for z in range(q):
        W_yz = 0
        for i in range(q):
          for j in range(q):
            W_yz += f[j][i] * cmath.exp(pi2j*float(i*y+j*z)/q).real 
            # since we only take the real part, our FT only works for symmetric functions
        if square_bool ==1:
            FT_array[z][y] = W_yz
        elif square_bool ==2:
            FT_array[z][y] = W_yz**2
    return FT_array

def FT3d(f, q):
    #input: f, output, hat f
    FT_array = np.array( [[[ 0.0 for i in range(q)] for j in range(q) ] for k in range(q) ] )
    for z in range(q):
      for y in range(q):
        for x in range(q):
            W_xyz = 0
            for k in range(q):
              for j in range(q):
                for i in range(q):
                  W_xyz += f[k][j][i] * cmath.exp(pi2j*float(i*x+j*y+k*z)/q).real 
            # since we only take the real part, our FT only works for symmetric functions
            FT_array[z][y][x] = W_xyz
    return FT_array


def FT_sphere(radius, thickness, V, q):

    extentg = [-(q)/2, (q)/2, -(q)/2, (q)/2]

    I_sphere = np.array( [[[ 0.0 for i in range(q)] for j in range(q) ] for k in range(q) ] )
    for k in range(-(q)/2, (q)/2):
      for j in range(-(q)/2, (q)/2):
        for i in range(-(q)/2, (q)/2):
            distance = math.sqrt(i**2+j**2+k**2)
            theta = 0
            if distance != 0:
                theta = np.arccos((i+0.0000)/distance)
                if j<0:
                    theta=-theta
            #print theta
            f_R_V = Gauss( (distance - radius)**2, V)
            hat_h_T = sp.jv(0, 2*3.1415*distance/thickness)
            I_sphere[k][j][i] = f_R_V #* hat_h_T * np.sin(16*theta)#*np.sign(i)
            #I_sphere[j][i]= I_sphere[j][i]**2

    I_sphere_show = np.array( [[ 0.0 for i in range(q)] for j in range(q) ] )
    for i in range(-(q)/2, (q)/2):
        for j in range(-(q)/2, (q)/2):
            I_sphere_show[j+(q)/2][i+(q)/2] = I_sphere[0][j][i]

    plt.imshow(I_sphere_show, extent = extentg, interpolation='none', cmap='RdBu')
    plt.colorbar()
    plt.show()
    plt.close()

    G_hat = FT3d(I_sphere, q)

    G_hat_show = np.array( [[ 0.0 for i in range(q)] for j in range(q) ] )
    for i in range(-(q)/2, (q)/2):
        for j in range(-(q)/2, (q)/2):
            G_hat_show[j+(q)/2][i+(q)/2] = G_hat[0][j][i]


    plt.imshow(G_hat_show, extent = extentg, interpolation='none', cmap='RdBu')
    plt.colorbar()
    plt.show()
    plt.close()

    
for i in range(4,5):

    FT_sphere(i, 4.0, 8.0, 10)
