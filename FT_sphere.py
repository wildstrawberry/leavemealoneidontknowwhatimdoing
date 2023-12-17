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

def sinc( B, q):
    """  QFT of bounded uniform of width 2B+1, exp(- |x|^2 ), it seems that we only need p = 2.0  """
    BU = np.array( [ 0.0 for i in range(q)] )
    for i in range(-B, B+1):
        BU[i] = 1.0
    sincf = FT(BU, q)

    '''print("sinc with sigma, q", B, q)
    plt.scatter( range(0, q), sinc )
    plt.show()
    plt.close()
    '''
    return sincf

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
    #input: f, output, (hat f)^2 over Z^2_q
    #if square_bool =2 then return the square of FT
    FT_array = np.array( [[ 0.0 for i in range(q)] for j in range(q) ] )
    for y in range(q):
      for z in range(q):
        W_yz = 0
        for i in range(q):
          for j in range(q):
            W_yz += f[i][j] * cmath.exp(pi2j*float(i*y+j*z)/q).real 
            # since we only take the real part, our FT only works for symmetric functions
        if square_bool ==1:
            FT_array[y][z] = W_yz
        elif square_bool ==2:
            FT_array[y][z] = W_yz**2
    return FT_array

def FT_window(r1, r2, q):

    sinc_funcr1 = sinc(r1, q)
    sinc_funcr2 = sinc(r2, q)

#    I_window = np.array( [[ 0.0 for i in range(q)] for j in range(q) ] )
#    for i in range((q)/2-r1, (q+1)/2+r1+1):
#        for j in range(-r2, r2+1):
#            I_window[i][j] = 1.0

    I_window = np.array( [[ Gauss( (i - float(q)/2)**2, r1)*(Gauss(j**2, r2) + Gauss((j-q)**2, r2)) for i in range(q)] for j in range(q) ] )
#    I_window = np.array( [[ sinc_funcr1[i - (q)/2]*sinc_funcr2[j] for i in range(q)] for j in range(q) ] )

    plt.imshow(I_window, interpolation='none')
    plt.colorbar()
    plt.show()
    plt.close()

    G_hat = FT2d(I_window, q, 1)
    plt.imshow(G_hat, interpolation='none')
    plt.colorbar()
    plt.show()
    plt.close()


def standardshow(I_sphere, q, balanced):
    """ if balanced=1, add a max point and min point so yhat  """
    extentg = [-(q)/2, (q)/2, -(q)/2, (q)/2]
    I_sphere_show = np.array( [[ 0.0 for i in range(q)] for j in range(q) ] )
    v_max = 0
    v_min = 0
    for i in range(-(q)/2, (q)/2):
        for j in range(-(q)/2, (q)/2):
            I_sphere_show[j+(q)/2][i+(q)/2] = I_sphere[j][i]
            if I_sphere[j][i]>v_max:
                v_max = I_sphere[j][i]
            if I_sphere[j][i]<v_min:
                v_min = I_sphere[j][i]

    if balanced == 1:
        #print v_max, v_min
        I_sphere_show[0][0] = v_max
        I_sphere_show[0][1] = v_min
        I_sphere_show[1][0] = -v_max
        I_sphere_show[1][1] = -v_min
        plt.imshow(I_sphere_show, extent = extentg, interpolation='none', cmap='RdBu')
        plt.colorbar()
        plt.show()
        plt.close()
    else:
        plt.imshow(I_sphere_show, extent = extentg, interpolation='none', cmap='Blues') #cmap='RdBu')
        plt.colorbar()
        plt.show()
        plt.close()

def FT_sphere(radius, thickness, V, r1, r2, q):

    extentg = [-(q)/2, (q)/2, -(q)/2, (q)/2]

    I_sphere = np.array( [[ 0.0 for i in range(q)] for j in range(q) ] )
    for i in range(-(q)/2, (q)/2):
        for j in range(-(q)/2, (q)/2):
            distance = math.sqrt(i**2+j**2)
            f_R_V = Gauss( (distance - radius)**2, V)
            hat_h_T = sp.jv(0, 2*3.1415*distance/thickness)
            I_sphere[i][j] = f_R_V*hat_h_T 
            # I_sphere[i][j] = Gauss( abs(i**2+j**2 - radius**2) , thickness)
            #if i**2+j**2 <= (radius+thickness)**2: # and i**2+j**2 >= (radius-thickness)**2:
            #if i**2+j**2 <= (radius+thickness)**2 and i**2+j**2 >= (radius-thickness)**2:
            #if i**2+j**2 <= (radius+1)**2 and i**2+j**2 >= (radius-1)**2:
            #    I_sphere[i][j] = 1.0
            #I_sphere[i][j] = Gauss( (i)**2, 5.0) + Gauss( (i-q/5)**2, 5.0)+ Gauss( (i+q/5)**2, 5.0)+ Gauss( (i-2*q/5)**2, 5.0)+ Gauss( (i+2*q/5)**2, 5.0)
            #I_sphere[i][j] *=  Gauss( (j+q/3)**2, 5.0)

    standardshow(I_sphere, q, 1)

    G_hat = FT2d(I_sphere, q, 1)

    standardshow(G_hat, q, 1)

    # put a rectangular window on G_hat
    '''
    W_G_hat = np.array( [[ 0.0 for i in range(q)] for j in range(q) ] )
    for i in range((q)/2-r1, (q+1)/2+r1+1):
      for j in range(-r2, r2+1):
        W_G_hat[i][j] = G_hat[i][j]
    for i in range(12, 12):
      for j in range(-2, 3):
        W_G_hat[i][j] = G_hat[i][j]
    '''
    
    
    #window on the x axis
    W_G_hat = np.array( [[ G_hat[i][j]*Gauss( (i - float(q)/6)**2, r1)*(Gauss(j**2, r2) + Gauss((j-q)**2, r2)) for i in range(q)] for j in range(q) ] )
    #window in the middle
    #W_G_hat = np.array( [[ G_hat[i][j] * Gauss( (i - float(q)/5)**2, r1) * Gauss( (j- float(q)/4)**2, r2 )  for i in range(q)] for j in range(q) ] )

    #
    #sinc_funcr1 = sinc(r1, q)
    #sinc_funcr2 = sinc(r2, q)
    #W_G_hat = np.array( [[ G_hat[i][j]*sinc_funcr1[i - (q)/2]*sinc_funcr2[j] for i in range(q)] for j in range(q) ] )
    standardshow(W_G_hat, q, 1)


    W_G = FT2d(W_G_hat, q, 1)

    #prepare W_G stats 
    W_G_normal = 0.0 # the normalization factor
#    critical_area = 0.0

    for y in range(q):
      for z in range(q):
        W_G_normal += W_G[y][z]
#        if (y<radius/2 or y > q - radius/2) and (z>2*radius/3 and z<q - 2*radius/3):
#            critical_area += W_G[y][z]

    print("radius/2", radius/2, "2radius/3", 2*radius/3)
#    print("normalization factor", W_G_normal)
#    print("critical area:", critical_area/W_G_normal)

    #W_G = W_G/W_G_normal

    standardshow(W_G, q, 1)

    var_y = 0.0
    var_z = 0.0
    for y in range(-q/2,q/2+1):
      for z in range(-q/2,q/2+1):
        #if y==0:
        #    print y, z, W_G[y][z]
        var_y += W_G[y][z]*y*y
        var_z += W_G[y][z]*z*z
    print("var_y, |||", var_y)
    print("var_z, ---", var_z)
    

    
for i in range(17,18):

    #FT_window(3, 3, 48)
    FT_sphere(i, 5.0, 10.0, 10, 3, 50)
    #FT_sphere(i, 2.0, 10.0, 6, 6, 40)
    #FT_sphere_another(i, 2.0, 10.0, 10, 6, 80)
    #cosine_sphere_intersection(i, 5.0, 20.0, 400)
    


def cosine_sphere_intersection(radius, thickness, V, q):
    # g_R,V,T = f_{R, V}\cdot hat(h_V)

    I_sphere = np.array( [[ 0.0 for i in range(q)] for j in range(q) ] )
    for i in range(-(q)/2, (q)/2):
        for j in range(-(q)/2, (q)/2):
            distance = math.sqrt(i**2+j**2)
            f_R_V = Gauss( (distance - radius)**2, V)
            hat_h_T = sp.jv(0, 2*3.1415*distance/thickness)
            I_sphere[i][j] = f_R_V*hat_h_T 


    plt.imshow(I_sphere, interpolation='none', cmap='Blues')
    plt.colorbar()
    plt.show()
    plt.close()

    range_x = 120
    xlist = np.linspace(0, range_x, range_x)
    ylist = []
    for d in xlist:
        sums = 0.0
        for i in range(q):
            for j in range(q):
                sums += I_sphere[i][j]*I_sphere[i-d][j]
        ylist.append(sums)

    plt.plot(xlist, ylist )
    plt.xlim((0, range_x))
    plt.ylim((-10, 24))
    plt.legend((  'intersection' ), loc = 0)
    plt.xlabel('$x$')
    plt.ylabel('intersection')
    plt.grid(True)
    plt.tight_layout(0.5)
    plt.savefig('intersection.pdf')



def FT_sphere_another(radius, thickness, V, r1, r2, q):
    # the purpose is to investigate different shapes of thin spheres, and their fourier transform

    I_sphere = np.array( [[ 0.0 for i in range(q)] for j in range(q) ] )
    for i in range(-(q)/2, (q)/2):
        for j in range(-(q)/2, (q)/2):
            # I_sphere[i][j] = Gauss( abs(i**2+j**2 - radius**2) , thickness)
            #if i**2+j**2 <= (radius+thickness)**2: # and i**2+j**2 >= (radius-thickness)**2:
            #if i**2+j**2 <= (radius+thickness)**2 and i**2+j**2 >= (radius-thickness)**2:
            #if i**2+j**2 <= (radius+1)**2 and i**2+j**2 >= (radius-1)**2:
            #    I_sphere[i][j] = 1.0
            distance = math.sqrt(i**2+j**2)
            f_R_V = Gauss( (distance - radius)**2, V)
            hat_h_T = sp.jv(8, 2*3.1415*distance/thickness)
            I_sphere[i][j] = f_R_V*hat_h_T 

    plt.imshow(I_sphere, interpolation='none', cmap='Blues')
    plt.colorbar()
    plt.show()
    plt.close()

    G_hat = FT2d(I_sphere, q, 1)
    plt.imshow(G_hat, interpolation='none', cmap='Blues')
    plt.colorbar()
    plt.show()
    plt.close()

    # put a Gaussian filter to filter out the high frequency domain
    for i in range(-(q)/2, (q)/2):
        for j in range(-(q)/2, (q)/2):
            distance = math.sqrt(i**2+j**2)
            #G_hat[i][j] *= 2*Gauss( abs(distance**2), q/radius) + Gauss( (distance-r1)**2, q/radius)  #*(distance**8)
            #G_hat[i][j] *= Gauss( abs(distance**2), q/r1)*(distance**2)
            G_hat[i][j] *= Gauss( (distance-r1)**2, r2) #/ ((distance+0.01)**0.5)
            #G_hat[i][j] = G_hat[i][j]*G_hat[i][j]

    plt.imshow(G_hat, interpolation='none', cmap='Blues')
    plt.colorbar()
    plt.show()
    plt.close()

    G_convol = FT2d(G_hat, q, 1)
    plt.imshow(G_convol, interpolation='none', cmap='Blues')
    plt.colorbar()
    plt.show()
    plt.close()
