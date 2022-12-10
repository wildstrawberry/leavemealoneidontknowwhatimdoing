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

def flip(p):
    return 1.0 if random.random() < p else 0.0

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


def leg(i, q):
    # compute the legendre symbol, assuming q is prime
    x = (i**((q-1)/2))%q
    if x == q-1:
        x = -1
    return x

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


def FT3d(f, q, square_bool):
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
            if square_bool ==1:
                FT_array[z][y][x] = W_xyz
            elif square_bool ==2:
                FT_array[z][y][x] = W_xyz**2
    return FT_array


def overlap(in1, in2, q):
    summm=0.0
    for i in range(q):
        for j in range(q):
            summm += in1[j][i]*in2[j][i]
    return summm


def FT_sphere(radius, q, bool_A):

    extentg = [-(q)/2, (q)/2, -(q)/2, (q)/2]

    I_sphere = np.array( [[ 0.0 for i in range(q)] for j in range(q) ] )
    for j in range(q):
        for i in range(q):
            distance = (i**2+j**2)%q
            #distance2 = ((i+1)**2+j**2)%q
            if distance == radius:
                I_sphere[j][i] = 1.0
            #if distance2 == radius:
            #    I_sphere[j][i]*= -1.0

    I_sphere_show = np.array( [[ 0.0 for i in range(q)] for j in range(q) ] )
    for i in range(-(q)/2, (q)/2):
        for j in range(-(q)/2, (q)/2):
            I_sphere_show[j+(q)/2][i+(q)/2] = I_sphere[j][i]

    plt.imshow(I_sphere_show, extent = extentg, interpolation='none', cmap='Blues') #cmap='RdBu')
    plt.colorbar()
    plt.show()
    plt.close()

    G_hat = FT2d(I_sphere, q, 2)
    if bool_A ==1:
        a_1 = 1
        a_2 = 4
        for i in range(q):
            for j in range(q):
                if (a_1 *j+a_2*i)%q != 0:
                    G_hat[j][i] = 0

    # add the restriction that a_1 y_1 + a_2 y_2 = 0 mod q 

    G_hat_show = np.array( [[ 0.0 for i in range(q)] for j in range(q) ] )
    for i in range(-(q)/2, (q)/2):
        for j in range(-(q)/2, (q)/2):
            G_hat_show[j+(q)/2][i+(q)/2] = G_hat[j][i]
            if i==0 and j==0:
                G_hat_show[j+(q)/2][i+(q)/2]=0#*dis
    #print(G_hat[0][0]) #, G_hat[1][0], G_hat[0][1], G_hat[-1][0], G_hat[0][-1])

    plt.imshow(G_hat_show, extent = extentg, interpolation='none', cmap='RdBu')
    plt.colorbar()
    plt.show()
    plt.close()

    G_hat_hat = FT2d(G_hat, q, 1)

    G_hat_hat_show = np.array( [[ 0.0 for i in range(q)] for j in range(q) ] )
    for i in range(-(q)/2, (q)/2):
        for j in range(-(q)/2, (q)/2):
            G_hat_hat_show[j+(q)/2][i+(q)/2] = G_hat_hat[j][i]
            if i==0 and j==0:
                G_hat_hat_show[j+(q)/2][i+(q)/2]=0#*dis
    #print(G_hat[0][0]) #, G_hat[1][0], G_hat[0][1], G_hat[-1][0], G_hat[0][-1])

    plt.imshow(G_hat_hat_show, extent = extentg, interpolation='none', cmap='Blues')
    plt.colorbar()
    plt.show()
    plt.close()

#FT_sphere(1, 19, 0)
#FT_sphere(1, 19, 1)



def FT_sphere_3d(radius, q, bool_A):
    print(" radius: ", radius, q, "mode:", bool_A )

    extentg = [-(q)/2, (q)/2, -(q)/2, (q)/2]

    I_sphere = np.array( [[[ 0.0 for i in range(q)] for j in range(q) ] for k in range(q) ] )
    for k in range(q):
      for j in range(q):
        for i in range(q):
            distance = (i**2+j**2+k**2)%q
            #distance2 = ((i+1)**2+j**2)%q
            if distance in radius:
                I_sphere[k][j][i] = 1.0
            #if distance2 == radius:
            #    I_sphere[j][i]*= -1.0

    I_sphere_show = np.array( [[ 0.0 for i in range(q)] for j in range(q) ] )
    for i in range(-(q)/2, (q)/2):
        for j in range(-(q)/2, (q)/2):
            I_sphere_show[j+(q)/2][i+(q)/2] = I_sphere[0][j][i]

    plt.imshow(I_sphere_show, extent = extentg, interpolation='none', cmap='Blues') #cmap='RdBu')
    plt.colorbar()
    plt.show()
    plt.close()

    I_sphere_show2 = np.array( [[ 0.0 for i in range(q*(q/2))] for j in range(q) ] )
    for k in range((q/2)):
      for j in range(q):
        for i in range(q):
            I_sphere_show2[j][k*q+i] = I_sphere[k][j][i]

    plt.imshow(I_sphere_show2, interpolation='none', cmap='Blues') #cmap='RdBu')
    plt.colorbar()
    plt.show()
    plt.close()

    G_hat = FT3d(I_sphere, q, 2)
    if bool_A ==1:
        a_1 = 1
        a_2 = 1
        a_3 = 2
        b_1 = 0
        b_2 = 1
        b_3 = 3
        for k in range(q):
            for j in range(q):
                for i in range(q):
                    if (a_1*k+a_2*j+a_3*i)%q != 0: #or (b_1*k+b_2*j+b_3*i)%q != 0:
                        G_hat[k][j][i] = 0

    # add the restriction that a_1 y_1 + a_2 y_2 = 0 mod q 
    '''
    G_hat_show = np.array( [[ 0.0 for i in range(q)] for j in range(q) ] )
    for i in range(-(q)/2, (q)/2):
        for j in range(-(q)/2, (q)/2):
            G_hat_show[j+(q)/2][i+(q)/2] = G_hat[0][j][i]
            #if i==0 and j==0:
            #    G_hat_show[j+(q)/2][i+(q)/2]=0#*dis

    plt.imshow(G_hat_show, extent = extentg, interpolation='none', cmap='RdBu')
    plt.colorbar()
    plt.show()
    plt.close()

    G_hat_show2 = np.array( [[ 0.0 for i in range(q*q)] for j in range(q) ] )
    for k in range(q):
      for j in range(q):
        for i in range(q):
            G_hat_show2[j][k*q+i] = G_hat[k][j][i]


    plt.imshow(G_hat_show2, interpolation='none', cmap='RdBu')
    plt.colorbar()
    plt.show()
    plt.close()
    '''
    stat = [ 0 for i in range(q) ]
    stat_weight = [ 0 for i in range(q) ]
    stat_r = [ 0 for i in range(q*q) ]
    stat_r_weight = [ 0 for i in range(q*q) ]


    G_hat_hat = FT3d(G_hat, q, 1)


    G_hat_hat_show2 = np.array( [[ 0.0 for i in range(q*q)] for j in range(q) ] )
    for k in range(-(q)/2, (q)/2):
      for j in range(-(q)/2, (q)/2):
        for i in range(-(q)/2, (q)/2):
            distance = (i**2+j**2+k**2)
            stat_r[distance]+=1
            stat_r_weight[distance]+=G_hat_hat[k][j][i]
            distance = distance%q
            stat[distance]+=1
            stat_weight[distance]+=G_hat_hat[k][j][i]

            G_hat_hat_show2[j+(q)/2][k*q+i+(q)/2] = G_hat_hat[k][j][i]
            if k==0 and j==0 and i==0:
                G_hat_hat_show2[j][k*q+i]=0
    #print stat_r
    #print stat_r_weight
    #print stat
    #print stat_weight

    print "(without mod q) distance avg"
    for i in range(q):
        if stat_r[i]>0:
            print i, stat_r_weight[i]/stat_r[i]

    print "mod q distance avg"
    for i in range(q):
        if stat[i]>0:
            print i, stat_weight[i]/stat[i]

    G_hat_hat_show = np.array( [[ 0.0 for i in range(q)] for j in range(q) ] )
    for i in range(-(q)/2, (q)/2):
        for j in range(-(q)/2, (q)/2):
            G_hat_hat_show[j+(q)/2][i+(q)/2] = G_hat_hat[0][j][i]
            if i==0 and j==0:
                G_hat_hat_show[j+(q)/2][i+(q)/2]=0#*dis
    #print(G_hat[0][0]) #, G_hat[1][0], G_hat[0][1], G_hat[-1][0], G_hat[0][-1])

    plt.imshow(G_hat_hat_show, extent = extentg, interpolation='none', cmap='Blues')
    plt.colorbar()
    plt.show()
    plt.close()

    plt.imshow(G_hat_hat_show2, interpolation='none', cmap='Blues') #cmap='RdBu')
    plt.colorbar()
    plt.show()
    plt.close()


#FT_sphere_3d([1,5,4], 29, 0)
FT_sphere_3d([4], 16, 1)





def FT_legendre_sphere(q):

    # mode = 0: assign amplitude; mode = 1: sample accordinh to the amplitude

    extentg = [-(q)/2, (q)/2, -(q)/2, (q)/2]

    #print -(q)/2, (q)/2

    I_sphere = np.array( [[ 0.0 for i in range(q)] for j in range(q) ] )
    for j in range(q):
        for i in range(q):
            distance = (i**2+j**2)%q
            #distance2 = ((i+1)**2+j**2)%q
            if leg(distance, q) == 1:
                I_sphere[j][i] = 1.0
            #if leg(distance, q) ==0:
                #I_sphere[j][i] = 1.0
            #elif leg(distance, q) == -1:
            #    I_sphere[j][i] = -1.0

    I_sphere_show = np.array( [[ 0.0 for i in range(q)] for j in range(q) ] )
    for i in range(-(q)/2, (q)/2):
        for j in range(-(q)/2, (q)/2):
            I_sphere_show[j+(q)/2][i+(q)/2] = I_sphere[j][i]

    plt.imshow(I_sphere_show, extent = extentg, interpolation='none', cmap='Blues') #cmap='RdBu')
    plt.colorbar()
    plt.show()
    plt.close()

    G_hat = FT2d(I_sphere, q, 1)


    '''a_1 = 1
    a_2 = 6
    for i in range(q):
        for j in range(q):
            if (a_1 *j+a_2*i)%q != 0:
                G_hat[j][i] = 0
    '''


    #G_hat = np.array( [[ 0.0 for i in range(q)] for j in range(q) ] )

    G_hat_show = np.array( [[ 0.0 for i in range(q)] for j in range(q) ] )
    for i in range(-(q)/2, (q)/2):
        for j in range(-(q)/2, (q)/2):
            G_hat_show[j+(q)/2][i+(q)/2] = G_hat[j][i]
            if i==0 and j==0:
                G_hat_show[j+(q)/2][i+(q)/2]=0#*dis
    #print(G_hat[0][0]) #, G_hat[1][0], G_hat[0][1], G_hat[-1][0], G_hat[0][-1])

    plt.imshow(G_hat_show, extent = extentg, interpolation='none', cmap='Blues')
    plt.colorbar()
    plt.show()
    plt.close()

    G_back = FT2d(G_hat, q, 1)

    #G_hat = np.array( [[ 0.0 for i in range(q)] for j in range(q) ] )

    G_back_show = np.array( [[ 0.0 for i in range(q)] for j in range(q) ] )
    for i in range(-(q)/2, (q)/2):
        for j in range(-(q)/2, (q)/2):
            G_back_show[j+(q)/2][i+(q)/2] = G_back[j][i]
            if i==0 and j==0:
                G_back_show[j+(q)/2][i+(q)/2]=0#*dis
    #print(G_hat[0][0]) #, G_hat[1][0], G_hat[0][1], G_hat[-1][0], G_hat[0][-1])

    plt.imshow(G_back_show, extent = extentg, interpolation='none', cmap='RdBu')
    plt.colorbar()
    plt.show()
    plt.close()


#FT_legendre_sphere(17)



