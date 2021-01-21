"""
Tested with python 27. 
"""

import math
import cmath
import random
import numpy as np

import pylab

IFPLOT = 1

pi2j = cmath.pi*2j
#print cmath.exp( 0 )
#print cmath.exp( pi2j/4 )

#prob_az_3dim(7, 7*7, 50)

# estimate by the formula.
def exp_es(p, r):
    ''' 
        This experiment is designed to check if 
            sum_x exp( 2pi/q * LWR_sample * z ) = 0 
        for all z such that Az neq 0 mod q, where q = p*r, LWR_sample = xA - [xa]%r. 
        The result shows that when Az = 0 mod p, the sum is non-zero, and there is a reason...
    '''
    c = 0
    q = p*r
    a1 = 4
    a2 = 1
    z1 = -1
    z2 = 4
    for x in range(q):
        #print x, cmath.exp( pi2j*(x*a1*z1+x*a2*z2 + ((x*a1)%r)*z1 +((x*a2)%r)*z2 )/q ) 
        c +=  cmath.exp( pi2j*(x*a1*z1+x*a2*z2 + ((x*a1)%r - (r-1)/2)*z1 +((x*a2)%r - (r-1)/2)*z2 )/q ) 
    print "c = ", c/q

#exp_es(7, 17)

def exp_es_coprime(q, r):
    ''' 
        This experiment is the same as above, 
        except that the modulus and rounding border are coprime. 
    '''
    c = 0
    a1 = 3
    a2 = 1
    z1 = 22
    z2 = 7
    for x in range(q):
        #print x, cmath.exp( pi2j*(x*a1*z1+x*a2*z2 + ((x*a1)%r)*z1 +((x*a2)%r)*z2 )/q ) 
        c +=  cmath.exp( pi2j*(x*a1*z1+x*a2*z2 + ((x*a1)%r - (r-1)/2)*z1 +((x*a2)%r - (r-1)/2)*z2 )/q ) 
    print "c = ", c/q

#exp_es_coprime(71, 7)

def Gaussian_pdf(x, mu, sigma):
    x = float(x - mu) / sigma
    return math.exp(-x*x/2.0) / math.sqrt(2.0*math.pi) / sigma

def CLT_Fourier_gaussian_1dim(s, q, T, rg):
    """ given input m = number of dimensions, s = std of gaussian, q = modulus, rg = range of |z|, 
        plot the q-fourier transform of Gaussian_s = Gaussian_{q/s} """
    mu = 0
    FT_final = [ 0 for i in range(-rg,rg)  ]
    for z in range(-rg,rg):
        Wz = 0
        for r in range(-T, T+1):
            Wz += Gaussian_pdf(r, mu, s) * cmath.exp(pi2j*r*z/q).real
        FT_final[z+rg] = Wz
        #FT_final[z+rg] = Gaussian_pdf(z, mu, s)

    print FT_final
    pylab.plot(range(-rg, rg), FT_final )
    pylab.show()

#CLT_Fourier_gaussian_1dim(4, 49, 49, 49)
#CLT_Fourier_gaussian_1dim(4, 49, 2, 49)



def sinc_square_pdf(x, B, p):
    """  sin(x/B)^p/(x/B)^p  """
    x = 2*math.pi *(x+0.000001) / B
    return (math.sin(x)/x)**p

def Fourier_sinc(B, preci, rg):
    FT_final = [ 0 for i in range(-rg,rg)  ]
    for z in range(-rg,rg):
        Wz = 0
        for r in range(-preci, preci+1):
            Wz += sinc_square_pdf(float(r), B) * cmath.exp(pi2j*float(r*z)/100).real
        FT_final[z+rg] = Wz
        #FT_final[z+rg] = sinc_square_pdf(float(z), B)

    #print FT_final
    pylab.plot(range(-rg, rg), FT_final )
    pylab.show()

def sinc_square_delta(B, p, preci):
    # using Gaussian to approximate the LPN noise in the quantum state.
    X = 0
    Y = 0   
    for i in range(-preci, preci+1):
        X+=sinc_square_pdf(2*i, B, p)
        Y+=sinc_square_pdf(i, B, p)
        if i<6 and i>-1:
            print i, sinc_square_pdf(i, B, p)
    rt =  X/Y
    print "B = ", B, "p = ", p, "X = ", X, "Y = ", Y, "rt =", rt, "d>", 1/math.log(rt,0.5)

#sinc_square_delta(4.0, 8.0, 160)


#Fourier_sinc(4.0, 10, 10)


def super_Gaussian_pdf(x, p, sigma):
    """  exp(- |x|^p )  """
    x = float(x) / sigma
    return math.exp(-(x*x)**(p/2.0))

def super_Gaussian_to_delta(s, p, preci):
    # using Gaussian to approximate the LPN noise in the quantum state.
#    X = 0
    Y = 1   
    for i in range(-preci, preci+1):
        #X+=super_Gaussian_pdf(2*i, p, s)#*sinc_square_pdf(2*i, 4.0)
        Y+=super_Gaussian_pdf(2*i+1, p, s)#*sinc_square_pdf(i, 4.0)
        if i<6 and i>-1:
            print i, super_Gaussian_pdf(i, p, s)
    rt =  1/Y
    print "s = ", s, "p = ", p, "d>", 1/math.log(rt,0.5), "T/m>", (s**p)/p, "T/n>", (s**p)/(p*math.log(rt,0.5))

#super_Gaussian_to_delta(0.9, 2.0, 40)

def Fourier_super_gaussian(s, q, p):
    """ s = width, q = modulus 
        plot the fourier transform of super_Gaussian_pdf(s, p) """
    FT_final = [ 0 for i in range(q)  ]
    for z in range(q):
        Wz = 0
        for x in range(-q, q):
            Wz += super_Gaussian_pdf(x, p, s) * cmath.exp(pi2j*float(x*z)/q).real #* sinc_square_pdf(float(r), 4.0)
        FT_final[z] = Wz
        #FT_final[z+rg] = super_Gaussian_pdf(z, p, s)

    #print FT_final
    pylab.plot(range(0, q), FT_final )
    pylab.show()

s_gl = 10
for p_gl in [0.5, 1.0, 1.5, 2.0]:
    Fourier_super_gaussian(s_gl, 40, p_gl)




def CLT_Fourier_uniform(s, q, rg):
    """ given input m = number of dimensions, s = std of gaussian, q = modulus, rg = range of |z|, 
        plot the q-fourier transform of Gaussian_s = Gaussian_{q/s} """
    mu = 0
    FT_final = [ 0 for i in range(-rg,rg)  ]
    
    # manually create a PDF of the distribution of sum of m iid numbers chosen uniformly random from [-s, s] 
    Sum_uniform_pdf = [ 0 for i in range(-rg,rg)  ]
    
    for z in range(-rg,rg):
        Wz = 0
        for r in range(-s, s+1):
            Wz += cmath.exp(pi2j*r*z/q).real / s
        FT_final[z+rg] = Wz
        #FT_final[z+rg] = Gaussian_pdf(z, mu, s)

    #print FT_final
    pylab.plot(range(-rg, rg), FT_final )
    pylab.show()

#CLT_Fourier_uniform(4, 81, 4*81)



def CLT(p, r, m, N):
    Size = p*r*4
    prod = [ 0 for i in range(2*Size+2)]
    Counter = 0
    c = 0.0
    for i in range(N):
        v = np.random.randint(-p/2+1, high= p/2+1, size=m)            
        z = np.random.randint(-r/2+1, high= r/2+1, size=m)         
        t = np.inner(v, z)
        #print v, z, t
        if t>-Size and t<Size:
            Counter+=1
            prod[t+Size] +=1   
    for i in range(-Size, Size+2):
        c+= prod[i+Size]*math.cos(math.pi*i*2/(p*r))
    print p*r, Counter, c/N

    pylab.plot(range(-Size, Size+2), prod)
    pylab.show()


def prod_uniform(p, r, m, N):
    """ old experiments about learning with rounding """
    Size = p*r*4
    prod = [ 0 for i in range(2*Size+2)]
    Counter = 0
    c = 0.0
    for i in range(N):
        v = np.random.randint(-p/2+1, high= p/2+1, size=m)            
        z = np.random.randint(-r/2+1, high= r/2+1, size=m)         
        t = np.inner(v, z)
        #print v, z, t
        if t>-Size and t<Size:
            Counter+=1
            prod[t+Size] +=1   
    for i in range(-Size, Size+2):
        c+= prod[i+Size]*math.cos(math.pi*i*2/(p*r))
    print p*r, Counter, c/N

    pylab.plot(range(-Size, Size+2), prod)
    pylab.show()


def Gaussian_MR04(x, sigma):
    # without normalization
    x = float(x) / sigma
    return math.exp(-x*x*math.pi)

def Gaussian_to_delta(s, preci):
    # using Gaussian to approximate the LPN noise in the quantum state.
    X = 0
    Y = 0   
    for i in range(-preci, preci+1):
        X+=Gaussian_MR04(2*i, s)
        Y+=Gaussian_MR04(i, s)
    rt =  X/Y
    return rt, 1/math.log(rt,0.5)

#print Gaussian_to_delta(1.4, 20)

def tau_from_Gaussian(s, preci):
    # using Gaussian to approximate the LPN noise in the quantum state.
    X = 0
    Y = 0   
    for i in range(-preci, preci+1):
        X+=Gaussian_MR04(2*i, s)
        Y+=Gaussian_MR04(i, s)
    tau = 1 - X**2/(X**2+( Y - X )**2)
    prob = 0.5+math.sqrt(tau)*math.sqrt(1-tau)
    return tau, prob, 1/math.log(prob,0.5)


def cosh_MS19(x, a):
    # p9 of MS19, 
    y = 2 * math.pi * float(a)* float(x) / math.sqrt(3)
    return 1/(1+2*math.cosh(y))


def delta_from_cosh(s, preci):
    # using cosh_ms19 to approximate the probability directly.
    X = 0
    Y = 0   
    for i in range(-preci, preci+1):
        X+=cosh_MS19(2*i, 1/s)
        Y+=cosh_MS19(i, 1/s)
    rt =  X/Y
    print s, "d>", 1/math.log(rt,0.5), "T/m>", s*0.3927633, "T/n>", s*0.3927633/math.log(rt,0.5)


def tau_from_cosh(a, preci):
    # using cosh_ms19 to approximate the LPN noise in the quantum state.
    X = 0
    Y = 0   
    for i in range(-preci, preci+1):
        X+=cosh_MS19(2*i, a)
        Y+=cosh_MS19(i, a)
    tau = 1 - X**2/(X**2+( Y - X )**2)
    prob = 0.5+math.sqrt(tau)*math.sqrt(1-tau)
    return tau, prob, 1/math.log(prob,0.5), math.log(prob,0.5)/0.3928

#print tau_from_cosh(2, 10)


def sech_MS19(x, a):
    # p9 of MS19, 
    y =  math.pi * float(a)* float(x) 
    return 1/math.cosh(y)

def tau_from_sech(a, preci):
    # using sech to approximate the LPN noise in the quantum state.
    X = 0
    Y = 0   
    for i in range(-preci, preci+1):
        X+=sech_MS19(2*i, a)
        Y+=sech_MS19(i, a)
    tau = 1 - X**2/(X**2+( Y - X )**2)
    prob = 0.5+math.sqrt(tau)*math.sqrt(1-tau)
    return tau, prob, 1/math.log(prob,0.5)

def H2(delta):
    return -delta*math.log(delta, 2) - (1-delta)*math.log( (1-delta), 2)

def tau_c_d(tau, c):
    # input tau to decide the applicable d and c for m = dn, T = n/c
    prob = 0.5+math.sqrt(tau)*math.sqrt(1-tau)
    d = 1/math.log(prob,0.5)
    prob_m = ((1-prob)/prob)**(1/c)
    prob_l = 2**(d*H2(1/(c*d)))
    # want to see (prob_l*prob_m*0.5)^n

    return tau, prob, "d = ", d, prob_l*prob_m/2

#print tau_c_d(0.1, 2.0)
#for mu in range(30):
#prod_uniform(17, 13, 100, 100000)



'''
# using the formula
    for x in range(-(r/2), r/2+1):
        e = ((x*a)%r)
        for i in range(m):
            if e[i]>r/2:
                e[i] = e[i]-r
        t = np.inner(e, z)
        #print e, t
        mu+= t
        ez+= cmath.exp( pi2j *t/(p*r) )

        if IFPLOT:
            prob_p[(t/r)%p] +=1
#            prob_p[ t%(p*r*m/4) ] +=1
'''

