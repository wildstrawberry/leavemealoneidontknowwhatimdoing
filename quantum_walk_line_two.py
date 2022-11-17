"""
Two walker quantum walk. 
Modified from https://gist.github.com/qubyte/3768908
"""

import math
import cmath
import random
import numpy as np
import scipy.special as sp
import pylab
import matplotlib.pyplot as plt


STEPS = 100

pi2j = cmath.pi*2j
def FT(f, q):
    #input: f, output, hat f over Z_q
    FT_array = np.array( [ 1j for i in range(q)] )
    for z in range(q):
        Wz = 0*1j
        for x in range(q):
            Wz += f[x] * cmath.exp(pi2j*float(x*z)/q) / math.sqrt(q) 
        FT_array[z] = Wz
    return FT_array

def Gauss(x2, sigma):
    """  exp(- |x|^2 ), it seems that we only need p = 2.0  """
    x2 = float(x2) / (sigma*sigma)
    return math.exp(-3.1415*x2)


def probabilities(posn):
    """Returns a list of the probabilies for each place."""
    return [sum([abs(amp) ** 2 for amp in place]) for place in posn]

def normalise(posn):
    """Normalise function to normalise an input 1D line."""
    N = math.sqrt(sum(probabilities(posn)))
    return [[amp / N for amp in place] for place in posn]


def Hadamard_Coin(posn):
    """At each timestep, apply a Hadamard gate on the coin state."""
    return normalise( [[x[0] + x[1], x[0] - x[1]] for x in posn] )

def General_Coin(posn):
    """A general 1-qubit gate on each element parameterized by rho, phi, and theta."""
    rho = 0.3
    phi = 0.3
    theta = 0.3
    return normalise( [ [ math.sqrt(rho)*x[0] + math.sqrt(1 - rho)*cmath.exp(cmath.pi*1j*phi)*x[1], math.sqrt(1 - rho)*cmath.exp(cmath.pi*1j*theta)*x[0] - math.sqrt(rho)*cmath.exp(cmath.pi*1j*(phi+theta))* x[1]] for x in posn] )

def shift(coin):
    """Shift the up elements leftwards and the down elements rightwards."""
    newposn = [[0, 0] for i in range(len(coin))]
    for j in range(1, len(posn) - 1):
        newposn[j + 1][0] += coin[j][0]
        newposn[j - 1][1] += coin[j][1]
    return normalise(newposn)


# Initialise lists.
min, max = -STEPS, STEPS+1
qq= 2*STEPS+1

# the coin register is embedded in the position

posn = [[0, 0] for i in range(min, max)]
#posn[-min-20] = [0.5, 1j / 2]
#posn[-min+20] = [0.5, -1j / 2]
#posn[-min] = [1,0]

posn = [ [Gauss( i**2, 20.0 )* cmath.exp(pi2j*float(0.5*i)), 0] for i in range(min, max)]


# Run for some steps...
for time in range(100):
    p = 2
    posn = shift(Hadamard_Coin(posn))
#    posn = shift(General_Coin(posn))
    if time%p==p-1:    # measure after N steps and assume the result is always head
        b = random.randint(0, 1)
        print time, b
        posn = [ [ posn[i][b], 0 ] for i in range(qq)]


posnup = [posn[i][0] for i in range(qq)]
posndo = [posn[i][1] for i in range(qq)]

ft_up = FT(posnup, qq)
ft_do = FT(posndo, qq)

#print posn

pylab.plot(range(min, max), posnup)
#pylab.plot(range(min, max), posndo)
pylab.plot(range(min, max), abs(ft_up))
#pylab.plot(range(min, max), abs(ft_do))
#pylab.plot(range(min, max), probabilities(posn))
pylab.show()
pylab.close()
