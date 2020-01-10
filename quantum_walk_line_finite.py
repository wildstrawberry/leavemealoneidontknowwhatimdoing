"""
Two walker quantum walk on a finite line
Modified from https://gist.github.com/qubyte/3768908
"""

import math
import cmath
import pylab

STEPS = 200
LEN = 2*STEPS+1

def probabilities(posn):
    """Returns a list of the probabilies for each place."""
    return [sum([abs(amp) ** 2 for amp in place]) for place in posn]

def normalise(posn):
    """Normalise function to normalise an input 1D line."""
    N = math.sqrt(sum(probabilities(posn)))
    return [[amp / N for amp in place] for place in posn]


def Hadamard_Coin(posn):
    """At each timestep, applying a Hadamard gate on each element."""
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
    for j in range(0, LEN ):
        newposn[(j + 1)%LEN][0] += coin[j][0]
        newposn[j - 1][1] += coin[j][1]
    return normalise(newposn)


# Initialise lists.
posn = [[0, 0] for i in range(0, LEN)]

posn[STEPS-40] = [0.5, 1j / 2]
posn[STEPS+40] = [0.5, 1j / 2]

# Run for some steps...
for period in range(120):
    for time in range(9):
        posn = shift(Hadamard_Coin(posn))
    #    posn = shift(General_Coin(posn))
    print period, sum(probabilities(posn)[STEPS-40: STEPS+40]), probabilities(posn)[STEPS-40], probabilities(posn)[STEPS+40]
    
pylab.plot(range(0, LEN), probabilities(posn))
pylab.show()
