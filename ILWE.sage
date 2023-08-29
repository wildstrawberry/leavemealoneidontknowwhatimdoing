from sage.stats.distributions.discrete_gaussian_integer import DiscreteGaussianDistributionIntegerSampler
from sage.modules.free_module_integer import IntegerLattice
import random

# proof of concept test of Integer LWE (LWE without mod q), and its homogeneous version.

def biased_coin(prob):
    """ output 1 with probability prob  """
    rand = random.random()
    if rand < prob: return 1
    else: return 0

def gaussianvector(n, sigma):
    D = DiscreteGaussianDistributionIntegerSampler(sigma = sigma)
    v = vector(ZZ, [ D() for i in range(n)] )
    return v

def unitvector(n, j):
    v = vector(ZZ, [ j==i for i in range(n)] )
    return v

def genA(n, m, sigma):
    D = DiscreteGaussianDistributionIntegerSampler(sigma = sigma)
    A = Matrix(ZZ, [[ D() for i in range(n)] for j in range(m)] )
    return A

def genlow_rank_A(sec, n, m, sigma_A, sigma_e):
    # generate A such that A* sec = 0
    D = DiscreteGaussianDistributionIntegerSampler(sigma = sigma_A)
    A = Matrix(ZZ, [[ 0 for i in range(n)] for j in range(m)] )

    j = 0
    while j<m:
        v = vector(ZZ, [ D() for i in range(n)] )
        if abs(v*sec) < sigma_e:
            A[j] = v
            j+=1

    return A

def roundM(M,param, n):
    N = Matrix(ZZ, [[ round(M[i][j]/param) for i in range(n)] for j in range(n)] )
    return N

def truncM(M, n):
    N = Matrix(ZZ, [[ round(M[i][j]) for i in range(n-1)] for j in range(n-1)] )
    return N


def main():

    n = 5
    m = 100
    sigma_s = 1.0
    sigma_e = 7.0
    sigma_A = 5.0
    param = 100
    sec = gaussianvector(n, sigma_s)
    print("The LWE secret: s", sec)

    A = genlow_rank_A(sec, n, m, sigma_A, sigma_e)
    print("LWE samples", A*sec)
    M = A.T*A
    print("M = A^T *A, and M * s:", M, M*sec)
    print("eigenvalues of M:", M.eigenvalues())

    '''
    N = roundM(M,param, n)
    print(N, N*sec)
    print(N.eigenvalues())
    '''

    u = unitvector(n, n-1)
    Mu = M*u
    print(u, Mu)

    trunc = Mu[0:n-1]
    print(trunc)

    N = truncM(M, n)
    inv = N.inverse()*trunc

    inv = vector(RR, [ inv[i] for i in range(n-1)] )
    print(inv)

main()
