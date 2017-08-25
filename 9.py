from sage.stats.distributions.discrete_gaussian_integer import DiscreteGaussianDistributionIntegerSampler
from sage.modules.free_module_integer import IntegerLattice

def test_LLL_toy():
    L = IntegerLattice( [ [2,8,0, 0],[3,1,0,0], [0,0,1,2], [0,0,32,1] ]  )
    L2 = IntegerLattice( [ [2,8],[32,72], [20,10] ]  )
    print "L = ", L2
    print L2.volume()
    print L2.shortest_vector()


def test_SVP(m):
    sigma1 = m
    sigma2 = 1.5*m
    sigma3 = 5*m*m
    N = m+2
    D1 = DiscreteGaussianDistributionIntegerSampler(sigma=sigma1)
    D2 = DiscreteGaussianDistributionIntegerSampler(sigma=sigma2)
    D3 = DiscreteGaussianDistributionIntegerSampler(sigma=sigma3)
    M = Matrix([[D2() for i in range(N)] for j in range(m)])
    E = Matrix([[D3() for j in range(m)]])
    S = Matrix([[D1() for j in range(N)]])
    T = E*M+S
    MT = block_matrix( [[M], [T]] )
    L = IntegerLattice( MT  )
    print "M = ", M
    print "E = ", E
    print "S = ", S[0][0:10], RR(S[0].norm())
    print "The norm of T = EM+S = ", RR(T[0].norm())
    print "L = ", L, "\n det(L):", L.volume()
    LLLB = L.LLL()
    for i in range(10):
        print i, RR(LLLB[i].norm()), LLLB[i]
    print "shortest vector (this possibly runs in superpoly time)", L.shortest_vector()[0:10]

# input m>=10 is the dimension
test_SVP(15)


def test_inverse_Gaussian_matrix(m):
    """ test the expected norm of C^{-1}, where C is an mxm matrix, with each entry from gaussian with sigma = m """
    sigma = m
    D = DiscreteGaussianDistributionIntegerSampler(sigma=sigma)
    M = Matrix([[D() for i in range(m)] for j in range(m)])
    MM = M^(-1)
    for y in range(8*m):
        M2 = Matrix([[D() for i in range(m)] for j in range(m)])
        MM = MM + M2^(-1)
    print m, sigma, M.norm(), (M^(-1)).norm(), MM.norm()

def preimage_toy(d, A, F):
    """ d: the dimension, A and F: function and image; output M s.t. AM = F"""
    D = []
    for i in range(d):
        dc = randcolumn(d)
        z = F[0][i]
        for j in range(1,d):
            z = z - A[0][j]*dc[j][0]
        dc[0] = z/A[0][0]
        D.append(dc)
    MinZ = block_matrix(ZZ, 1,d,[ a for a in D ])
    return MinZ
