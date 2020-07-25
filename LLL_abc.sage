from sage.stats.distributions.discrete_gaussian_integer import DiscreteGaussianDistributionIntegerSampler
from sage.modules.free_module_integer import IntegerLattice
import random

B_small = 5000
dim = 130  # dim >= 10
dimprint = 3  # dimprint < dim
dim2 = 5#dim
dim3 = dim+5
B_large = (B_small*dim)^3#^(dim/2)
print B_small, dim, B_large

def biased_coin(prob):
    """ output 1 with probability prob  """
    rand = random.random()
    if rand < prob: return 1
    else: return 0

def rand_full_rank_B(k,l,B):
    """ generate a random full rank matrix whose entries are bounded by [-B, B] """
    M = Matrix(ZZ,[[random.randint(-B, B) for i in range(l)] for j in range(k)])
    while M.rank()<k and M.rank()<l:
        M = Matrix(ZZ,[[random.randint(-B, B) for i in range(l)] for j in range(k)])
    return M

def rand_matrix_Hamming(k,l,prob):
    """ generate a random matrix with hamming weight prob* dim """
    M = Matrix(ZZ,[[biased_coin(prob) for i in range(l)] for j in range(k)])
    return M

def special(k,B):
    """ generate a special basis """
    M = Matrix(ZZ,[[ (int(i==j) + random.randint(-B, B)*int(i==k-1) ) for i in range(k)] for j in range(k)])
    return M

def attack_knapsack(n,w):
    """ assume the knapsack instance is given by sum_{ai, xi} = 1, of n entries, each ai is of weight w """
    """ using the LLL attack on knapsack of , see https://web.eecs.umich.edu/~cpeikert/lic13/lec05.pdf """
    A = [ random.randint(1, w) for i in range(n) ]
    print A, sum(A)
    A = A + [sum(A)]
    B = n*2^(n/2)
    M = Matrix(ZZ,[[ (int(i==j) + (-B*A[j])*int(i==n and j!=n) ) + (B*A[j]-1)*int(i==j and i ==n)   for i in range(n+1)] for j in range(n+1)])
    print M
    L = IntegerLattice( M )
    print "L = ", L # already LLL reduced, oops
    print L.volume()^2
    x = L.shortest_vector()
    print x

#attack_knapsack(4, 2^16)

def test_LLL_toy():
    #L = IntegerLattice( [ [2,8,0, 0],[3,1,0,0], [0,0,1,2], [0,0,32,1] ]  )
    #L = IntegerLattice( [ [2,8],[32,72], [20,10] ]  )
    #Basis = rand_matrix_Hamming(40, 60, 0.5)
    Basis = special(dim, B_large)
    #print Basis
    L = IntegerLattice( Basis )
    print "L = ", L # already LLL reduced, oops
    print L.volume()^2
    print L.shortest_vector()
    #aList = [3,41,5,6,6,26,7]
    #aList.sort()
    #print aList
test_LLL_toy()

def lattice_intersection():
    """ produce the intersection of two lattices by taking the dual, then the union, then the dual...oops there is a button """
    L1 = IntegerLattice( [ [2, -1, 0, 0, 0],[0, 2, -1, 0, 0], [0,0,2, -1, 0 ],[0,0,0,2, -1 ], [1,1,1,1,1] ]  )
    L2 = IntegerLattice( [ [2, -1, 0,0, 0],[0, 2, -1,0, 0], [0,0,0, 1,0], [0,0,0,0, 1] ]  )
    #    L2 = IntegerLattice( [ [1, 0, 0,0, 0],[0, 0, 1,0, 0],[0, 1, 0,0, 0]]  )
    print "L1 = ", L1 # the basis in the output is already LLL reduced
    print "L2 = ", L2
    print L1.intersection(L2)
    B12 = L2.intersection(L1).basis()
    M12 = Pw*Matrix(B12)
    M12T = M12.transpose()
    M12Tdual = M12T*(( M12T.transpose()*M12T ).inverse())
    M12dual = M12Tdual.transpose()
    print M12dual

    B12dual = IntegerLattice( 21*31*M12dual )
    print "L12 dual:", B12dual
    print "The basis of L1 intersect L2:", M12, "the norm:", M12.norm()
    G, M = M12.gram_schmidt()
    print "G of the basis:", G, "the norm:", G.norm()
    print "M of the basis:", M, "the norm:", M.norm()
    
lattice_intersection()

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

def GS_Gadget():
    """ Gram-Schmidt of Gadget matrices   """
    b = 17
    B1 = Matrix(ZZ,[ [b-13,b-12,b-23,b-30,1],[b, -1, 0, 0, 0],[0, b, -1, 0, 0], [0,0,b, -1, 0 ], [0,0,0,b, 0 ] ])
    L1 = IntegerLattice( B1 )
    print "L1 = ", L1 # the basis in the output is already LLL reduced
    G, M = B1.gram_schmidt()
    for x in G:
        print x, RR(x.norm())
    B_LLL = Matrix(L1.basis())
    print B_LLL
    G_LLL, M_LLL = B_LLL.gram_schmidt()
    for x in G_LLL:
        print x, RR(x.norm())
