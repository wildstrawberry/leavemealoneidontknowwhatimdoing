#from sage.stats.distributions.discrete_gaussian_integer import DiscreteGaussianDistributionIntegerSampler
#from sage.modules.free_module_integer import IntegerLattice
import random
R = ZZ  # the base ring
F2 = FiniteField(2)

def rand_full_rank_MatR(k,l):
    """ generate a random full rank kxl matrix in R """
    M = Matrix(R,[[R.random_element() for i in range(l)] for j in range(k)])
    while M.rank()<k and M.rank()<l:
        M = Matrix(R,[[R.random_element() for i in range(l)] for j in range(k)])
    return M

def rand_rank_one_MatR(k,l):
    """ generate a random kxl matrix of rank 1 """
    M1 = Matrix(R,[[R.random_element() for i in range(1)] for j in range(k)])
    M2 = Matrix(R,[[R.random_element() for i in range(l)] for j in range(1)])
    return M1*M2

def rand_full_rank_B(k,l,B):
    """ generate a random full rank matrix whose entries are bounded by [-B, B] """
    M = Matrix(ZZ,[[random.randint(-B, B) for i in range(l)] for j in range(k)])
    while M.rank()<k and M.rank()<l:
        M = Matrix(ZZ,[[random.randint(-B, B) for i in range(l)] for j in range(k)])
    return M

def SVD_binary(dim):
    """ experiments with ABB10 difference  """
    M = rand_full_rank_MatR(dim, dim)
    print M, M.det()
    print M.eigenvalues()
    print M.change_ring(RDF).SVD()

SVD_binary(4)

def funfacts_symmatrix():
    """ to verify Nick's equation """
    R = Matrix(RR, [ [1, 1, 0], [2, 0, 1] ] )
    M = R*(R.transpose())
    print R
    print M
    I2 = matrix.identity(2)
    A = M+I2
    B = M+2*I2
    C = B.inverse()
    print A
    print C
    return A*A*(C) - (C)*A*A

def funfacts_unitary():
    R = Matrix(RR, [ [1, 1, 1, 1], [1, 1, -1, -1], [1, -1, 1, -1], [1, -1, -1, 1] ] )
    M = R*(R.transpose())
    print R
    print M

def ABB10(u0, u1, u2, u3):
    """ the Full-rank-difference mapping in ABB10 """
    M = Matrix(ZZ, [[u0, u1, u2, u3], [u3, u0-u3, u1, u2], [u2, u3-u2, u0-u3, u1], [u1, u2-u1, u3-u2, u0-u3]] )
    return M

def diff_ABB10(B):
    """ experiments with ABB10 difference  """
    M1 = ABB10(1, random.randint(-B, B), random.randint(-B, B), 0)
    M2 = ABB10(1, random.randint(-B, B), random.randint(-B, B), 0)
    print M1, M1.det()
    print M2, M2.det()
    print M1-M2, (M1-M2).det()

def toy():
    w = 3
    Iw = matrix.identity(w)
    Pw = Matrix(ZZ,[[ZZ(((i+1)%w)==j) for i in range(w)] for j in range(w)])  # some non-identity matrix P of dimension w
    Pinvw = Matrix(ZZ,[[ZZ(((i-1)%w)==j) for i in range(w)] for j in range(w)])  # P^{-1}
    M1 = rand_full_rank_B(w,w,1)
    M1inv = M1^(-1)
    M2 = rand_full_rank_B(w,w,1)
    M2inv = M2^(-1)
    M3 = rand_full_rank_B(w,w,1)
    print M1, M1.charpoly()
    print M1inv, M1inv.charpoly()
    print M2, M2.charpoly()
    print M2inv, M2inv.charpoly()
    print (M3*M1*M2).charpoly()
    print (M3*M2*M1).charpoly()
