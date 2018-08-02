from sage.stats.distributions.discrete_gaussian_integer import DiscreteGaussianDistributionIntegerSampler
from sage.modules.free_module_integer import IntegerLattice
import random

R = ZZ  # the base ring
F2 = FiniteField(2)

w=4

Iw = matrix.identity(w)
Pw = Matrix(ZZ,[[ZZ(((i+1)%w)==j) for i in range(w)] for j in range(w)])  # some non-identity matrix P of dimension w
Pinvw = Matrix(ZZ,[[ZZ(((i-1)%w)==j) for i in range(w)] for j in range(w)])  # P^{-1}

print "I =", Iw
print "P =", Pw
print "Pinverse =", Pinvw
print "P*P^{-1} =", Pw*Pinvw

Iw2 = matrix.identity(w)
Pw2 = Matrix(F2,[[ZZ(((i+1)%w)==j) for i in range(w)] for j in range(w)])  # some non-identity matrix P of dimension w
Pinvw2 = Matrix(F2,[[ZZ(((i-1)%w)==j) for i in range(w)] for j in range(w)])  # P^{-1}

print "I =", Iw2
print "P =", Pw2
print "Pinverse =", Pinvw2
print "P*P^{-1} =", Pw2*Pinvw2

print Iw2.charpoly(), Iw2.eigenvalues()
print Pw2.charpoly(), Pw2.eigenvalues()


#MMMM = block_matrix( [[Iw, Iw],[Iw, Iw]] )
#NNNN = block_matrix( [[Iw, Pinvw],[Pw, Iw]] )


def rand_full_rank_B(k,l,B):
    """ generate a random full rank matrix whose entries are bounded by [-B, B] """
    M = Matrix(ZZ,[[random.randint(-B, B) for i in range(l)] for j in range(k)])
    while M.rank()<k and M.rank()<l:
        M = Matrix(ZZ,[[random.randint(-B, B) for i in range(l)] for j in range(k)])
    return M

M1 = rand_full_rank_B(w,w,1)
M1inv = M1^(-1)
M2 = rand_full_rank_B(w,w,1)
M2inv = M2^(-1)

print M1, M1.charpoly()
print M1inv, M1inv.charpoly()
print M2, M2.charpoly()
print M2inv, M2inv.charpoly()

print (M1*M2).charpoly()

print (M2*M1).charpoly()
