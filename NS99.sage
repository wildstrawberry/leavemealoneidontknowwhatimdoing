# "Visually" a lattice vector from a basis B is v = xB, where v and x are visually row vectors
# Of course the correct way to look at them is to take each row vector as a column vector.
from sage.stats.distributions.discrete_gaussian_integer import DiscreteGaussianDistributionIntegerSampler
from sage.modules.free_module_integer import IntegerLattice
import random

B_small = 10
dim = 50  # dim >= 10
dimprint = 3  # dimprint < dim
dim2 = 5#dim
dim3 = dim+5
B_large = (B_small)^(dim)
print B_small, dim, B_large

w = 4
Pw = Matrix(ZZ,[[ZZ(((-i)%w)==j) for i in range(w)] for j in range(w)])  # some non-identity matrix P of dimension w
print Pw

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
    
def NS99exp():
    """
        From an m-dim vector t = X \cdot s,
        where X\in ZZ^{m\times n} has small entries,
        s is an n-dimensional vector (not necessarily have small entries), extension: allow s to be a matrix
        trying to recover X and s from b using Nyugen-Stern 99
    """
    X = rand_full_rank_B(dim, dim3, B_small)
    X[1] = X[0]  # creating a linearly dependent vector
    #X[2] = X[0]
    #X[2][0] = X[2][0]+1  # creating a slightly linearly independent vector

    #print X
    s = rand_full_rank_B(dim2, dim, B_large)
    #print s
    t = s * X
    #print t
    L = IntegerLattice( X  )
    print "L = ", L.basis()[0:dimprint] # already LLL reduced
    print "volume^2 of L:", L.volume()^2, len(str(L.volume()^2)) # volumn = determinant

    Lt = IntegerLattice( t  )
#    print "Lt = ", Lt.basis()[0:dimprint] # already LLL reduced
    print "volume^2 of Lt:", Lt.volume()^2, len(str(Lt.volume()^2)) # volumn = determinant


    tbot = kernel(t.transpose()).matrix()
    print "rank of the kernal of t", tbot.rank()
    Lbot = IntegerLattice( tbot  )
    print "Lbot = ", Lbot
    Bbot_firstfew = Lbot.basis()[0:(dim3-dim)]
#    print "the basis of the first few vectors of Lbot:", Bbot_firstfew
    Lbot_firstfew = IntegerLattice( Bbot_firstfew )
    print "volume^2 of the first few vectors of Lbot ", Lbot_firstfew.volume()^2, len(str(Lbot_firstfew.volume()^2))

#    Xbar = kernel(Matrix(Bbot_firstfew).transpose()).matrix()
#    print Xbar
#    Lbar = IntegerLattice( Xbar  )
#    print "Lbar = ", Lbar.basis()[0:dimprint]
#    print "volume^2 of Lbar:", Lbar.volume()^2, len(str(Lbar.volume()^2))

#NS99exp()
