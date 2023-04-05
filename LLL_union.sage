from sage.stats.distributions.discrete_gaussian_integer import DiscreteGaussianDistributionIntegerSampler
from sage.modules.free_module_integer import IntegerLattice
import random

def biased_coin(prob):
    """ output 1 with probability prob  """
    rand = random.random()
    if rand < prob: return 1
    else: return 0

def LA_basis_gen(R, P, q, n, m):
    """ generate a basis for L(A) = { y | y = sAP mod R }, for a random A """
    k = m-n
    M = Matrix(ZZ,[[random.randint(0, q) for i in range(k)] for j in range(n)])
    I_n = matrix.identity(n) 
    PI_n =P*I_n 
    RI_k =R*matrix.identity(k) 
    Z = Matrix(ZZ,[[0 for i in range(n)] for j in range(k)])

    errb = P
    Pnoise1 = Matrix(ZZ,[[random.randint(-errb, errb) for i in range(m-n)] for j in range(n)])
    Pnoise2 = Matrix(ZZ,[[random.randint(-errb, errb) for i in range(n)] for j in range(n)])

    PI_nwide = block_matrix( [ [PI_n, PI_n, PI_n, PI_n] ] )

    MT = block_matrix( [ [RI_k, Z] , [P*M, PI_n], [ Pnoise1, Pnoise2] ] )
    #print(MT)
    #Gram = MT.T * MT
    #print(Gram.det())
    L = IntegerLattice(MT)
    return L

#The classical Gran-Schmidt algorithm (in truth a very common variant, known for its stability)
#This algorithm performs operations over the field F. One can typically take F to be QQ or RR.
def Stable_GS(A,F):
    N = A.nrows();
    M = A.ncols();
    B = matrix(F,N,M);
    v = vector(F,N);
    for i in range(N):
        B[i] = A[i];
        for j in range(i):
            B[i] = B[i] - (B[i]*B[j])/(v[j])*B[j];
        v[i] = B[i]*B[i];
    return B;

def GS_length(q, G):
    print("x,   G[x].norm(), G[x].norm()^2 " )
    for x in range(q):
        print(x, RR(G[x].norm()), RR(G[x].norm())**2 )

def main():
    P = 100
    q = 27
    R = P*q

    n = 3
    m = n*3
    L = LA_basis_gen(R, P, q, n, m)
    print(L)

    LM = Matrix(L.basis())
    print("determinant:  ", LM.det())
    print("det predicted:",  (P**(m-n)) * q**(m-2*n))
    G = Stable_GS(LM, RR)
    GS_length(m, G)
main()
