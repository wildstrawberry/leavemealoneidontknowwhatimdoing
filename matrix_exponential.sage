
import random
import cmath
pi2j = I

def mexp(dim):
    M = Matrix(RDF,[[(i+j)%2 for i in range(dim)] for j in range(dim)])
    print(M)
    print(pi2j*M) 
    U = exp(pi2j*M)

    print(U)
    print(U.conjugate_transpose())
    print(U*U.conjugate_transpose())
    SVD_matrix = M.SVD()[1]
    for i in range(dim):
        print(i, SVD_matrix[i][i])

mexp(2)
