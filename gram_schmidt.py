import numpy as np
import math

# Gram-Schmidt Orthogonalization
def gs(X, row_vecs=True, norm = False):
    if not row_vecs:
        X = X.T
    Y = X[0:1,:].copy()
    for i in range(1, X.shape[0]):
        proj = np.diag((X[i,:].dot(Y.T)/np.linalg.norm(Y,axis=1)**2).flat).dot(Y)
        Y = np.vstack((Y, X[i,:] - proj.sum(0)))
    if norm:
        Y = np.diag(1/np.linalg.norm(Y,axis=1)).dot(Y)
    if row_vecs:
        return Y
    else:
        return Y.T

def GS_prod_check(GS):
    prod = 1
    for gsv in GS.T:
        prod*=np.linalg.norm(gsv)
    return prod


def Gaussian_pdf(x, mu, sigma):
    x = float(x - mu) / sigma
    return math.exp(-x*x/2.0) / math.sqrt(2.0*math.pi) / sigma


def Gaussian_error_state_matrix(sigma, n):
    M = np.matrix( [[(Gaussian_pdf(i, j, sigma)+Gaussian_pdf(i-n, j, sigma)+Gaussian_pdf(i+n, j, sigma)) for i in range(n)] for j in range(n)])
    print(sigma, n)
    print(M)
    return M

if __name__ == "__main__":
    B = 1.0
    n = 17
    GS = gs( Gaussian_error_state_matrix( B, n) )
    print("norm of Gram_Schmidt(M)", [ np.linalg.norm(gsv) for gsv in GS ])
    print("norm of GS(first row)/norm of GS(last row) ", np.linalg.norm(GS[0])/np.linalg.norm(GS[-1]) )

