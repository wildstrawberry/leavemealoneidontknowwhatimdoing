import numpy as np
import math

def Gaussian_pdf(x, mu, sigma):
    x = float(x - mu) / sigma
    return math.exp(-x*x/2.0) / math.sqrt(2.0*math.pi) / sigma

def Gaussian_shift_n(x, mu, sigma, q):
    # simulating Gaussian over integer mod q
    return Gaussian_pdf(x, mu, sigma)+Gaussian_pdf(x-q, mu, sigma)+Gaussian_pdf(x+q, mu, sigma)+Gaussian_pdf(x-2*q, mu, sigma)+Gaussian_pdf(x+2*q, mu, sigma)

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
    
def Gaussian_error_state_matrix(sigma, q):
    M = np.matrix( [[ Gaussian_shift_n(i, j, sigma, q) for i in range(q)] for j in range(q)])
    print(sigma, q)
    #print(M)
    return M

if __name__ == "__main__":
    for db in range(1, 10):
        for q in range(17, 18):
            GS = gs( Gaussian_error_state_matrix( 1.5+0.0000000000014*db, q) )
            print("norm of GS(first row)/norm of GS(last row) ", np.linalg.norm(GS[0])/np.linalg.norm(GS[-1]) )
            print("norm of Gram_Schmidt(M)", [ np.linalg.norm(gsv) for gsv in GS ])
            #print("norm of Gram_Schmidt(M), first 2 and last 2", np.linalg.norm(GS[0]), np.linalg.norm(GS[1]), np.linalg.norm(GS[-2]), np.linalg.norm(GS[-1]))
            
