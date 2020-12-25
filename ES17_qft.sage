import random

def Gen_random_matrix_mod_q(n,m, q):
    return Matrix(ZZ,[[ Integers(q).random_element() for i in range(m)] for j in range(n)])
    
def ES17_Fourier(n, m, q):
    print("n,m,q:", n,m,q)
    a = Integers(q).random_element()
    b = Integers(q).random_element()
    A = Matrix(ZZ,[[ 1, 1, 0] ])  #Gen_random_matrix_mod_q(n,m, q)  
    B = Matrix(ZZ,[[ 1, 0, 1] ])
    C = A 
    D = B #+ Matrix(ZZ,[[ 1, 1, 0] ])
    print("A=", A)
    print("B=", B)
    print("C=", C)
    print("D=", D)
    M = Matrix( [[A*C.T, A*D.T], [B*C.T, B*D.T]] )
    print("M = ", M, "determinant of M:", M.determinant() )
    U = Matrix( CDF, [[ exp(2*pi*I*  ( ( ((i-i%q)/q)*A + (i%q)*B )*( ((j-j%q)/q)*C.T + (j%q)*D.T) )[0][0] /q)/q for i in range(q^2)] for j in range(q^2)])
    
    print("the determinant of U:", U.determinant())
    print("U^2:")
    print(U*U)
    list_eigenvalues = [0,0,0,0]
    eigen = U.eigenvalues()
    print("eigenvalues:")
    for evs in eigen:
        if "%0.4f"%evs[0]=='1.0000': list_eigenvalues[0]+=1
        elif "%0.4f"%evs[0]=='-1.0000': list_eigenvalues[2]+=1
        elif "%0.4f"%evs[1]=='1.0000': list_eigenvalues[1]+=1
        elif "%0.4f"%evs[1]=='-1.0000': list_eigenvalues[3]+=1
        #print("%0.4f"%evs[0], "%0.4f"%evs[1])
    print(list_eigenvalues)
    #eigenvectors = U.eigenvectors_right()
    #for evs in eigenvectors:
    #    print(evs[0].n(prec=10), evs[1][0].n(prec=10))

ES17_Fourier(1, 3, 17)

  
def normal_QFT(q):
    print("q=", q)
    U = Matrix( CDF, [[ exp(2*pi*I*i*j/q)/sqrt(q) for i in range(q)] for j in range(q)])
    print("U:")
    print(U)
    print("U^2:")
    print(U*U)
    print("the determinant of U:", U.determinant())
    #eigen = U.eigenvalues()
    #print("eigenvalues:")
    #for evs in eigen:
    #    print("%0.4f"%evs[0], "%0.4f"%evs[1])
    eigenvectors = U.eigenvectors_right()
    for evs in eigenvectors:
        print(evs[0].n(prec=10), evs[1][0].n(prec=10))
    
#normal_QFT(4)
