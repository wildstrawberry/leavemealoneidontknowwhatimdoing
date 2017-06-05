# Weil Pairing Example by elliptic division sequence
# https://eprint.iacr.org/2006/392.pdf
# E: y^2 = x^3 + 30x + 34 mod 631
# P = E((36, 60))
# Q = E((121, 387))

q = 631
A = 30
B = 34
E = EllipticCurve(GF(q), [30, 34])

TH = 40

FF = FiniteField(q)
R.<x> = PolynomialRing(FF)

def double(V, W20, W11, W_11, W2_1):
    X = [0]+V[0]
    Y = V[1]
    Vnew = [[0,0,0,0,0,0,0,0],[0,0,0]]
    Vnew[0][0] = X[4]*X[2]^3 - X[1]*X[3]^3
    Vnew[0][2] = X[5]*X[3]^3 - X[2]*X[4]^3
    Vnew[0][4] = X[6]*X[4]^3 - X[3]*X[5]^3
    Vnew[0][6] = X[7]*X[5]^3 - X[4]*X[6]^3
    Vnew[0][1] = X[3]*(X[5]*X[2]^2-X[1]*X[4]^2)/W20
    Vnew[0][3] = X[4]*(X[6]*X[3]^2-X[2]*X[5]^2)/W20
    Vnew[0][5] = X[5]*(X[7]*X[4]^2-X[3]*X[6]^2)/W20
    Vnew[0][7] = X[6]*(X[8]*X[5]^2-X[4]*X[7]^2)/W20

    Vnew[1][0] = (Y[2]*Y[0]*X[3]^2 - X[4]*X[2]*Y[1]^2)/W11
    Vnew[1][1] = (Y[0]*Y[2]*X[4]^2 - X[3]*X[5]*Y[1]^2)
    Vnew[1][2] = (Y[0]*Y[2]*X[5]^2 - X[4]*X[6]*Y[1]^2)/W_11
    print "double:", Vnew
    return Vnew

def doubleadd(V, W20, W11, W_11, W2_1):
    X = [0]+V[0]
    Y = V[1]
    Vnew = [[0,0,0,0,0,0,0,0],[0,0,0]]
    Vnew[0][1] = X[5]*X[3]^3 - X[2]*X[4]^3
    Vnew[0][3] = X[6]*X[4]^3 - X[3]*X[5]^3
    Vnew[0][5] = X[7]*X[5]^3 - X[4]*X[6]^3
    Vnew[0][7] = X[8]*X[6]^3 - X[5]*X[7]^3
    Vnew[0][0] = X[3]*(X[5]*X[2]^2-X[1]*X[4]^2)/W20
    Vnew[0][2] = X[4]*(X[6]*X[3]^2-X[2]*X[5]^2)/W20
    Vnew[0][4] = X[5]*(X[7]*X[4]^2-X[3]*X[6]^2)/W20
    Vnew[0][6] = X[6]*(X[8]*X[5]^2-X[4]*X[7]^2)/W20

    Vnew[1][0] = (Y[0]*Y[2]*X[4]^2 - X[3]*X[5]*Y[1]^2)
    Vnew[1][1] = (Y[0]*Y[2]*X[5]^2 - X[4]*X[6]*Y[1]^2)/W_11
    Vnew[1][2] = (X[5]*X[7]*Y[1]^2 - Y[0]*Y[2]*X[6]^2)/W2_1
    print "doubleadd:", Vnew
    return Vnew

def katepairing(m, x1, y1, x2, y2):
    mbits = bin(m)[3:]  # bin(m)=0b... string type
    n = len(mbits)

    # the range are set to be larger to incoorperate -1 indices, W is consistant with the notation in https://eprint.iacr.org/2006/392.pdf
    W = [[FF(0) for x in range(3)] for y in range(m+5)]
    W[1][0] = FF(1)
    W[2][0] = FF(2*y1)
    W[3][0] = FF(3*(x1^4) + 6*A*(x1^2) + 12*B*x1 - A^2)
    W[4][0] = FF(4*y1*(x1^6 + 5*A*(x1^4) + 20*B*(x1^3) - 5*A^2*x1^2 - 4*A*B*x1 - 8*B^2 - A^3 ))
    W[0][1] = FF(1)
    W[1][1] = FF(1)
    W[2][1] = FF(2*x1+x2-((y2-y1)/(x2-x1))^2)
    W[-1][1] = FF( x1-x2 )
    W[2][-1] = FF( (y1+y2)^2-(2*x1+x2)*(x1-x2)^2 )

    # initialize a block centered at W[2][0]
    V = [[-W[2][0], -1, 0, 1, W[2][0], W[3][0], W[4][0], ((W[2][0])^3)*W[4][0] - (W[3][0])^3] , [1, W[1][1], W[2][1]]]

    for i in mbits:
        if i=='0':
            V = double(V, W[2][0], W[1][1], W[-1][1], W[2][-1])
        else:
            V = doubleadd(V, W[2][0], W[1][1], W[-1][1], W[2][-1])

    return V[0][3], V[1][1], V[1][1]/V[0][3]

def katepairingmanual(m, x1, y1, x2, y2):
    # compute inefficiently
    W = [[FF(0) for x in range(3)] for y in range(m+5)]
    W[1][0] = FF(1)
    W[2][0] = FF(2*y1)
    W[3][0] = FF(3*(x1^4) + 6*A*(x1^2) + 12*B*x1 - A^2)
    W[4][0] = FF(4*y1*(x1^6 + 5*A*(x1^4) + 20*B*(x1^3) - 5*A^2*x1^2 - 4*A*B*x1 - 8*B^2 - A^3 ))
    W[0][1] = FF(1)
    W[1][1] = FF(1)
    W[2][1] = FF(2*x1+x2-((y2-y1)/(x2-x1))^2)
    W[-1][1] = FF( x1-x2 )
    W[2][-1] = FF( (y1+y2)^2-(2*x1+x2)*(x1-x2)^2 )

    for i in range(5,m+1):
        b = i%2
        k = int((i+b)/2)
        print i, k, b
        if b==0:
            W[i][0] = W[k][0]*(W[k+2][0]*W[k-1][0]^2 - W[k-2][0]*W[k+1][0]^2)/W[2][0]
        else:
            W[i][0] = W[k+1][0]*W[k-1][0]^3 - W[k-2][0]*W[k][0]^3
            
    for i in range(3,m+1):
        b = i%2
        k = int((i-b)/2)
        print i, k, b
        if b==1:
            W[i][1] = (W[k-1][1]*W[k+1][1]*W[k+1][0]^2 - W[k][0]*W[k+2][0]*W[k][1]^2)/W[-1][1]
        else:
            W[i][1] = W[k-1][1]*W[k+1][1]*W[k][0]^2 - W[k-1][0]*W[k+1][0]*W[k][1]^2

    return W, W[m][0], W[m][1], W[m][1]/W[m][0]

print katepairing(6, 121, 387, 36, 60)
print katepairingmanual(6, 121, 387, 36, 60)
print katepairing(6, 36, 60, 121, 387)
print katepairingmanual(6, 36, 60, 121, 387)
