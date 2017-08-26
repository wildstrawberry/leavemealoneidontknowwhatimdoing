# Hyperelliptic nets
q = 113
FF = FiniteField(q)
R.<x> = PolynomialRing(FF)
hy = x^5 + x^2 + 53*x + 11
C = HyperellipticCurve(hy)
Jac = C.jacobian()
X = Jac(FF)
ch = C.frobenius_polynomial()
print "Frobenius:", ch, " #Jac(C) = ", ch(1), "=", factor(ch(1))

# y2 + (b2x2 + b1x + b0)y = x5 + a4x4 + a3x3 + a2x2 + a1x + a0, default = missing(b2x2 + b1x + b0)
a4 = 0
a3 = 0
a2 = 1
a1 = 53
a0 = 11

P = X([x^2 + 58*x + 94, 68*x + 56])
Q = X([x^2 + 4*x + 24, 38*x + 65])
R = X([x^2 + 4*x + 63, 10*x + 26])

for i in range(8):
    print i, i*P, i*Q, i*R
#print P[0], P[0][1], P[0][0], P[1][1], P[1][0]

def Fg(Pin, Qin):
    """ Uchida, Uchiyama page 222, proposition 9: Let P = (P_1, ..., P_n)
        W[v+w](P)W[v-w](P)/W[v]^2(P)W[w]^2(P) = F_g( [v1]P_1 + ... + [vn]P_n , [w1]P_1 + ... + [wn]P_n ) 
        Uchida, Uchiyama page 229 gives the expression for F_2(P,Q)
    """
    u11, u12, v11, v12 = Pin[0][1], Pin[0][0], Pin[1][1], Pin[1][0]
    u21, u22, v21, v22 = Qin[0][1], Qin[0][0], Qin[1][1], Qin[1][0]
    return -(v11^2 - u11*u12 + u11^3 - a4*u11^2 + a3*u11) + (v21^2 - u21*u22 + u21^3 - a4*u21^2 + a3*u21) - u12*u21 + u11*u22


def W2(u, v, Pin, Qin):
    """ W_{u,v}(P, Q),   for torison groups of size>3, u>v"""

    u11, u12, v11, v12 = Pin[0][1], Pin[0][0], Pin[1][1], Pin[1][0]
    u21, u22, v21, v22 = Qin[0][1], Qin[0][0], Qin[1][1], Qin[1][0]

    # initialize the table, the dimension is set as u+v+7 because the current algorithm first compute W[*][j-1] then W[*][j]; also wanna leave space for W[-1][*]
    W = [[FF(0) for j in range(u+v+7)] for i in range(u+v+7)]

    # compute the column of W[*][0]
    W[1][0] = FF(1)
    W[2][0] = FF( (-4*u12+6*u11^2+(-4*a4)*u11+2*a3)*v12 + 2*v11^3+((-8*u11+4*a4)*u12+2*u11^3 + (-2*a4)*u11^2+(2*a3)*u11-2*a2)*v11 )
    for i in range(3, u+v+5):
        delta = 1
        k = i - delta
        while ( W[k][0]==0 or W[k-delta][0]==0 or W[delta][0]==0 ):
            delta+=1
            k = i - delta
        W[k+delta][0] = Fg( k*Pin , delta*Pin )*W[k][0]^2*W[delta][0]^2/W[k-delta][0]

    # compute the column of W[*][1]
    W[0][1] = FF(1)
    W[1][1] = FF(1)
    for i in range(2, u+v+5):
        delta = 1
        k = i - delta
        while ( W[k][1]==0 or W[k-delta][1]==0 or W[delta][0]==0 ):
            delta+=1
            k = i - delta
        W[k+delta][1] = Fg( k*Pin+Qin , delta*Pin )*W[k][1]^2*W[delta][0]^2/W[k-delta][1]

    W[-1][1] = Fg( Qin, Pin )     # v = (0,1), w = (1,0), so W[-1][1] = Fg( Qin, Pin )*W[0][1]^2*W[1][0]^2/W[1][1]
    # compute the column of W[*][2]
    W[0][2] = Fg( Pin+Qin, Qin-Pin )*W[1][1]^2*W[-1][1]^2/W[2][0]   # v = (1,1), w = (-1,1)
    W[1][2] = Fg( 2*Pin+Qin, Qin-Pin )*W[2][1]^2*W[-1][1]^2/W[3][0]  # v = (2,1), w = (-1,1)
    for i in range(2, u+v+5):
        delta = 1
        k = i - delta
        while ( W[k][2]==0 or W[k-delta][2]==0 or W[delta][0]==0 ):
            delta+=1
            k = i - delta
        W[k+delta][2] = Fg( k*Pin+2*Qin , delta*Pin )*W[k][2]^2*W[delta][0]^2/W[k-delta][2]

    # compute the columns of W[*][j], for j>=3;
    # first compute W[0][j] = Fg( Pin+(j-1)*Qin, Qin-Pin )*W[1][j-1]^2*W[-1][1]^2/W[2][j-2]   # v = (1,j-1), w = (-1,1)
    #  then compute W[1][j] = Fg( 2*Pin+(j-1)*Qin, Qin-Pin )*W[2][j-1]^2*W[-1][1]^2/W[3][j-2]   # v = (2,j-1), w = (-1,1)
    #  then compute W[k+delta][j] = Fg( 2*Pin+(j-1)*Qin, Qin-Pin )*W[2][j-1]^2*W[-1][1]^2/W[3][j-2]   # v = (2,j-1), w = (-1,1)

    for j in range(3,u+v+5):
        delta = 1
        k = j - delta
        while ( W[0][k]==0 or W[0][k-delta]==0 or W[0][delta]==0 ):
            delta+=1
            k = j - delta
        W[0][k+delta] = Fg( k*Qin , delta*Qin )*W[0][k]^2*W[0][delta]^2/W[0][k-delta]   # v = [0][k], w = [0][delta]
        W[1][j] = Fg( 2*Pin+(j-1)*Qin, Qin-Pin )*W[2][j-1]^2*W[-1][1]^2/W[3][j-2]   # v = (2,j-1), w = (-1,1)
        W[2][j] = Fg( 3*Pin+(j-1)*Qin, Qin-Pin )*W[3][j-1]^2*W[-1][1]^2/W[4][j-2]   # v = (3,j-1), w = (-1,1)

        for i in range(3,u+v+5):
            delta = 1
            k = i - delta
            while ( W[k][j]==0 or W[k-delta][j]==0 or W[delta][0]==0 ):
                delta+=1
                k = i - delta
            W[k+delta][j] = Fg( k*Pin+j*Qin , delta*Pin )*W[k][j]^2*W[delta][0]^2/W[k-delta][j]

#    for i in W:
#        print i

    return (W[u][v])**16

#W2(8,8,P,Q)

for i in range(1,7):
    for j in range(1,7):
        print i, j, W2(7, 2, i*P, j*Q)/W2(7, 1, i*P, j*Q)/W2(0, 2, i*P, j*Q)


def Tri(u,v,w, Pin, Qin, Rin):

    u11, u12, v11, v12 = Pin[0][1], Pin[0][0], Pin[1][1], Pin[1][0]
    u21, u22, v21, v22 = Qin[0][1], Qin[0][0], Qin[1][1], Qin[1][0]
    u31, u32, v31, v32 = Rin[0][1], Rin[0][0], Rin[1][1], Rin[1][0]

    W = [[[FF(0) for x in range(3)] for y in range(3)] for z in range(m+5)]
    W[1][0][0] = FF(1)
    W[2][0][0] = FF( (-4*u12+6*u11^2+(-4*a4)*u11+2*a3)*v12 + 2*v11^3+((-8*u11+4*a4)*u12+2*u11^3 + (-2*a4)*u11^2+(2*a3)*u11-2*a2)*v11 )
    for i in range(3, m+2):
        b = i%2
        k = int((i-b)/2)  # i = 2k+b
        #print i, k, b
        if b==0:  # FIXME later
            W[i][0][0] = Fg( (k+1)*Pin, (k-1)*Pin )*W[k+1][0][0]^2*W[k-1][0][0]^2/W[2][0][0]
            #print i, Fg( (k+1)*Pin, (k-1)*Pin ), W[i][0][0]
        else:
            W[i][0][0] = Fg( (k+1)*Pin, (k)*Pin )*W[k+1][0][0]^2*W[k][0][0]^2
            #print i, Fg( (k+1)*Pin, (k)*Pin ), W[i][0][0]

    W[0][1][0] = FF(1)
    W[1][1][0] = FF(1)
    for i in range(2, m+2):
        b = i%2
        k = int((i-b)/2)  # i = 2k+b
        #print i, k, b
        if b==0:  # FIXME later
            W[i][1][0] = Fg( k*Pin + Qin, k*Pin )*W[k][1][0]^2*W[k][0][0]^2
            #print i, Fg( k*Pin + Qin, k*Pin ), W[i][1][0]
        else:
            W[i][1][0] = Fg( (k+1)*Pin + Qin, k*Pin )*W[k+1][1][0]^2*W[k][0][0]^2
            #print i, Fg( (k+1)*Pin + Qin, (k)*Pin ), W[i][0][0]

    W[0][0][1] = FF(1)
    W[1][0][1] = FF(1)
    for i in range(2, m+2):
        b = i%2
        k = int((i-b)/2)  # i = 2k+b
        #print i, k, b
        if b==0:  # FIXME later
            W[i][0][1] = Fg( k*Pin + Rin, k*Pin )*W[k][0][1]^2*W[k][0][0]^2
            #print i, Fg( k*Pin + Rin, k*Pin ), W[i][0][1]
        else:
            W[i][0][1] = Fg( (k+1)*Pin + Rin, k*Pin )*W[k+1][0][1]^2*W[k][0][0]^2
            #print i, Fg( (k+1)*Pin + Rin, (k)*Rin ), W[i][0][1]

    W[0][1][1] = FF(1)
    W[1][1][1] = FF(1)  #this is not correct
    for i in range(2, m+2):
        b = i%2
        k = int((i-b)/2)  # i = 2k+b
        #print i, k, b
        if b==0:  # FIXME later
            W[i][1][1] = Fg( k*Pin + Qin + Rin, k*Pin )*W[k][1][1]^2*W[k][0][0]^2
            #print i, Fg( k*Pin + Qin + Rin, k*Pin ), W[i][1][1]
        else:
            W[i][1][1] = Fg( (k+1)*Pin + Qin + Rin, k*Pin )*W[k+1][1][1]^2*W[k][0][0]^2
            #print i, Fg( (k+1)*Pin + Qin + Rin, (k)*Rin ), W[i][1][1]

    #print "Manually compute W given m=", m, ", P = ", Pin, ", Q = ", Qin, ", R = ", Rin
    for i in range(8):
        print i, W[i][0][0], W[i][1][0], W[i][0][1], W[i][1][1]
    return (W[7][1][0])**16
