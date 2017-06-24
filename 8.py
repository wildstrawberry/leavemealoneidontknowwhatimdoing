# Compute hyperelliptic psi function Psi_{n,m}(P, Q), when Q is possibly in \Theta

q = 113^2
FF = FiniteField(q)
RR.<x> = PolynomialRing(FF)
hy = x^5 + x^2 + 53*x + 11
C = HyperellipticCurve(hy)
NumberofpointsonC = C.count_points(1)
print "number of points:", NumberofpointsonC, "factors:", factor(NumberofpointsonC[0])
J = C.jacobian()(FF)
X = J
ch = C.frobenius_polynomial()
print "Frobenius:", ch, " #J(C) = ", ch(1), factor(ch(1))

# y2 + (b2x2 + b1x + b0)y = x5 + a4x4 + a3x3 + a2x2 + a1x + a0, default = missing(b2x2 + b1x + b0)
a4 = 0
a3 = 0
a2 = 1
a1 = 53
a0 = 11

P = X([x^2 + 58*x + 94, 68*x + 56])
Q = X([x^2 + 4*x + 24, 38*x + 65])
R = X([x^2 + 4*x + 63, 10*x + 26])

def extractpoints(MumDiv):
    """ input a divisor in mumford representation, output points in the divisor """
    ros = (MumDiv[0]).roots(multiplicities=True)
    #TBD: the general case
    if len(ros)==2:
        return (ros[0][0], MumDiv[1](ros[0][0])), (ros[1][0], MumDiv[1](ros[1][0]))
    elif len(ros)==1:
        return (ros[0][0], MumDiv[1](ros[0][0])), (ros[0][0], MumDiv[1](ros[0][0]))

for i in range(7):
    print i,   extractpoints(i*P),   extractpoints(i*Q),   extractpoints(i*R)

def threesome(XXX, YYY, ZZZ):
    """ create Y1Z1, X1Z2, X2Y2 ... """
    X1, X2 = extractpoints(XXX)
    Y1, Y2 = extractpoints(YYY)
    Z1, Z2 = extractpoints(ZZZ)
    return J(C(Y1)) + J(C(Z1)), J(C(X1)) + J(C(Z2)), J(C(X2)) + J(C(Y2)),     J(C(Y2)) + J(C(Z2)), J(C(X2)) + J(C(Z1)), J(C(X1)) + J(C(Y1)),     J(C(Y1)) - J(C(Z2)), - J(C(X1)) + J(C(Z1)), J(C(X2)) - J(C(Y1)),     J(C(Y2)) - J(C(Z1)), - J(C(X2)) + J(C(Z2)), J(C(X1)) - J(C(Y2))

Q1R1 = []
P1R2 = []
P2Q2 = []

Q2R2 = []
P2R1 = []
P1Q1 = []

Q1R2 = []
P1R1 = []
P2Q1 = []

Q2R1 = []
P2R2 = []
P1Q2 = []

for ii in range(7):
    Q1R1.append([])
    P1R2.append([])
    P2Q2.append([])
    Q2R2.append([])
    P2R1.append([])
    P1Q1.append([])
    Q1R2.append([])
    P1R1.append([])
    P2Q1.append([])
    Q2R1.append([])
    P2R2.append([])
    P1Q2.append([])
    for jj in range(7):
        Q1R1[ii].append(0)
        P1R2[ii].append(0)
        P2Q2[ii].append(0)
        Q2R2[ii].append(0)
        P2R1[ii].append(0)
        P1Q1[ii].append(0)
        Q1R2[ii].append(0)
        P1R1[ii].append(0)
        P2Q1[ii].append(0)
        Q2R1[ii].append(0)
        P2R2[ii].append(0)
        P1Q2[ii].append(0)

for i in range(1,7):
    for j in range(1,7):
        for k in range(1,7):
            Q1R1[j][k], P1R2[i][k], P2Q2[i][j],  Q2R2[j][k], P2R1[i][k], P1Q1[i][j],  Q1R2[j][k], P1R1[i][k], P2Q1[i][j],  Q2R1[j][k], P2R2[i][k], P1Q2[i][j] = threesome(i*P, j*Q, k*R)
            #print i, j, k, Q1R1[j][k], P1R2[i][k], P2Q2[i][j]

def Fg(Pin, Qin):
    u11, u12, v11, v12 = Pin[0][1], Pin[0][0], Pin[1][1], Pin[1][0]
    u21, u22, v21, v22 = Qin[0][1], Qin[0][0], Qin[1][1], Qin[1][0]
    return -(v11^2 - u11*u12 + u11^3 - a4*u11^2 + a3*u11) + (v21^2 - u21*u22 + u21^3 - a4*u21^2 + a3*u21) - u12*u21 + u11*u22

def Opairing(m, Pin, Qin):

    u11, u12, v11, v12 = Pin[0][1], Pin[0][0], Pin[1][1], Pin[1][0]
    u21, u22, v21, v22 = Qin[0][1], Qin[0][0], Qin[1][1], Qin[1][0]
    #print u21, u22, v21, v22

    W = [[FF(0) for x in range(3)] for y in range(m+5)]
    W[1][0] = FF(1)
    W[2][0] = FF( (-4*u12+6*u11^2+(-4*a4)*u11+2*a3)*v12 + 2*v11^3+((-8*u11+4*a4)*u12+2*u11^3 + (-2*a4)*u11^2+(2*a3)*u11-2*a2)*v11 )
    for i in range(3, m+2):
        b = i%2
        k = int((i-b)/2)  # i = 2k+b
        #print i, k, b
        if b==0:  # FIXME later
            W[i][0] = Fg( (k+1)*Pin, (k-1)*Pin )*W[k+1][0]^2*W[k-1][0]^2/W[2][0]
        else:
            W[i][0] = Fg( (k+1)*Pin, (k)*Pin )*W[k+1][0]^2*W[k][0]^2

    W[0][1] = FF(1)
    W[1][1] = FF(1)
    for i in range(2, m+2):
        b = i%2
        k = int((i-b)/2)  # i = 2k+b
        #print i, k, b
        if b==0:  # FIXME later
            W[i][1] = Fg( k*Pin + Qin, k*Pin )*W[k][1]^2*W[k][0]^2
            #print "===", i, k, b, k*Pin + Qin, Fg( k*Pin + Qin, k*Pin )
        else:
            W[i][1] = Fg( (k+1)*Pin + Qin, k*Pin )*W[k+1][1]^2*W[k][0]^2
            #print "===", i, k, b, (k+1)*Pin + Qin, Fg( (k+1)*Pin + Qin, k*Pin )

    for i in range(0):
        print i, W[i][0], W[i][1]
    return (W[8][1]/W[8][0])^((113^2-1)/7)
