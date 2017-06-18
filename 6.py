# Hyperelliptic nets

# Examples of curves
# E: y^2 = x^3 + 30x + 34 mod 631
# P = E((36, 60))
# Q = E((121, 387))

# q = 61
# FF = FiniteField(q)
# R.<x> = PolynomialRing(FF)
# hy = x**5 + 12

# q = 113
# FF = FiniteField(q)
# R.<x> = PolynomialRing(FF)
# hy = x^5 + x^2 + 53*x + 11
# P, Q, R: (x^2 + 4*x + 63, y + 10*x + 26) (x^2 + 4*x + 24, y + 38*x + 65) (x^2 + 58*x + 94, y + 68*x + 56)

q = 113
FF = FiniteField(q)
R.<x> = PolynomialRing(FF)
hy = x^5 + x^2 + 53*x + 11
C = HyperellipticCurve(hy)
Jac = C.jacobian()
X = Jac(FF)
print "X:", X   # fake set
ch = C.frobenius_polynomial()
print "Frobenius:", ch, " #Jac(C) = ", ch(1), "=", factor(ch(1))

# y2 + (b2x2 + b1x + b0)y = x5 + a4x4 + a3x3 + a2x2 + a1x + a0, default = missing(b2x2 + b1x + b0)
a4 = 0
a3 = 0
a2 = 1
a1 = 53
a0 = 11

P = X([x^2 + 4*x + 63, 10*x + 26])
Q = X([x^2 + 4*x + 24, 38*x + 65])
R = X([x^2 + 58*x + 94, 68*x + 56])

def Fg(Pin, Qin):
    u11, u12, v11, v12 = Pin[0][1], Pin[0][0], Pin[1][1], Pin[1][0]
    u21, u22, v21, v22 = Qin[0][1], Qin[0][0], Qin[1][1], Qin[1][0]
    return -(v11^2 - u11*u12 + u11^3 - a4*u11^2 + a3*u11) + (v21^2 - u21*u22 + u21^3 - a4*u21^2 + a3*u21) - u12*u21 + u11*u22

def Opairing(m, Pin, Qin, Rin):

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
    W[1][1][1] = FF(1)  #?
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

    return (W[7][1][1])**16
