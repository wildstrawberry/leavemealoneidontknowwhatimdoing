# Weil Pairing Example, paired by Cantor's algorithm
# Example 8.3 from Silverman
# E: y^2 = x^3 + 30x + 34 mod 631
# E(K) = 5x130

TH = 10000  # threshold order

q = 631
FF = FiniteField(q)
R.<x> = PolynomialRing(FF)
f = x**3 + 30*x + 34
C = HyperellipticCurve(f)
n = C.count_points(1)
genus = C.genus()
print "genus:", genus, "number of points:", n, "factor n:", factor(n[0])

J = C.jacobian()
X = J(FF)
print "X:", X   # fake set
print "Frobenius:", C.frobenius_polynomial()

# P = E((36, 60))
# Q = E((121, 387))
# n = 5

a1 = x - 36
b1 = 60
a2 = x - 121
b2 = 387

D1 = X([a1,b1])
D2 = X([a2,b2])
ide = X([1,0])
print "D1 = ", D1, "D2 = ", D2

def MillerH(u1, v1, u2, v2, uE, vE):
    """ Input: D1 = [u1, v1], D2 = [u2, v2], E = [uE, vE].
        Output: Reduced divisor ρ(D1 + D2) and evaluation h^norm_{D1,D2}(E), represented by [h˜1(x), h˜2(x),h3].
        H=0
    """
    d1, e1, e2 = xgcd(u1,u2)
    d, c1, c2 = xgcd(d1,v1+v2)
    h1 = d % uE
    print "h1, uE, d:", h1, uE, d
    h2 = 1
    h3 = 1
    s1 = c1*e1
    s2 = c1*e2
    s3 = c2
    dinv = d**(-1)
    u = R(u1*u2*dinv*dinv)
    #print "u,d", type(u), type(u1), type(d)
    v = R((s1*u1*v2 + s2*u2*v1 + s3*(v1*v2+f ))*dinv)%u
    print "u:",u,"v:",v
    while (u.degree()>genus):
        ut = R((f - v*v)/u)
        print "ut:", ut
        ut = ut/ut.lc()
        print "ut:", ut
        vt = (-v) % ut
        h1 = (h1*(vE - v)) % uE
        h2 = (h2*ut) % uE
        if v.degree()>genus:
            h3 = -(v.lc())*h3
        u = ut
        v = vt
    return u, v, h1, h2, h3

#MillerH(a1, b1, a2, b2, a1, b1)

def FD(m, u1, v1, u2, v2):
    """ miller's algorithm for hyperelliptic curve pairing
    Input: m, D1 = [u1, v1], D2 = [u2, v2].
    Output: Pairing value f_{m,D1}(eff(D2))
    """
    mbits = bin(m)[3:]  # bin(m)=0b...
    n = len(mbits)
    u = u1
    v = v1
    f = 1
    f1 = 1
    f2 = 1
    f3 = 1
    for i in range(n):
        f1 = (f1*f1) % u2
        f2 = (f2*f2) % u2
        f3 = f3*f3
        u, v, h1, h2, h3 = MillerH(u, v, u, v, u2, v2)
        f1 = (f1*h1) % u2
        f2 = (f2*h2) % u2
        f3 = f3*h3
        if int(mbits[i])==1:
            u, v, h1, h2, h3 = MillerH(u, v, u1, v1, u2, v2)
            f1 = (f1*h1) % u2
            f2 = (f2*h2) % u2
            f3 = f3*h3
    f = u2.resultant(f1)/( (f3**(u2.degree()))*u2.resultant(f2) )
    return f

print FD(5, a1, b1, a2, b2)
print FD(5, a2, b2, a1, b1)



def Jacobian_order(D):
    """ input a divisor, output its order """
    ide = X([1,0])
    i=2
    while (i<TH):
        if i*D == ide:
            return i
        else:
            i=i+1
