# HE pairing according to http://homes.esat.kuleuven.be/~fvercaut/papers/pairing2007.pdf

q = 13
# f = x**5 + 2*x**3 + (16)*x,   q=13, seems to have torsion point of rank 3, 
# (x^2 + x + 9, y + 5*x + 8) (x^2 + 12*x + 9, y + 12*x + 12)
# f = x**5 + 4*x**3 + (2401)*x, q=31, seems to be 12x30x120, not useful
# f = x**5 + 4*x**3 + (2401)*x, q=13, seems to be 2x10x20, not useful

TH = 10000  # threshold order

FF = FiniteField(q)
genus = 2
R.<x> = PolynomialRing(FF)
# the text book curve f = x**5 + 1184*x**3 + 1846*x**2 + 956*x + 560
f = x**5 + 2*x**3 + (16)*x    # chosen according to the special case mentioned in https://eprint.iacr.org/2011/604.pdf
C = HyperellipticCurve(f)
n = C.count_points(1)
print "number of points:", n, "factor n:", factor(n[0])

J = C.jacobian()
X = J(FF)
print "X:", X   # fake set
print "Frobenius:", C.frobenius_polynomial()

def Jacobian_order(D):
    """ input a divisor, output its order """
    ide = X([1,0])
    i=2
    while (i<TH):
        if i*D == ide:
            return i
        else:
            i=i+1

def Jacobian_enumeratepoints():
    """ Enumerate the points on jacobian by enumerate the coefficients """
    points = [] # list of all the points

    for i in range(q):
        for j in range(q):
            a = x**2 + i*x + j
            for k in range(q):
                for l in range(q):
                    b = k*x + l  # TBD: compute b by taking square root mod a
                    try:
                        D = X([a,b])
                        points.append(D)
                    except ValueError:
                        continue
    print "Number of points:", len(points)
    return points

#JJJ = Jacobian_enumeratepoints()

# q = 13, f = x**5 + 2*x**3 + (16)*x
# two divisors that seem to stay in free tate module of rank 2:
# 1 (x^2 + x + 9, y + 5*x + 8) (x^2 + 12*x + 9, y + 12*x + 12)
# 2 (x^2 + x + 9, y + 8*x + 5) (x^2 + 12*x + 9, y + x + 1)
# 3 (1) (1)

a1 = x^2 + x + 9
b1 = 8*x + 5
a2 = x^2 + 12*x + 9
b2 = x + 1

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
    print "h1:", h1
    h2 = 1
    h3 = 1
    s1 = c1*e1
    s2 = c1*e2
    s3 = c2
    u = u1*u2/(d*d)
    print "u,d", u, d
    v = (s1*u1*v2 + s2*u2*v1 + s3*(v1*v2 + f))/d
    print "v=", v
    v = v%u  # we have gcd, so they must work
    print "u,v:",u,v
    while (u.degree()>genus):
        ut = (f - v*v)/u
        ut = ut/ut.lc()
        vt = (-v) % ut
        h1 = (h1*(vE - v)) % uE
        h2 = (h2*ut) % uE
        if v.degree()>g:
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

print FD(3, a1, b1, a2, b2)
