# HE pairing according to http://homes.esat.kuleuven.be/~fvercaut/papers/pairing2007.pdf

q = 61

TH = 2000  # threshold order

FF = FiniteField(q)
R.<x> = PolynomialRing(FF)
hy = x**5 + 0*x**4 + 0*x**3 + 0*x**2 + 0*x + 12
# x^5 + x^2 + 53*x + 11, q=113

C = HyperellipticCurve(hy)
genus = C.genus()

J = C.jacobian()
X = J(FF)
print "X:", X   # fake set
ch = C.frobenius_polynomial()
print "Frobenius:", ch, " #J(C) = ", ch(1), factor(ch(1))

def Jacobian_order(D):
    """ input a divisor, output its order """
    ide = X([1,0])
    i=2
    while (i<TH):
        if i*D == ide:
            return i
        else:
            i=i+1

def iterrule(i,j,k,l):
    l+=1
    if l==q:
        l=0
        k+=1
        if k==q:
            k=0
            j+=1
            if j==q:
                j=0
                i+=1
    return i,j,k,l

def Jacobian_find_torsion():
    """ Find group structure on jacobian by enumerating the coefficients """
    points = [] # list of all the points enumerated
    i, j, k, l = 0, 0, 0, 0
    while(len(points)<5 or i==q):
        a = x**2 + i*x + j
        b = k*x + l  # TBD: compute b by taking square root mod a
        i,j,k,l = iterrule(i,j,k,l)
        try:
            D = X([a,b])
            points.append(D)
            print Jacobian_order(D), a, b
            if gcd(5,Jacobian_order(D))==5:
                print int(Jacobian_order(D)/5)*D
        except ValueError:
            continue
    #print "Number of points:", len(points)
    return points

print Jacobian_find_torsion()


a1 = x^2 + 11*x + 16
b1 = x + 36
a2 = x^2 + 11*x + 16
b2 = 60*x + 25
a3 = x^2 + 38*x + 15
b3 = 11
a4 = x^2 + 2*x + 11
b4 = 21*x + 40

D1 = X([a1,  b1])
D2 = X([a2,  b2])
D3 = X([a3,  b3])
D4 = X([a4,  b4])
#print "D1, D2, D3, and their orders:", D1, Jacobian_order(D1), D2, Jacobian_order(D2),  D3, Jacobian_order(D3),  D4, Jacobian_order(D4)
for i in range(5):
    print i, i*D1, i*D2, i*D3, i*D4

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
    v = R((s1*u1*v2 + s2*u2*v1 + s3*(v1*v2+hy ))*dinv)%u
    print "u:",u,"v:",v,"s3:",s3
    while (u.degree()>genus):
        ut = R((hy - v*v)/u)
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
        print "u:", u, "v:", v, "h1:", h1, "h2:", h2, "h3:", h3
        f1 = (f1*h1) % u2
        f2 = (f2*h2) % u2
        f3 = f3*h3
        if int(mbits[i])==1:
            u, v, h1, h2, h3 = MillerH(u, v, u1, v1, u2, v2)
            f1 = (f1*h1) % u2
            f2 = (f2*h2) % u2
            f3 = f3*h3
    print "u2:", u2, "f1:", f1, "f3:", f3, "f2:", f2, u2.resultant(f2), u2.resultant(f1)  # sometimes both  u2.resultant(f2) and u2.resultant(f1) are zero ...
    f = u2.resultant(f1)/( (f3**(u2.degree()))*u2.resultant(f2) )
    return f

num = FD(5, a1, b1, a4, b4)
den = FD(5, a4, b4, a1, b1)

print "### num, den, e = -num/den:", num, den, -num/den

num2 = FD(5, a2, b2, a4, b4)
print num2

den2 = FD(5, a4, b4, a2, b2)
print den2, -num2/den2
