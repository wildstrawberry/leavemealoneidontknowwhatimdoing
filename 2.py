# lookslike 2x10x20
# HE pairing according to http://homes.esat.kuleuven.be/~fvercaut/papers/pairing2007.pdf

q = 13
TH = 52  # threshold order

FF = FiniteField(q)
R.<x> = PolynomialRing(FF)
# the text book curve f = x**5 + 1184*x**3 + 1846*x**2 + 956*x + 560
f = x**5 + 4*x**3 + (2401)*x    # chosen according to the special case mentioned in https://eprint.iacr.org/2011/604.pdf
C = HyperellipticCurve(f)
n = C.count_points(1)
print "number of points:", n, "factor n:", factor(n[0])

J = C.jacobian()
X = J(FF)
print "X:", X   # fake set
print "Frobenius:",C.frobenius_polynomial()


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
    listorder = []  # list the number of points of order i, indicating the group structure
    for i in range (TH):
        listorder.append(0)

    for i in range(q):
        for j in range(q):
            a = x**2 + i*x + j
            for k in range(q):
                for l in range(q):
                    b = k*x + l  # TBD: compute b by taking square root mod a
                    try:
                        D = X([a,b])
                        if (len(points)<40):
                            print "ijkl:",i,j,k,l
                            print D
                        if D not in points:
                            points.append(D)
                        cz = 2
                        while (cz<TH):
                            czD = cz*D
                            if czD==D:
                                listorder[cz-1]+=1
                                break
                            else:
                                cz=cz+1
                                if czD not in points:
                                    points.append(czD)
                            if cz == TH:
                                print "order>",TH
                    except ValueError:
                        continue
    print listorder
    print points, len(points)
    return points

# q = 13, f = x**5 + 4*x**3 + (2401)*x
# two divisors that seem to stay in free tate module of rank 2:
# (x^2 + 9*x + 1, y + 2*x + 7) (x^2 + 8*x + 3, y + x + 10)
# (x^2 + 4*x + 1, y + 10*x + 4) (x^2 + 9*x + 3, y + 1)
# (x^2 + 4*x + 1, y + 3*x + 9) (x^2 + 9*x + 3, y + 12)
# (x^2 + 9*x + 1, y + 11*x + 6) (x^2 + 8*x + 3, y + 12*x + 3)
# (1) (1)

a1 = x^2 + 9*x + 1
b1 = 11*x + 6
a2 = x^2 + 8*x + 3
b2 = 12*x + 3

D1 = X([a1,b1])
D2 = X([a2,b2])
ide = X([1,0])
print D1, D2, ide

E = D1+D2
print Jacobian_order(E)

def MillerH(u1, v1, u2, v2, uE, vE):
    """ Input: D1 = [u1, v1], D2 = [u2, v2], E = [uE, vE].
        Output: Reduced divisor ρ(D1 + D2) and evaluation h^norm_{D1,D2}(E), represented by [h˜1(x), h˜2(x),h3].  """
    d1, e1, e2 = xgcd(u1,u2)


MillerH(a1, b1, a2, b2, a1, b1)

def miller(m, P):
    """ miller's algorithm for hyperelliptic curve pairing  """
    m = bin(m)[3:]
    n = len(m)

#JJJ = Jacobian_enumeratepoints()
