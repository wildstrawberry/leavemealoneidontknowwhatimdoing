
q = 61

TH = 2000  # threshold order

FF = FiniteField(q)
R.<x> = PolynomialRing(FF)
f = x**5 + 13*x**4 + 2*x**3 + 4*x**2 + 11*x + 1
C = HyperellipticCurve(f)
n = C.count_points(1)
genus = C.genus()
print "genus:", genus, "number of points:", n, "factor n:", factor(n[0])

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

def my_mind_is_blowing():
    mind = []
    i, j, k, l = 0, 0, 0, 0
    while(len(mind)<50 or i==q):
        try:
            fff = x**5 + 0*x**4 + i*x**3 + j*x**2 + k*x + l
            i,j,k,l = iterrule(i,j,k,l)
            C = HyperellipticCurve(fff)
            print "C:", C
            mind.append(C)
            J = C.jacobian()
            X = J(FF)
            ch = C.frobenius_polynomial()
            print "Frobenius:", ch, " #J(C) = ", ch(1), factor(ch(1))
        except ValueError:
            continue
    return mind

print my_mind_is_blowing()

def Jacobian_find_torsion():
    """ Find group structure on jacobian by enumerating the coefficients """
    points = [] # list of all the points enumerated
    i, j, k, l = 0, 0, 0, 0
    while(len(points)<10 or i==q):
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

#print Jacobian_find_torsion()
