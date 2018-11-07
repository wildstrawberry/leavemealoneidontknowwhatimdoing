TH = 2000  # threshold order

q = 1009
FF = FiniteField(q)
RR.<x> = PolynomialRing(FF)

# examples from https://hal.archives-ouvertes.fr/hal-01088933/document
# "Computing functions on Jacobians and their quotients"

hC = x*(x-1)*(x-2)*(x-3)*(x-85)
C = HyperellipticCurve(hC)
JC = C.jacobian()(FF)
charpolyC = C.frobenius_polynomial()
print "Frobenius:", charpolyC, " #Jac(C) = ", charpolyC(1), factor(charpolyC(1))

hD = x*(x-513)*(x-51)*(x-243)*(x-987)
D = HyperellipticCurve(hD)
JD = D.jacobian()(FF)
charpolyD = D.frobenius_polynomial()
print "Frobenius:", charpolyD, " #Jac(D) = ", charpolyD(1), factor(charpolyD(1))

SXN = 354*(x^5 + 647*x^4 + 931*x^3 + 597*x^2 + 73*x + 361)
PXN = 50*(x^5 + 262*x^4 + 812*x^3 + 770*x^2 + 868*x + 314)
XD = x^5 + 832*x^4 + 811*x^3 + 215*x^2 + 420*x

RYN = 304*(x^6 + 437*x^5 + 623*x^4 + 64*x^3 + 194*x^2 + 3*x + 511)
TYN = 678*(x^6 + 697*x^5 + 263*x^4 + 895*x^3 + 859*x^2 + 204*x + 130)
YD = x^8 + 239*x^7 + 983*x^6 + 800*x^5 + 214*x^4 + 489*x^3 + 191*x^2

print XD, XD.factor()
print YD, YD.factor()

T1 = x^2 + 247*x + 67
T2 = x^2 + 903*x + 350
print T1.factor(), T2.factor()

Y1 = JD([ x^2 - SXN(-510)/XD(-510)*x + PXN(-510)/XD(-510), (599-261*510)*(TYN(-510) + x*RYN(-510))/YD(-510)  ])
Y2 = JD([ x^2 - SXN(-746)/XD(-746)*x + PXN(-746)/XD(-746), (599-261*746)*(TYN(-746) + x*RYN(-746))/YD(-746)  ])
print Y1, Y2, Y1+Y2


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

def factor_torsion(n):
    """ return 1 if n has factors 3, 5 with multiplicity >3 """
    return (gcd(n, 343) == 343)

def my_mind_is_blowing():
    mind = []
    i, j, k, l = 0, 0, 0, 0
    while(len(mind)<40 or i==q):
        try:
            fff = x**5 + 0*x**4 + i*x**3 + j*x**2 + k*x + l
            i,j,k,l = iterrule(i,j,k,l)
            C = HyperellipticCurve(fff)
            J = C.jacobian()
            X = J(FF)
            ch = C.frobenius_polynomial()
            if factor_torsion(ch(1)):
                mind.append(C)
                print "C:", C, "Frobenius:", ch, " #J(C) = ", ch(1), factor(ch(1))
        except ValueError:
            continue
    return mind

#print my_mind_is_blowing()

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

#D = X([x,24])
#print D
