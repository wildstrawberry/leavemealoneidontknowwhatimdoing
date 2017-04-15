# folk examples from http://doc.sagemath.org/html/en/reference/curves/sage/schemes/hyperelliptic_curves/jacobian_generic.html

q = 7
FF = FiniteField(q)
R.<x> = PolynomialRing(FF)
f = x**5 + 1184*x**3 + 1846*x**2 + 956*x + 560
C = HyperellipticCurve(f)
n = C.count_points(1)
print "number of points:", n, "factor n:", factor(n[0])

#rpt = C.rational_points()
#print rpt
J = C.jacobian()
X = J(FF)
print "X:",X

counterJ=0
# try to find points on jacobian by enumeration
for i in range(q):
    for j in range(q):
        a = x**2 + i*x + j
        for k in range(q):
            for l in range(q):
                b = k*x + l
                try:
                    D = X([a,b])
                    print "ijkl:",i,j,k,l
                    print "D:", D
                    counterJ=counterJ+1
                except ValueError:
                    continue
print counterJ
