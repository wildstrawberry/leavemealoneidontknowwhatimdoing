# folk examples from http://doc.sagemath.org/html/en/reference/curves/sage/schemes/hyperelliptic_curves/jacobian_generic.html

q = 13
TH = 13
listorder = []
for i in range (TH):
    listorder.append(0)

FF = FiniteField(q)
R.<x> = PolynomialRing(FF)
# the text book curve f = x**5 + 1184*x**3 + 1846*x**2 + 956*x + 560
f = x**5 + 4*x**3 + (2401)*x
C = HyperellipticCurve(f)
n = C.count_points(1)
print "number of points:", n, "factor n:", factor(n[0])

#rpt = C.rational_points()
#print rpt
J = C.jacobian()
X = J(FF)
print "X:",X
print "Frobenius:",C.frobenius_polynomial()

counterJ=0  #this is a fake counter
# try to find points on jacobian by enumeration
for i in range(q):
    for j in range(q):
        a = x**2 + i*x + j
        for k in range(q):
            for l in range(q):
                b = k*x + l  # TBD: compute square root
                try:
                    D = X([a,b])
                    counterJ=counterJ+1
                    cz = 2
                    while (cz<TH):
                        czD = cz*D
                        if czD==D:
                            print "ijkl:",i,j,k,l
                            print cz-1, D  # (cz-1)*D == (1)
                            listorder[cz-1]+=1
                            break
                        else:
                            cz=cz+1
#                    if cz == TH:
 #                       print "order>",TH
                except ValueError:
                    continue
print counterJ, factor(counterJ)
print listorder
