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

# the numerators and denominators of the isogeny
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

Y11 = JD([ x^2 - SXN(-510)/XD(-510)*x + PXN(-510)/XD(-510), (599-261*510)*(TYN(-510) + x*RYN(-510))/YD(-510)  ])
Y12 = JD([ x^2 - SXN(-746)/XD(-746)*x + PXN(-746)/XD(-746), (599-261*746)*(TYN(-746) + x*RYN(-746))/YD(-746)  ])
print "evalutaions on the roots of the kernel:", Y11, Y12, Y11+Y12
