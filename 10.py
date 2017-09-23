# from http://doc.sagemath.org/html/en/reference/curves/sage/schemes/elliptic_curves/ell_curve_isogeny.html

E = EllipticCurve(GF(37), [0,0,0,1,8])
print "j-inv of E: ", E.j_invariant()
print "all the points on E: ",E.points()
print "Number of points on E: ",E.count_points(1)
R.<x> = GF(37)[]

f = x^3 + x^2 + 28*x + 33  # kernel polynomial
f.roots()
phi = EllipticCurveIsogeny(E, f)  # generate an isogeny from kernel polynomial
phi_hat = phi.dual()
print phi #, phi.rational_maps()
print phi_hat #, phi_hat.rational_maps()
(X, Y) = phi.rational_maps()
(Xhat, Yhat) = phi_hat.rational_maps()

P = E(25,9)
print "P and its order: ", P, P.order(), 2*P, 3*P

phi2 = EllipticCurveIsogeny(E, P)

phi3 = E.multiplication_by_m(3)
print phi3
#phi5 = E.multiplication_by_m(5)


E2 = EllipticCurve(GF(1009), [0,0,0,1,3])
factor(1060)
print E2
print "Number of points on E2: ",E2.count_points(1)
R2.<x> = GF(1009)[]

g = x^5 + 270*x^4 + 289*x^3 + 659*x^2 +533*x + 399
g.roots()

phiE2 = EllipticCurveIsogeny(E2, g)  # generate an isogeny from kernel polynomial
print phiE2
