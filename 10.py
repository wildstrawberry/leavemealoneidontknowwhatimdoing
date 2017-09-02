# from http://doc.sagemath.org/html/en/reference/curves/sage/schemes/elliptic_curves/ell_curve_isogeny.html

E = EllipticCurve(GF(37), [0,0,0,1,8])
E.points()
E.count_points(1)
R.<x> = GF(37)[]
f = x^3 + x^2 + 28*x + 33  # kernel polynomial
phi = EllipticCurveIsogeny(E, f)
phi_hat = phi.dual()
print phi, phi.rational_maps()
print phi_hat, phi_hat.rational_maps()

P = E(3,1)
P.order()

Q = phi_hat(phi(P))
phi2 = EllipticCurveIsogeny(E, Q)
phi2_hat = phi2.dual()
print phi2_hat

print phi2_hat(phi2(Q))

phi7 = E.multiplication_by_m(7)
phi5 = E.multiplication_by_m(5)

print "multiply-by-7 map", phi7
print "multiply-by-5 map", phi5
