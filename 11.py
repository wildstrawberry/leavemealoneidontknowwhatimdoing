# from http://doc.sagemath.org/html/en/reference/curves/sage/schemes/elliptic_curves/ell_curve_isogeny.html
from sage.schemes.elliptic_curves.ell_curve_isogeny import compute_isogeny_starks, compute_sequence_of_maps

E = EllipticCurve(GF(11), [0,0,0,0,4])
print "j-inv of E: ", E, E.j_invariant(), E.count_points(1)
E2 = EllipticCurve(GF(11), [0,0,0,5,0])
print "j-inv of E2: ", E2, E2.j_invariant(), E2.count_points(1)
listpoints =  E.points()
for pos in listpoints:
    print pos, pos.order()
R.<x> = GF(11)[]
f = x^2 -16*x + 64
f.roots()
phi = EllipticCurveIsogeny(E, f)  # generate an isogeny from kernel polynomial
#phi_hat = phi.dual()
#print "the isogeny:", phi, "\n the kernel poly:", phi.kernel_polynomial(), "\n The rational maps ", phi.rational_maps()

#P = E(3,36)
#print P, P.order(), 7*P
