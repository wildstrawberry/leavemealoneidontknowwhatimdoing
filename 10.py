# from http://doc.sagemath.org/html/en/reference/curves/sage/schemes/elliptic_curves/ell_curve_isogeny.html
from sage.schemes.elliptic_curves.ell_curve_isogeny import compute_isogeny_starks, compute_sequence_of_maps

E = EllipticCurve(GF(37), [0,0,0,1,8])
print "j-inv of E: ", E.j_invariant()
listpoints =  E.points()
#for pos in listpoints:
#    print pos, pos.order()
print "Number of points on E: ",E.count_points(1)
R.<x> = GF(37)[]

for A in range(7,15):
    for B in range(3,15):
        try:
            E2 = EllipticCurve(GF(37), [0,0,0,A,B])
            print E2, "j-inv:", E2.j_invariant()
            for degree in range(4,8):
                fd = compute_isogeny_starks(E, E2, degree)
                print degree, fd, fd.roots()
        except (ArithmeticError)and(ValueError):
            print "Oops, a signular curve, or no isogeny of prescribed degree"
            continue


E2 = EllipticCurve(GF(1009), [0,0,0,1,3])
factor(1060)
print E2
print "Number of points on E2: ",E2.count_points(1)
R2.<x> = GF(1009)[]

g = x^5 + 270*x^4 + 289*x^3 + 659*x^2 +533*x + 399
g.roots()

phiE2 = EllipticCurveIsogeny(E2, g)  # generate an isogeny from kernel polynomial
print phiE2
