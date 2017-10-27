# from http://doc.sagemath.org/html/en/reference/curves/sage/schemes/elliptic_curves/ell_curve_isogeny.html

E = EllipticCurve(GF(11), [0,0,0,0,4])
print E, "j-inv of E: ", E.j_invariant(), "#(E)=", E.count_points(1)
E2 = EllipticCurve(GF(11), [0,0,0,5,0])
print E2, "j-inv of E2: ", E2.j_invariant(), "#(E2)=", E2.count_points(1)
R.<x> = GF(11)[]
kerf = x - 8
kerf.roots()
phi = EllipticCurveIsogeny(E, kerf)  # generate an isogeny from kernel polynomial
#phi_hat = phi.dual()
print "the isogeny:", phi, "\n the kernel poly:", phi.kernel_polynomial(), "\n The rational maps ", phi.rational_maps()

listpoints =  E.points()
for PP in listpoints:
    print PP, PP.order(), phi(PP), phi(PP).order()

def Kohel_bruteforce( psi ):
    """ generate the map (for K where char(K) is odd) directly from the kernel psi """
    return 1/psi^2

Kohel_bruteforce( kerf )
