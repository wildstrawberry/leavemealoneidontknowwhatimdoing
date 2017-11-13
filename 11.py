# from http://doc.sagemath.org/html/en/reference/curves/sage/schemes/elliptic_curves/ell_curve_isogeny.html

RRR = GF(11)

E = EllipticCurve(RRR, [0,0,0,0,4])
print E, "j-inv of E: ", E.j_invariant(), "#(E)=", E.count_points(1)
E2 = EllipticCurve(RRR, [0,0,0,5,0])
print E2, "j-inv of E2: ", E2.j_invariant(), "#(E2)=", E2.count_points(1)
R.<x> = RRR[]
kerf = x - 8
kerf.roots()
phi = EllipticCurveIsogeny(E, kerf)  # generate an isogeny from kernel polynomial
#phi_hat = phi.dual()
print "the isogeny:", phi, "\n the kernel poly:", phi.kernel_polynomial(), "\n The rational maps ", phi.rational_maps()

listpoints =  E.points()
for PP in listpoints:
    print PP, PP.order(), phi(PP), phi(PP).order()
    
def modular3(X, Y):
    return (X+Y)^4 - X^3*Y^3 + 2232*X^2*Y^2*(X+Y)+36864000*(X+Y)^3 - 1069960*X*Y*(X + Y)^2 + 2590058000*X^2*Y^2 + 8900112384000*X*Y*( X + Y ) + 452984832000000*(X + Y)^2 - 771751936000000000*X*Y + 1855425871872000000000*(X + Y)

modular3(RRR(0),RRR(1))
