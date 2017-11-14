# References:
# http://doc.sagemath.org/html/en/reference/curves/sage/schemes/elliptic_curves/ell_curve_isogeny.html
# http://www.math.uwaterloo.ca/~mrubinst/modularpolynomials/phi_l.html
# https://math.mit.edu/~drew/ClassicalModPolys.html
from sage.schemes.elliptic_curves.ell_curve_isogeny import compute_isogeny_starks, compute_sequence_of_maps

FS = 37

RRR = GF(FS)

def printgroupmorphism(EE):
    listpoints =  EE.points()
    for PP in listpoints:
        print PP, PP.order(), phi(PP), phi(PP).order()

def modular3(X, Y):
    """ from Kohel's thesis, Psi function  """
    return (X+Y)^4 - X^3*Y^3 + 2232*X^2*Y^2*(X+Y)+36864000*(X+Y)^3 - 1069960*X*Y*(X + Y)^2 + 2590058000*X^2*Y^2 + 8900112384000*X*Y*( X + Y ) + 452984832000000*(X + Y)^2 - 771751936000000000*X*Y + 1855425871872000000000*(X + Y)

def modular5(X, Y):
    """ X_0(5) from BLS  """
    return 141359947154721358697753474691071362751004672000 + 53274330803424425450420160273356509151232000*(X+Y) - 264073457076620596259715790247978782949376*X*Y + 6692500042627997708487149415015068467200*(X^2+Y^2) +  36554736583949629295706472332656640000*(X^2*Y+X*Y^2) + 5110941777552418083110765199360000*(X^2*Y^2) +  280244777828439527804321565297868800*(X^3 + Y^3) -192457934618928299655108231168000*(X^3*Y+X*Y^3) + 26898488858380731577417728000*(X^3*Y^2+X^2*Y^3) - 441206965512914835246100*X^3*Y^3 + 1284733132841424456253440*(X^4 + Y^4) + 128541798906828816384000*(X^4*Y+X*Y^4) + 383083609779811215375*(X^4*Y^2+X^2*Y^4)+ 107878928185336800*(X^4*Y^3+X^3*Y^4) + 1665999364600*X^4*Y^4 + 1963211489280*(X^5 + Y^5) - 246683410950*(X^5*Y+X*Y^5) + 2028551200*(X^5*Y^2+X^2*Y^5) -4550940*(X^5*Y^3+X^3*Y^5) + 3720*(X^5*Y^4+X^4*Y^5) - X^5*Y^5 + (X^6+Y^6)

def check_connectivity():
    """ enumerate all the j variants, to see if they are connected """
    for j1 in xrange(0,FS):
        for j2 in xrange(0,FS):
            rc3 = modular3(RRR(j1), RRR(j2))
            if rc3==0:
                print j1, j2, "3"
            rc5 = modular5(RRR(j1), RRR(j2))
            if rc5==0:
                print j1, j2, "5"

#check_connectivity()

def finddots():
    """ find curves with specific j-inv's  """
    listcurve1 = []
    listcurve2 = []
    for A in xrange(0,FS):
        for B in xrange(0,FS):
            try:
                E1 = EllipticCurve(RRR, [0,0,0,A,B])
                if E1.j_invariant()==9 or E1.j_invariant()==18 or E1.j_invariant()==21 or E1.j_invariant()==34:
                    if E1.count_points(1)==36:
                        listcurve1.append(E1)
                    elif E1.count_points(1)==40:
                        listcurve2.append(E1)
                    #print "E1:", E1, "j(E1):", E1.j_invariant(), "#(E1)=", E1.count_points(1)
            except (ArithmeticError):
                continue
    return listcurve1, listcurve2

candidate36, candidate40 = finddots()
ll = len(candidate40)

# index 0 -> 42 -> 50 -> 45 -> 52 -> 4 -> 64
# j-inv 18->  9 -> 34 -> 21 -> 18 -> 9 -> 34

E1 = candidate40[11]
E2 = candidate40[56]

(isom1, isom2, E1pr, E2pr, ker_poly) = compute_sequence_of_maps(E1, E2, 5)
phi12 = EllipticCurveIsogeny(E1pr, ker_poly)
for P in E1.points():
    print P, P.order(), phi12(P), phi12(P).order()
#(isom1, isom2, E1pr, E2pr, ker_poly) = compute_sequence_of_maps(E1, E2, 5)
#print compute_sequence_of_maps(E2, E3, 5)
#print compute_sequence_of_maps(E3, E4, 5)
#print compute_sequence_of_maps(E4, E5, 5)


for i in xrange(0,0):
    for j in xrange(0,ll):
        E1 = candidate40[i]
        E2 = candidate40[j]
        if E1.j_invariant()!=E2.j_invariant():
            try:
                (isom1, isom2, E1pr, E2pr, ker_poly) = compute_sequence_of_maps(E1, E2, 5)
                #phi = EllipticCurveIsogeny(E, fd)
                if ker_poly.degree() ==2 and len(ker_poly.roots())==2:
                    phi = EllipticCurveIsogeny(E1pr, ker_poly)
                    print 5, i,j," : " ,E1.j_invariant(), E2.j_invariant(), ker_poly, ker_poly.roots()
                    print phi
            except (ValueError):
                #print "Oops, no isogenies with degree", degree
                continue
