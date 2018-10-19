qq = 113
FF = FiniteField(qq)
RR.<x> = PolynomialRing(FF)

ff0 = x^5 + x^2 + 53*x + 11
CC0 = HyperellipticCurve(ff0)
NOP_CC0 = CC0.count_points(1)
print "number of points on C0:", NOP_CC0, "factors:", factor(NOP_CC0[0])
Jac0 = CC0.jacobian()(FF)
ch0 = CC0.frobenius_polynomial()
print "Frobenius:", ch0, " #J(C) and its factorization:", ch0(1), factor(ch0(1))

D_0 = Jac0([x,89])
for i in range(10):
    print D_0*i
    
Lev = 4  # level of theta coordinates

def convert_mumford_to_theta():
    return
