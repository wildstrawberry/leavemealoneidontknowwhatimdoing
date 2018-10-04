PP = 83
FFP = FiniteField(PP)
RingFFP.<x> = FFP[]
print PP

db = HilbertClassPolynomialDatabase()
f251 = db[-251]

print RingFFP(f251).roots()

def j_to_A_B(j):
    """ from the j invariant, generate the coefficients A, B for Weierstrass form """
    return j/(48*(1728 - j)), j/(864*(1728-j))
#    return 3*j*(1728 - j), 2*j*(1728-j)^2
