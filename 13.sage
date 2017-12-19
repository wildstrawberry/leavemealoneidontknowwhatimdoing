p = next_prime(12)
PP = FiniteField(p)
RR.<x> = PolynomialRing(PP)
print p

AA = 1
BB = 3

def F1(x, y):
    """ eval f1(i,j) mod m  """
    return AA*(x^2 + y^2)

def F2(x, y):
    """ eval f2(i,j) mod m  """
    return BB*(x^3 + y^3) + x*y

for i in PP:
    for j in PP:
        if F1(i,j)==0:
            print "F1", i, j
        if F2(i,j)==0:
            print "F2", i, j
