#https://arxiv.org/pdf/1710.05147.pdf
def example1():
    q = 23
    FF = FiniteField(q)
    RR.<x> = PolynomialRing(FF)
    hy = x^5 + x^4 + 3*x^3 + 22*x^2 + 19*x
    C = HyperellipticCurve(hy)
    print C.frobenius_polynomial()
    K.<z> = NumberField(C.frobenius_polynomial())
    #print K.subfields()
    print K.maximal_totally_real_subfield()

# https://hal.archives-ouvertes.fr/hal-01520262v2/document
def example2():
    q = 2333
    FF = FiniteField(q)
    RR.<x> = PolynomialRing(FF)
    hy = 356*x^6 + 116*x^5 + 1589*x^4 + 986*x^3 + 178*x^2 + 1094*x + 1229
    C = HyperellipticCurve(hy)
    print C.frobenius_polynomial()
    K.<z> = NumberField(C.frobenius_polynomial())
    #print K.subfields()
    print K.maximal_totally_real_subfield()
example1()

# https://martindale.info/research/Thesis.pdf
# K.<z> = NumberField(x^4 - 605104*x^3 - 5215893977257194*x^2 - 1622371429548014920304*x + 7188537318834090069340399032601)
#K.<z> = NumberField(x^4 + 37*x^2 + 281)
#print K.maximal_totally_real_subfield()
# K.<z> = CyclotomicField(5)
#print K.subfields()
# where is Q(sqrt(5))?
