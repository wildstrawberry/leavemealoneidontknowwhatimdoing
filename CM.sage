#https://arxiv.org/pdf/1710.05147.pdf

q = 23
FF = FiniteField(q)
RR.<x> = PolynomialRing(FF)
hy = x^5 + x^4 + 3*x^3 + 22*x^2 + 19*x
C = HyperellipticCurve(hy)
print C.frobenius_polynomial()
K.<z> = NumberField(C.frobenius_polynomial())
#print K.subfields()
print K.maximal_totally_real_subfield()

# I guess K_0 = x^2 - 32

# https://martindale.info/research/Thesis.pdf
# K.<z> = NumberField(x^4 - 605104*x^3 - 5215893977257194*x^2 - 1622371429548014920304*x + 7188537318834090069340399032601)
#K.<z> = NumberField(x^4 + 37*x^2 + 281)
#print K.maximal_totally_real_subfield()
# K.<z> = CyclotomicField(5)
#print K.subfields()
# where is Q(sqrt(5))?
