db = HilbertClassPolynomialDatabase()

#for d in [-4, -8..-400]:

for k in range(251,270):
    d = -k
    if ( d%4 == 0 )or( d%4 == 1 ): #is_fundamental_discriminant(d): # and is_prime(4*k+1):
        FF = QuadraticField(d, 'a')
        h = FF.class_number()
        #if is_prime(h) and h>1:
        print d, is_fundamental_discriminant(d), h, db[d].degree()
        #Cl = FF.class_group()
        #print [c.representative_prime() for c in Cl]

def insidecl(d):
    Hil = db[-d]
    QD = QuadraticField(-d, 'a')
    h = QD.class_number()
    #if is_prime(h) and h>1:
    print "D = ",d, "class number:", h
    Cl = QD.class_group()
    for ele in Cl:
        print ele, ele.representative_prime()
    for pp in range(80, 90):
        if is_prime(pp):
            FFP = FiniteField(pp)
            PolyP.<x> = PolynomialRing(FFP)
            print pp, PolyP(Hil).roots()
insidecl(251)
