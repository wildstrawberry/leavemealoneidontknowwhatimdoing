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
    if is_fundamental_discriminant(-d):
        FF = QuadraticField(-d, 'a')
        h = FF.class_number()
        #if is_prime(h) and h>1:
        print d, h, is_prime(h)
        Cl = FF.class_group()
        for ele in Cl:
            print ele, ele.representative_prime()

insidecl(255)
