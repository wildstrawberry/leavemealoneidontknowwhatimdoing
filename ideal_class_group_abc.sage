# http://doc.sagemath.org/html/en/reference/number_fields/sage/rings/number_field/class_group.html

for d in [-1, -2..-6]:
    if is_fundamental_discriminant(d):
        h = QuadraticField(d, 'a').class_number()
        #if is_prime(h) and h>1:
        print d, h, is_prime(h)


def findniceclassgroups():
    for i in range(5):
      for j in range(5):
        for k in range(4):
            Disc = 2**i*3**j*5**k+1
            K.<a> = NumberField(x^2+Disc)
            CN = K.class_number()
            if is_prime(CN):
                print i,j,k,Disc, CN
                Cl = K.class_group()
                print [c.representative_prime() for c in Cl]

findniceclassgroups()
