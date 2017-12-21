# http://doc.sagemath.org/html/en/reference/number_fields/sage/rings/number_field/class_group.html

for d in [-3, -6..-30]:
    if is_fundamental_discriminant(d):
        h = QuadraticField(d, 'a').class_number()
        #if is_prime(h) and h>1:
        print d, h, is_prime(h)


def findniceclassgroups():
    for i in range(4):
      for j in range(4):
        for k in range(3):
          for l in range(1):
            Disc = (2^i)*(3^j)*(5^k)*(7^l)-1
            if is_fundamental_discriminant(-Disc):
                FF = QuadraticField(-Disc, 'a')
                print i, j, k, l, Disc, FF.class_number()
                Cl = FF.class_group()
                print [c.representative_prime() for c in Cl]

findniceclassgroups()
