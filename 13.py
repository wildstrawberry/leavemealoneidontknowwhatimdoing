# reference: http://doc.sagemath.org/html/en/reference/number_fields/sage/rings/number_field/class_group.html

for d in [-1, -2..-600]:
    if is_fundamental_discriminant(d):
        h = QuadraticField(d, 'a').class_number()
        if is_prime(h) and h>6:
            print d, h, is_prime(h)

K.<a> = NumberField(x^2+167)
K.class_number()
Cl = K.class_group()
[c.representative_prime() for c in Cl]
