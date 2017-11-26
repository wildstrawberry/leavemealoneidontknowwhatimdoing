for d in [-1, -2..-300]:
    if is_fundamental_discriminant(d):
        h = QuadraticField(d, 'a').class_number()
        if is_prime(h) and h>10:
            print d, h, is_prime(h)
