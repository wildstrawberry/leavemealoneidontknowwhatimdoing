# http://doc.sagemath.org/html/en/reference/number_fields/sage/rings/number_field/class_group.html

#for d in [-4, -8..-400]:
for k in range(0):
    d = -4*(4*k+1)
    if is_fundamental_discriminant(d) and is_prime(4*k+1):
        FF = QuadraticField(d, 'a')
        h = FF.class_number()
        #if is_prime(h) and h>1:
        print k, d, h, is_prime(h)
        Cl = FF.class_group()
        print [c.representative_prime() for c in Cl]

def quadratic_character(Disc):
    """ chi_F(n) = Kron( D_F / n ) , D_F is negative for IQF """
    QR_sum = 0
    NQR_sum = 0
    QR = []
    NQR = []
    for x in range(-Disc):
        if kronecker(Disc, x)==1:
            QR_sum += x
            QR.append(x)
        elif kronecker(Disc, x)==-1:
            NQR_sum -= x
            NQR.append(x)
    print Disc
#    print QR, len(QR)
#    print NQR, len(NQR)
    print QR_sum, NQR_sum, (QR_sum+NQR_sum)/Disc

for k in range(0):
    quadratic_character(-4*(4*k+1))
    
def findniceclassgroups():
    for i in range(4):
      for j in range(4):
        for k in range(4):
          for l in range(4):
            Disc = (2^i)*(3^j)*(5^k)*(7^l)*11*13*17
            if is_fundamental_discriminant(-Disc):
                FF = QuadraticField(-Disc, 'a')
                h = FF.class_number()
                print "class number", h, i,j,k,l,"Disc", Disc
                #Cl = FF.class_group()
                #print [c.representative_prime() for c in Cl]

PRIMESET = Primes()
PSS = PRIMESET[1:400]  # prime subset
PSS[1:10]
PSS[11:20]

def findsmoothclassgroups(z):
    zbits = bin(z)[3:]  # bin(z)=0b...
    Disc = 4
    for i in range(z, z+4):
        Disc = Disc * PSS[i]
    if is_fundamental_discriminant(-Disc):
        FF = QuadraticField(-Disc, 'a')
        h = FF.class_number()
        print z, "class number", h, h.factor(), "Disc", Disc, Disc.factor()
        #Cl = FF.class_group()
        #print [c.representative_prime() for c in Cl]

for z in range(40, 100):
    findsmoothclassgroups(z)
    
"""
41 class number 21504 2^10 * 3 * 7 Disc 5780560756 2^2 * 191 * 193 * 197 * 199
42 class number 69760 2^7 * 5 * 109 Disc 6385855076 2^2 * 193 * 197 * 199 * 211
44 class number 53760 2^9 * 3 * 5 * 7 Disc 8502100676 2^2 * 199 * 211 * 223 * 227
46 class number 16960 2^6 * 5 * 53 Disc 10803938788 2^2 * 223 * 227 * 229 * 233
47 class number 62208 2^8 * 3^5 Disc 11579109284 2^2 * 227 * 229 * 233 * 239
49 class number 55168 2^7 * 431 Disc 13474249268 2^2 * 233 * 239 * 241 * 251
50 class number 52096 2^7 * 11 * 37 Disc 14862154772 2^2 * 239 * 241 * 251 * 257
51 class number 76496 2^4 * 7 * 683 Disc 16354588724 2^2 * 241 * 251 * 257 * 263
52 class number 49392 2^4 * 3^2 * 7^3 Disc 18254706916 2^2 * 251 * 257 * 263 * 269
53 class number 85120 2^7 * 5 * 7 * 19 Disc 19709265236 2^2 * 257 * 263 * 269 * 271
54 class number 40976 2^4 * 13 * 197 Disc 21243060196 2^2 * 263 * 269 * 271 * 277
56 class number 82656 2^5 * 3^2 * 7 * 41 Disc 23878212164 2^2 * 271 * 277 * 281 * 283
58 class number 35408 2^4 * 2213 Disc 28612693492 2^2 * 281 * 283 * 293 * 307
60 class number 31408 2^4 * 13 * 151 Disc 35024400772 2^2 * 293 * 307 * 311 * 313
61 class number 51584 2^7 * 13 * 31 Disc 37893293668 2^2 * 307 * 311 * 313 * 317
62 class number 43648 2^7 * 11 * 31 Disc 40855635844 2^2 * 311 * 313 * 317 * 331
64 class number 40960 2^13 * 5 Disc 49080233812 2^2 * 317 * 331 * 337 * 347
65 class number 250080 2^5 * 3 * 5 * 521 Disc 54034705364 2^2 * 331 * 337 * 347 * 349
67 class number 116224 2^9 * 227 Disc 61388079524 2^2 * 347 * 349 * 353 * 359
68 class number 74752 2^10 * 73 Disc 64926297364 2^2 * 349 * 353 * 359 * 367
69 class number 55488 2^6 * 3 * 17^2 Disc 69391143028 2^2 * 353 * 359 * 367 * 373
72 class number 108816 2^4 * 3 * 2267 Disc 84247380916 2^2 * 373 * 379 * 383 * 389
73 class number 71536 2^4 * 17 * 263 Disc 89668123924 2^2 * 379 * 383 * 389 * 397
75 class number 51200 2^11 * 5^2 Disc 101313607588 2^2 * 389 * 397 * 401 * 409
78 class number 90544 2^4 * 5659 Disc 124381757284 2^2 * 409 * 419 * 421 * 431
79 class number 92112 2^4 * 3 * 19 * 101 Disc 131680442308 2^2 * 419 * 421 * 431 * 433
80 class number 142400 2^6 * 5^2 * 89 Disc 137965904948 2^2 * 421 * 431 * 433 * 439
82 class number 112960 2^6 * 5 * 353 Disc 151238539636 2^2 * 433 * 439 * 443 * 449
83 class number 126848 2^7 * 991 Disc 159621276244 2^2 * 439 * 443 * 449 * 457
86 class number 75408 2^4 * 3 * 1571 Disc 182211166468 2^2 * 457 * 461 * 463 * 467
88 class number 87696 2^4 * 3^3 * 7 * 29 Disc 201754085332 2^2 * 463 * 467 * 479 * 487
89 class number 232176 2^4 * 3 * 7 * 691 Disc 213955196324 2^2 * 467 * 479 * 487 * 491
90 class number 99840 2^9 * 3 * 5 * 13 Disc 228615937828 2^2 * 479 * 487 * 491 * 499
91 class number 184128 2^6 * 3 * 7 * 137 Disc 240070598596 2^2 * 487 * 491 * 499 * 503
93 class number 242688 2^10 * 3 * 79 Disc 266246573732 2^2 * 499 * 503 * 509 * 521
94 class number 261472 2^5 * 8171 Disc 279052020164 2^2 * 503 * 509 * 521 * 523
96 class number 325120 2^9 * 5 * 127 Disc 322540306964 2^2 * 521 * 523 * 541 * 547
97 class number 247760 2^4 * 5 * 19 * 163 Disc 344827161188 2^2 * 523 * 541 * 547 * 557
98 class number 151744 2^6 * 2371 Disc 371200175428 2^2 * 541 * 547 * 557 * 563
99 class number 250624 2^8 * 11 * 89 Disc 390412014452 2^2 * 547 * 557 * 563 * 569
"""
