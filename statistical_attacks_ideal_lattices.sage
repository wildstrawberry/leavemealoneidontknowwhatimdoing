#  experiments of Gentry-Szydlo averaging attack
#  the avergating term is alpha_i/alpha_j

n=4    # 2x the degree of cyclotomic polynomial
B=1    # bound for "small"
Big=14  # bound for "big"
Bound_secret = 1  # bound for "big"

print "secret Bound =", Bound_secret, "small Bound =", B, "Big bound =", B, "n=", n

#G = K.galois_group(); G

K.<x>=CyclotomicField(2*n)
OK = K.ring_of_integers()
R = ZZ["x"]
f = R.gen()**n + 1

#KK.<z>, L = K.subfield(x^2)
#print KK

#Kc=K.complex_conjugation()
#print Kc

def randVecOK(l):
    return Matrix(OK,[OK.random_element(x=-B,y=B+1) for i in range(l)])

def randMatOK(l):
    return Matrix(OK,[[OK.random_element(x=-B,y=B+1) for i in range(l)] for j in range(l)])

def toy_examples():
    """ inverse, norm ,.... in R or K """
    """ Z[z_n] is Euclidean iff n \in {1, 3, 4, 5, 7, 8, 9, 11, 12, 13, 15, 16, 17, 19, 20, 21, 24, 25, 27, 28, 32, 33, 35, 36, 40, 44, 45, 48, 60, 84 } """
    print "generating random elements f and g"
    f = OK.random_element(x=-B,y=B+1)
    g = OK.random_element(x=-B,y=B+1)
    #units = K.S_units([])
    #print units
    #print units[2]*units[2]

    invf, invg = 1/f, 1/g
    MF, MG = f.matrix(), g.matrix()
    charMF, charMG = MF.charpoly(), MG.charpoly()
    MFinv = 1/MF
    print "f=", f, "\n trace of f", f.trace(), "\n norm of f", f.norm(), "\n 1/f=", invf
    #print "fu=", f*units[1]
    print "g=", g, "\n norm of g", g.norm(), "\n 1/g=", invg
    print "MF=", MF, "\n 1/f=", invf.matrix()
    #print "MF eigenvalues", MF.eigenvalues()
    #print "MF eigenvectors", MF.eigenmatrix_right()
    #print "MG=", MG, "\n 1/g=", invg.matrix()

    fdg = f*invg
    Mfdg = fdg.matrix()
    print "f/g:", fdg
    print "g.norm() * f/g:", OK(g.norm()*fdg)
    I = OK.ideal([f])
    IJ = OK.ideal(OK(g.norm()*fdg))
    print "I = OK.ideal([f])", I.norm(), I
    print "I.is_principal()", I.is_principal()
    print IJ.norm(), IJ
    print IJ.is_principal()
#toy_examples()

def directavg():
    """ test whether 1/aj converges, without worrying about conjuagation.  """

    C = 100
    aj = [1/OK.random_element(x=-B,y=B+1) for i in range(C)]
    #aj = [float( (OK.random_element(x=-B,y=B+1)/OK.random_element(x=-B,y=B+1)).list()[0] ) for i in range(C)]
    #print aj
    summ = sum(aj)
    #print summ/C
    for i in range(4):
        print float(summ.list()[i])/C

def attack_rejection():

    Nsam = 1000
    threshold = 10
    sum_cinv = OK(0)
    sumzcinv = OK(0)
    sumz = OK(0)
    counter_success = 0
    s = OK.random_element(x = - B, y = B + 1)   # the msk
    print "the secret", s
    for i in range(Nsam):
        c = OK.random_element(x=-B,y=B+1)
        r = OK.random_element(x=-Big, y=Big+1)
        z = r+c*s
        if max(z.list()) <= threshold and min(z.list()) >= - threshold and (max(c.list())>0 or min(c.list())<0):
            #print max(z.list()), min(z.list())
            counter_success+=1
            sum_cinv += 1/c
            sumzcinv += z/c
            sumz += z
    print "number of success samples:", counter_success
    print "below are the averages over the things that are not rejected"
    for j in range(n):
        print "avg of z/c, the ", j, "th entry:", float(sumzcinv.list()[j])/counter_success
        print "avg of 1/c", float(sum_cinv.list()[j])/counter_success
        print "avg of z", float(sumz.list()[j])/counter_success

def attack_rejection_float():

    Nsam = 100000
    threshold = 9
    sum_cinv = [ float(0) for i in range(n) ]
    sumzcinv = [ float(0) for i in range(n) ]
    sumz = [ float(0) for i in range(n) ]
    counter_success = 0
    s = OK.random_element(x = - Bound_secret, y = Bound_secret + 1)   # the msk
    print "the secret", s
    for i in range(Nsam):
        c = OK.random_element(x=-B,y=B+1)
        r = OK.random_element(x=-Big, y=Big+1)
        z = r+c*s
        if max(z.list()) <= threshold and min(z.list()) >= - threshold and (max(c.list())>0 or min(c.list())<0):
            #print max(z.list()), min(z.list())
            counter_success+=1
            cinv = 1/c
            zcinv = z*cinv
            for j in range(n):
                sum_cinv[j]+= float( (cinv).list()[j] )
                sumzcinv[j]+= float( (zcinv).list()[j] )
                sumz[j]+= float( z.list()[j] )

    print "number of success samples:", counter_success
    print "below are the averages over the things that are not rejected"
    print "avg of z/c", sumzcinv
    print "avg of 1/c", sum_cinv
    print "avg of z", sumz

def averaging_GS():
    """ a\in R, suppose you get a1/a2, a1a3/a2a4, a1a3a5/a2a4a6, ... find a representation of a1*bar{a1} a.k.a. Norm_{K/KK}(a) where   """
    C=10000
    a2p =[K(OK.random_element(x=-B,y=B+1)) for i in range(C)]
    print "a0 = ", a2p[0]
    print "conj(a0)=", Kc(a2p[0])
    print "norm of a0", a2p[0].norm()
    a0ca0 = Kc(a2p[0])*a2p[0]
    print "a0*conj(a0)=", a0ca0, a0ca0.norm()
    print "inv of a0*conj(a0)=", 1/a0ca0, (1/a0ca0).norm()
    v=0
    for a in a2p:
        v=v+1/(Kc(a)*a)
    print float((v/C).list()[0])
    print float((v/C).list()[1])
    print float((v/C).list()[2])
    #print KK(v/C)/(a2p[0].norm(KK))

    
def averaging_GS2():
    """ difference from the above is, instead of assuming you have a/b1, a/b2, averaging over the inverse of them and get 1/a times 1/bar(a)   """
    C=1000
    a2p =[K(OK.random_element(x=-B,y=B+1)) for i in range(C)]
    a0 = K(OK.random_element(x=-B,y=B+1))
    print "a0 = ", a0
    print "conj(a0)=", Kc(a0)
    print "norm of a0", a0.norm()
    a0ca0 = Kc(a0)*a0
    print "a0*conj(a0)=", a0ca0, a0ca0.norm()
    print "inv of a0*conj(a0)=", 1/a0ca0, (1/a0ca0).norm()
    v=0
    for a in a2p:
        v=v+(Kc(a)*a)
    print float((v/C).list()[0])
    print float((v/C).list()[1])
    print float((v/C).list()[2])
    #print KK(v/C)/(a2p[0].norm(KK))

#for i in range(3):
#    averaging_GS2()
#    averaging_GS()


def GSexercises():
    """ test whether ai/aj converges  """

    C=400
    a2p =[K(OK.random_element(x=-B,y=B+1)) for i in range(C)]
    a0 = K(OK.random_element(x=-B,y=B+1))
    print "a0 = ", a0
    print "conj(a0)=", Kc(a0)
    print "norm of a0", a0.norm()
    a0ca0 = Kc(a0)*a0
    print "a0*conj(a0)=", a0ca0, a0ca0.norm()
    print "inv of a0*conj(a0)=", 1/a0ca0, (1/a0ca0).norm()
    v=0
    for a in a2p:
        for b in a2p:
            v=v+(Kc(a)*a)/(Kc(b)*b)
    print float((v/(C*C)).list()[0])
    print float((v/(C*C)).list()[1])
    print float((v/(C*C)).list()[2])
