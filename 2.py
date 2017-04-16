# lookslike 2x10x20
q = 13
TH = 52

FF = FiniteField(q)
R.<x> = PolynomialRing(FF)
# the text book curve f = x**5 + 1184*x**3 + 1846*x**2 + 956*x + 560
f = x**5 + 4*x**3 + (2401)*x
C = HyperellipticCurve(f)
n = C.count_points(1)
print "number of points:", n, "factor n:", factor(n[0])

J = C.jacobian()
X = J(FF)
print "X:", X   # fake set
print "Frobenius:",C.frobenius_polynomial()


def quicker_enumeratepoints():
    points = [] # list of all the points
    listorder = []  # the group structure
    for i in range (TH):
        listorder.append(0)

    for i in range(q):
        for j in range(q):
            a = x**2 + i*x + j
            for k in range(q):
                for l in range(q):
                    b = k*x + l  # TBD: compute square root
                    try:
                        D = X([a,b])
                        if D not in points:
                            points.append(D)
                        cz = 2
                        while (cz<TH):
                            czD = cz*D
                            if czD==D:
                                print "ijkl:",i,j,k,l
                                print cz-1, D  # (cz-1)*D == (1)
                                listorder[cz-1]+=1
                                break
                            else:
                                cz=cz+1
                                if czD not in points:
                                    points.append(czD)
                            if cz == TH:
                                print "order>",TH
                    except ValueError:
                        continue
    print listorder
    print points, len(points)


def bruteforce_enumeratepoints():
    counterJ=0  #this is a fake counter
    # try to find points on jacobian by enumeration
    listorder = []
    for i in range (TH):
        listorder.append(0)

    for i in range(q):
        for j in range(q):
            a = x**2 + i*x + j
            for k in range(q):
                for l in range(q):
                    b = k*x + l  # TBD: compute square root
                    try:
                        D = X([a,b])
                        counterJ=counterJ+1
                        cz = 2
                        while (cz<TH):
                            czD = cz*D
                            if czD==D:
                                print "ijkl:",i,j,k,l
                                print cz-1, D  # (cz-1)*D == (1)
                                listorder[cz-1]+=1
                                break
                            else:
                                cz=cz+1
    #                    if cz == TH:
     #                       print "order>",TH
                    except ValueError:
                        continue
    print counterJ, factor(counterJ)
    print listorder

quicker_enumeratepoints()
