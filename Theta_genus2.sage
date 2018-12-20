# Theta coordinates of Abelian Varieties
# References:
# [1] Isogeny https://arxiv.org/pdf/1001.2016.pdf
# [2] modular correspondences https://hal.archives-ouvertes.fr/hal-00426338v2/document
# [3] Convert Theta and Mumford http://www.normalesup.org/~robert/pro/publications/articles/niveau.pdf
# http://www.normalesup.org/~robert/pro/publications/slides/2011-02-Marseille_theta.pdf
# http://www.normalesup.org/~robert/pro/publications/slides/2010-07-Phd-Nancy.pdf

TH = 200  # a threshold for enumeration
QQ = 109
FF = FiniteField(QQ)
RR.<x> = PolynomialRing(FF)

LAM, MU, NU = FF(-33), FF(-61), FF(-98)  # Example 4.4 of [3]
FF0 = x*(x-1)*(x-LAM)*(x-MU)*(x-NU)   #Rosenhaim form
#FF0 = x^5 + 82*x^4 + 24*x^3 + 95*x^2 + 16*x  mod 109
CC0 = HyperellipticCurve(FF0)
print CC0, FF0.factor()
JAC0 = CC0.jacobian()(FF)
CH0 = CC0.frobenius_polynomial() #charpoly of Frobenius
#DD0 = JAC0([x^2 + 11, 15*x + 23])  # order = 12
print "#J(C) and its factorization:", CH0(1), factor(CH0(1))#, "a point on the Jac", DD0

Lev = 2  # level of theta coordinates

def convert_Rosenhaim_to_theta_null():
    """ Input lam, mu, nu, output the theta null. [3, p13]  """
    theta0 = FF(2)
    theta1 = FF(97)
    theta2 = FF(70)  # this is the 4th number of Example 4.4
    theta3 = FF(44)  # this is the 3rd number of Example 4.4
    print "reference group", (theta1/theta0)^4, (theta2/theta0)^4, (theta3/theta0)^4

    theta_4_0_P4 = MU/(LAM*NU)
    theta_1_0_P4 = MU*(NU-1)*(LAM-1)/(LAM*NU*(MU-1))   #order of nu and lam doesn't matter for 1/0
    theta_8_0_P4 = MU*(NU-1)*(LAM-MU)*MU/(NU*(MU-1)*(LAM-NU))
    theta_2_0_P4 = MU*(LAM-1)*(NU-MU)/(LAM*(MU-1)*(NU-LAM))
    print "init power of 4 theta ratios", theta_1_0_P4, theta_2_0_P4, theta_4_0_P4, theta_8_0_P4
    theta_1_0_P2 = theta_1_0_P4.sqrt()
    theta_4_0_P2 = theta_4_0_P4.sqrt()
    theta_6_P2 = theta2^2/(NU*theta_4_0_P2)
    theta_3_P2 = (NU-1)*theta_6_P2*(theta_4_0_P2)/(theta_1_0_P2)
    print theta_4_0_P2, theta_3_P2, (theta_3_P2/(theta0^2))^2
    return theta0, theta1, theta2, theta3

THETA0, THETA1, THETA2, THETA3 = convert_Rosenhaim_to_theta_null()
print THETA0, THETA1, THETA2, THETA3

def Jacobian_order(D):
    """ input a divisor, output its order """
    ide = JAC0([1,0])
    i=2
    while (i<TH):
        if i*D == ide:
            return i
        else:
            i=i+1

def iterrule(i,j,k,l):
    l+=1
    if l==QQ:
        l=0
        k+=1
        if k==QQ:
            k=0
            j+=1
            if j==QQ:
                j=0
                i+=1
    return i,j,k,l

def Jacobian_find_torsion(Num):
    """ Find group structure on jacobian by enumerating the coefficients """
    points = [] # list of all the points enumerated
    i, j, k, l = 0, 0, 0, 0
    while(len(points)<Num or i==QQ):
        a = x**2 + i*x + j
        b = k*x + l  # TBD: compute b by taking square root mod a
        i,j,k,l = iterrule(i,j,k,l)
        try:
            D = JAC0([a,b])
            points.append(D)
            print Jacobian_order(D), a, b
            #if gcd(11,Jacobian_order(D))==1:
            #    print int(Jacobian_order(D)/5)*D
        except ValueError:
            continue
    #print "Number of points:", len(points)
    return points
#Jacobian_find_torsion(10)
