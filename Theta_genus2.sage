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
print CC0 #, FF0.factor()
JAC0 = CC0.jacobian()(FF)
CH0 = CC0.frobenius_polynomial() #charpoly of Frobenius
print "#J(C) and its factorization:", CH0(1), factor(CH0(1))#, "a point on the Jac", DD0

def Jacobian_order(D):
    """ input a divisor, output its order """
    ide = JAC0([1,0])
    i=2
    while (i<TH):
        if i*D == ide:
            return i
        else:
            i=i+1

def Rosenhaim_to_theta_null():
    """ Input lam, mu, nu, output the theta null. [3, p13]  """
    theta0 = FF(2)   # fix a theta 0 first
    theta1 = FF(97)
    theta2 = FF(70)  # this is the 4th number of Example 4.4
    theta3 = FF(44)  # this is the 3rd number of Example 4.4
    print "reference points, 4th power of 1/0, 2/0, 3/0", (theta1/theta0)^4, (theta2/theta0)^4, (theta3/theta0)^4

    theta_4_0_P4 = MU/(LAM*NU)
    theta_1_0_P4 = MU*(NU-1)*(LAM-1)/(LAM*NU*(MU-1))   #order of nu and lam doesn't matter for 1/0
    theta_8_0_P4 = MU*(NU-1)*(LAM-MU)*MU/(NU*(MU-1)*(LAM-NU))
    theta_2_0_P4 = MU*(LAM-1)*(NU-MU)/(LAM*(MU-1)*(NU-LAM))
    print "init 4th power of 1/0, 2/0, 4/0, 8/0", theta_1_0_P4, theta_2_0_P4, theta_4_0_P4, theta_8_0_P4

    theta_1_0_P2 = theta_1_0_P4.sqrt()
    theta_4_0_P2 = theta_4_0_P4.sqrt()
    theta_6_P2 = theta2^2/(NU*theta_4_0_P2) #need theta2
    #theta_12_P2 = theta2^2/(NU*theta_4_0_P2) #need theta2
    theta_3_P2 = (NU-1)*theta_6_P2*(theta_4_0_P2)/(theta_1_0_P2)
    print (theta_3_P2/(theta0^2))^2
    #wrap up stage, maybe not needed since we are at level 2
    #theta4 = (theta_4_0_P2.sqrt())*theta0
    return theta0, theta1, theta2, theta3

def Mumford_to_theta():
    """ Input mumford coordinate, output theta coordinate of level 2  """
    return 0

THETA0, THETA1, THETA2, THETA3 = Rosenhaim_to_theta_null()
print THETA0, THETA1, THETA2, THETA3

D_P = JAC0([x^2 + 53*x + 28, 29*x ])  # order = 3
D_Q = JAC0([x^2 + 32*x + 1, 52*x+41 ])  # order = 3
print Jacobian_order(D_P), Jacobian_order(D_Q)
