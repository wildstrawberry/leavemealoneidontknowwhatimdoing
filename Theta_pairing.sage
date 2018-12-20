# Pairing using the Theta coordinates of Abelian Varieties
# References:
# [0] Pairing http://www.normalesup.org/~robert/pro/publications/articles/pairings.pdf
# [1] Isogeny https://arxiv.org/pdf/1001.2016.pdf
# [2] modular correspondences https://hal.archives-ouvertes.fr/hal-00426338v2/document
# [3] Convert Theta and Mumford http://www.normalesup.org/~robert/pro/publications/articles/niveau.pdf
# http://www.normalesup.org/~robert/pro/publications/slides/2011-02-Marseille_theta.pdf
# http://www.normalesup.org/~robert/pro/publications/slides/2010-07-Phd-Nancy.pdf

TH = 200  # a threshold for enumeration
PP = 331
FF = FiniteField(PP)
RR.<x> = PolynomialRing(FF)

LAM, MU, NU = PP - 35, PP - 42, PP - 128   #Rosenhaim form
HF = x * (x - 1) * (x - LAM) * (x - MU) * (x - NU)    #x * (x + 35) * (x + 42) * (x + 128) * (x + 330)
CC = HyperellipticCurve(HF)
print CC #, "factors of f:", HF.factor()
JAC = CC.jacobian()(FF)
CH = CC.frobenius_polynomial() #charpoly of Frobenius
print "#J(C) and its factorization:", CH(1), factor(CH(1))#, "a point on the Jac", DD0

Lev = 2  # level of theta coordinates
#ThetaNull = [328, 213, 75, 1]  # why?
#P = [255, 89, 30, 1]  # need the convertion between theta and mumford to see what happened...
#DD0 = JAC([x^2 + 11, 15*x + 23])  # order = ?
#Q = []

def Jacobian_order(D):
    """ input a divisor, output its order """
    ide = JAC([1,0])
    i=2
    while (i<TH):
        if i*D == ide:
            return i
        else:
            i=i+1
