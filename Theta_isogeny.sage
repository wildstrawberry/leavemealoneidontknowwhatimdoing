# Isogeny computation using the Theta coordinates of Abelian Varieties
# References:
# [1] Isogeny https://arxiv.org/pdf/1001.2016.pdf
# [2] modular correspondences https://hal.archives-ouvertes.fr/hal-00426338v2/document
# [3] Convert Theta and Mumford http://www.normalesup.org/~robert/pro/publications/articles/niveau.pdf
# http://www.normalesup.org/~robert/pro/publications/slides/2011-02-Marseille_theta.pdf
# http://www.normalesup.org/~robert/pro/publications/slides/2010-07-Phd-Nancy.pdf

PP = 83^2
FF = FiniteField(PP)
RR.<x> = PolynomialRing(FF)

Ros = [0, FF(0), FF(1), FF(3), FF(15), FF(20)]
#LAM, MU, NU = FF(-35), FF(-42), FF(-128)   #Rosenhaim form
HF = 1
for i in Ros[1:]:
    HF = HF*(x-i)
CC = HyperellipticCurve(HF)
print CC #, "factors of f:", HF.factor()
JAC = CC.jacobian()(FF)
CH = CC.frobenius_polynomial() #charpoly of Frobenius
print "#J(C) and its factorization:", CH(1), factor(CH(1))#, "a point on the Jac", DD0

def Rosenhaim_to_theta_null():
    """ according to rosenhain.m """
    SqTheta = [ FF(0) for i in range(17) ]

    SqTheta[1] = 1

    SqTheta[3] = ((Ros[1]-Ros[4])*(Ros[3]-Ros[2])*(Ros[5]-Ros[2])
                    /(Ros[2]-Ros[4])/(Ros[3]-Ros[1])/(Ros[5]-Ros[1])).sqrt()
    SqTheta[4] = ((Ros[1]-Ros[4])*(Ros[3]-Ros[2])*(Ros[5]-Ros[4])
                    /(Ros[1]-Ros[3])/(Ros[4]-Ros[2])/(Ros[5]-Ros[3])).sqrt()
    SqTheta[5] = ((Ros[1]-Ros[2])*(Ros[1]-Ros[4])
                    /(Ros[1]-Ros[3])/(Ros[1]-Ros[5])).sqrt()
    SqTheta[7] = ((Ros[1]-Ros[4])*(Ros[3]-Ros[4])*(Ros[5]-Ros[2])
                    /(Ros[1]-Ros[5])/(Ros[3]-Ros[5])/(Ros[4]-Ros[2])).sqrt()

    if not (SqTheta[3]).is_square(): SqTheta[3]*=-1
    if not (SqTheta[4]).is_square(): SqTheta[4]*=-1
    if not (SqTheta[5]).is_square(): SqTheta[5]*=-1
    if not (SqTheta[7]).is_square(): SqTheta[7]*=-1

    SqTheta[6]=SqTheta[1]*SqTheta[4]/SqTheta[5]/Ros[5];
    SqTheta[8]=SqTheta[1]*SqTheta[7]/SqTheta[5]/Ros[3];
    SqTheta[2]=SqTheta[5]*SqTheta[6]/SqTheta[3]*(Ros[5]-1);
    SqTheta[9]=SqTheta[5]*SqTheta[8]/SqTheta[3]*(Ros[3]-1);
    SqTheta[10]=(SqTheta[1]*SqTheta[2]-SqTheta[3]*SqTheta[4])/SqTheta[8];
    return SqTheta
SQTHETA = Rosenhaim_to_theta_null()

print SQTHETA
((SQTHETA[2]-27)/29)^2982

U_P = (x-43)*(x-10)
#V_P = F.1^954*X + F.1^2518
#D_P = JAC([U_P, V_P])

def Jacobian_order(D):
    """ input a divisor, output its order """
    TH = 200  # a threshold for enumeration
    ide = JAC([1,0])
    i=2
    while (i<TH):
        if i*D == ide:
            return i
        else:
            i=i+1
