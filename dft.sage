#http://sage-doc.sis.uta.fi/reference/calculus/sage/calculus/transforms/dft.html
import random
poww = 6     # power
N = 2**poww  # number of inputs

def biased_coin(prob):
    """ output 1 with probability prob  """
    rand = random.random()
    if rand < prob: return 1
    else: return 0

def int_to_bits_tlength(si, tlength):
    """ input an int si, turn into bits, and pad 00..0 ahead of it so that the length of the whole string is tlength """
    sbits = list(bin(si)[2:])
    sbits = [int(i) for i in sbits]
    if len(sbits)>= tlength:
        return sbits
    else:
        return [0 for i in range(tlength - len(sbits))] + sbits

def SF_detector(Lin, delta):
    """ input the list of fourier coefficients, output the significant fourier bit, and anything that is above (delta/100)% to that guy """
    max_percentage = max(Lin)
    print "the max probability is: ", max_percentage, "the following indices are over", delta
    for i in range(len(Lin)):
        if Lin[i] > delta*max_percentage:
            print i, Lin[i]

J = list(range(N))
# F = [i%2 for i in range(N)]  # LSB
# F = [-1 for i in range(N/2)] + [1 for i in range(N/2)]  # half function
# F = [ZZ(i**2+2*i) for i in J]   # polynomial function

def LPN_secret(secret):
    """ input a secret < N, output the DFT of LPN """
    sv = int_to_bits_tlength(secret, poww)
    LPN = []
    for i in range(N):
        iv = int_to_bits_tlength(i, poww)
        y = 0
        for j in range(poww):
            y = y + sv[j]*iv[j]
        LPN.append((-1)**(y%2 + biased_coin(0.05) ))
    return LPN

F = LPN_secret(3)   # LPN with secret

seq = IndexedSequence(F,J)
seq
seqhat = seq.dft()
listhat = seqhat.list()
#print "fourier: ", lshat

listhatnorm = []
for i in listhat:
    listhatnorm.append( (i.norm())**(2/N) )  # hope to fix the norm to be the L2 norm; however Im not sure if (i.norm())**(2/N) is correct
#print "norm:", lshatnorm #, lshatnorm[3]/lshatnorm[4]
normalizednorm = [float(i/sum(listhatnorm)) for i in listhatnorm]
print "normalized", normalizednorm
SF_detector(normalizednorm, 0.8)

seqnorm = IndexedSequence(normalizednorm, J)
P = seqnorm.plot()
show(P)
