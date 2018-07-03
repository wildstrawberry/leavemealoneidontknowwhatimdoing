#http://sage-doc.sis.uta.fi/reference/calculus/sage/calculus/transforms/dft.html
import random
poww = 4
N = 2**poww


def biased_coin(prob):
    """ output 1 with probability prob  """
    ra = random.random()
    if ra<prob: return 1
    else: return 0

def int_to_bits_tlength(si, tlength):
    """ input an int si, turn into bits, and pad 00..0 ahead of it so that the length of the whole string is tlength """
    sbits = list(bin(si)[2:])
    sbits = [int(i) for i in sbits]
    if len(sbits)>= tlength:
        return sbits
    else:
        return [0 for i in range(tlength - len(sbits))] + sbits

J = list(range(N))
# F = [i%2 for i in range(N)]  # LSB
# F = [-1 for i in range(N/2)] + [1 for i in range(N/2)]  # half function
# F = [ZZ(i**2+2*i) for i in J]   # polynomial function

def LPN_secre():
    """ input a bit string as a list, pad 00..0 ahead of it so that the length of the whole string is tlength """
    sv = int_to_bits_tlength(4, poww)
    LPN = []
    for i in range(N):
        iv = int_to_bits_tlength(i, poww)
        y = 0
        for j in range(poww):
            y = y + sv[j]*iv[j]
        LPN.append((-1)**(y%2 + biased_coin(0.5) ))
    return LPN

F = LPN_secre()   # LPN with secret s without noise

seq = IndexedSequence(F,J)
seq
seqhat = seq.dft()
lshat = seqhat.list()
print "fourier: ", lshat

lshatnorm = []
for i in lshat:
    lshatnorm.append(i.norm())
print "norm:", lshatnorm #, lshatnorm[3]/lshatnorm[4]
normedlshatnorm = [float(i/sum(lshatnorm)) for i in lshatnorm]
print "normalized", normedlshatnorm

seqnorm = IndexedSequence(lshatnorm, J)
P = seqnorm.plot()
show(P)
