#http://sage-doc.sis.uta.fi/reference/calculus/sage/calculus/transforms/dft.html

N = 8

J = list(range(N))
F = [ZZ(i**2+2*i) for i in J]
seq = IndexedSequence(F,J)
seq
seqhat = seq.dft()
lshat = seqhat.list()
print "fourier: ", lshat
lshatnorm = []
for i in lshat:
    lshatnorm.append(i.norm())
print "norm:", lshatnorm, lshatnorm[3]/lshatnorm[4]

#I = list(range(100))
#A = [ZZ(i^2)+1 for i in I]
#s = IndexedSequence(A,I)

seqnorm = IndexedSequence(lshatnorm[1:N],J[1:N])
P = seqnorm.plot()
show(P)
