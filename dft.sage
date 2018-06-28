#http://sage-doc.sis.uta.fi/reference/calculus/sage/calculus/transforms/dft.html

J = list(range(8))
A = [ZZ(8-i) for i in J]
s = IndexedSequence(A,J)
s
s.dft()

P = s.plot()
show(P)
