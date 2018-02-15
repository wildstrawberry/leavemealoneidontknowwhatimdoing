R = ZZ  # the base ring
F2 = FiniteField(2)

w=4

Iw = matrix.identity(w)
Pw = Matrix(ZZ,[[ZZ(((i+1)%w)==j) for i in range(w)] for j in range(w)])  # some non-identity matrix P of dimension w
Pinvw = Matrix(ZZ,[[ZZ(((i-1)%w)==j) for i in range(w)] for j in range(w)])  # P^{-1}

print "I =", Iw
print "P =", Pw
print "Pinverse =", Pinvw
print "P*P^{-1} =", Pw*Pinvw

Iw2 = matrix.identity(w)
Pw2 = Matrix(F2,[[ZZ(((i+1)%w)==j) for i in range(w)] for j in range(w)])  # some non-identity matrix P of dimension w
Pinvw2 = Matrix(F2,[[ZZ(((i-1)%w)==j) for i in range(w)] for j in range(w)])  # P^{-1}

print "I =", Iw2
print "P =", Pw2
print "Pinverse =", Pinvw2
print "P*P^{-1} =", Pw2*Pinvw2

print Iw2.charpoly(), Iw2.eigenvalues()
print Pw2.charpoly(), Pw2.eigenvalues()


#MMMM = block_matrix( [[Iw, Iw],[Iw, Iw]] )
#NNNN = block_matrix( [[Iw, Pinvw],[Pw, Iw]] )
