DIM = 3

M = []

for i in xrange(0,DIM):
    ROW = []
    for j in xrange(0,DIM):
        ROW.append(i*j-2*j+4)
    M.append( ROW )

print M

def grabmat(tens, d, jstar):
    """ cut a submatrix out of a whole by eliminating the 0-th row and jstar-th column"""
    mat = []
    for i in xrange(1,d):
        rowmat = []
        for j in xrange(0,d):
            if j<>jstar:
                rowmat.append(tens[i][j])
        mat.append(rowmat)
    return mat

def hyperdet(tens, d):
    """ input a 2(maybe more)-tensor, output the hyperdeterminant of dim d """

    if d==1:
        return tens[0][0]
    else:
        det = 0
        for i in xrange(0, d):
            subM = grabmat(tens, d, i)
            det = det + (-1)^i * tens[0][i] * hyperdet(subM, d-1)
        return det

print hyperdet(M, DIM)
