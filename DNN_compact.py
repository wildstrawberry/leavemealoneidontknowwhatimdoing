# -*- coding: utf-8 -*-
"""
https://pytorch.org/tutorials/beginner/examples_tensor/two_layer_net_numpy.html#sphx-glr-beginner-examples-tensor-two-layer-net-numpy-py
https://iamtrask.github.io/2015/07/12/basic-python-network/
"""
import numpy as np

def nonlin(x,deriv=False):
    if(deriv==True):
        return x*(1-x)
    return 1/(1+np.exp(-x))

#np.random.seed(1)

# N is batch size; D_in is input dimension; D_out is output dimension; the middle of H are hidden dimensions;
N, H = 1000, [10, 8, 4, 2]
D_in, D_out, depthH = H[0], H[-1], len(H)
print D_in, D_out, depthH

ITER = 40000
ITER_print = ITER/10

# Create random input and output data
x = np.random.randint(2, size=(N, D_in))
y = [ [ xx[0], 1-xx[0] ] for xx in x ]
#y = [ [ xx[0], xx[1] ] for xx in x ]
#y = np.random.randint(2, size=(N, D_out))

w = []
# Randomly initialize weights
for i in range(len(H)-1):
    w.append( np.random.randn(H[i], H[i+1]) )

learning_rate = 1e-3

#print x
#print y

for t in range(ITER+1):
    # Forward pass: compute predicted y
    htemp = [ x ]
    herror = [ ]
    hdelta = [ ]
    for i in xrange(1,depthH):
        htemp.append( nonlin(np.dot( htemp[i-1] , w[i-1] )) )
    # Compute and print loss
    loss = np.square(htemp[-1] - y).sum()
    herror.append(y - htemp[-1])

    if (t<5) or (t%ITER_print == 0):
        print t, "total loss in the l2 norm:", loss, "; mean loss in the l1 norm:", np.mean(np.abs(herror[0]))

    for i in xrange(1,depthH):
        hdelta.append( herror[i-1]*nonlin(htemp[depthH-i],deriv=True) )
        herror.append( hdelta[i-1].dot(w[depthH-1-i].T) )
        w[depthH-1-i] += learning_rate * htemp[depthH-1-i].T.dot(hdelta[-1])

#print w[0], w[1]
#test group
def test():
    x = np.random.randint(2, size=(N, D_in))
    y = [ [ xx[0], 1-xx[0] ] for xx in x ]
    h_test = x
    for i in xrange(1,depthH):
        h_test =  nonlin(np.dot( h_test , w[i-1] )) 
    # Compute and print loss
    loss = np.square(h_test - y).sum()
    print "test => loss in the l2 norm:", loss, "; mean loss in the l1 norm:", np.mean(np.abs(y - h_test))
test()
