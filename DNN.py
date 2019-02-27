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

# N is batch size; D_in is input dimension;
# H is hidden dimension; D_out is output dimension.
N, H = 1000, [2, 8, 8, 2]
D_in, D_out = H[0], H[-1]
print D_in, D_out

ITER = 4000
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

learning_rate = 1e-4

#print x
#print y

for t in range(ITER+1):
    # Forward pass: compute predicted y
    h1 = nonlin(np.dot(x,w[0]))
#    h1 = x.dot(w[0])
#    h1_relu = np.maximum(h1, 0)
    h2 = nonlin(np.dot(h1,w[1]))
#    h2 = h1_relu.dot(w[1])
#    h2_relu = np.maximum(h2, 0)
    h3 = nonlin(np.dot(h2,w[2]))
#    y_pred = h2_relu.dot(w[2])

    # Compute and print loss
    loss = np.square(h3 - y).sum()

    h3_error = y - h3

    if (t<5) or (t%ITER_print == 0):
        print t, "total loss in the l2 norm:", loss, "; mean loss in the l1 norm:", np.mean(np.abs(h3_error))

    # in what direction is the target value?
    # were we really sure? if so, don't change too much.
    h3_delta = h3_error*nonlin(h3,deriv=True)

    # how much did each l2 value contribute to the l3 error (according to the weights)?
    h2_error = h3_delta.dot(w[2].T)
    h2_delta = h2_error*nonlin(h2,deriv=True)
    h1_error = h2_delta.dot(w[1].T)
    h1_delta = h1_error*nonlin(h1,deriv=True)

    w[2] += learning_rate * h2.T.dot(h3_delta)
    w[1] += learning_rate * h1.T.dot(h2_delta)
    w[0] += learning_rate * x.T.dot(h1_delta)

#print w[0], w[1]
