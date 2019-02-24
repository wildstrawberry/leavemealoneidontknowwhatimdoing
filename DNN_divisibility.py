# -*- coding: utf-8 -*-
"""
https://pytorch.org/tutorials/beginner/examples_tensor/two_layer_net_numpy.html#sphx-glr-beginner-examples-tensor-two-layer-net-numpy-py
Training to learn odd or even.
"""
import numpy as np

# N is batch size; D_in is input dimension;
# H is hidden dimension; D_out is output dimension.
N, H = 500, [2, 10, 2]
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

learning_rate = 1e-6

#print x
#print y

for t in range(ITER+1):
    # Forward pass: compute predicted y
    h = x.dot(w[0])
    h_relu = np.maximum(h, 0)
    y_pred = h_relu.dot(w[1])

    # Compute and print loss
    loss = np.square(y_pred - y).sum()
    if (t<10) or (t%ITER_print == 0):
        print(t, loss)#, y_pred

    # Backprop to compute gradients of w1 and w2 with respect to loss
    grad_y_pred = 2.0 * (y_pred - y)
    grad_w1 = h_relu.T.dot(grad_y_pred)
    grad_h_relu = grad_y_pred.dot(w[1].T)
    grad_h = grad_h_relu.copy()
    grad_h[h < 0] = 0
    grad_w0 = x.T.dot(grad_h)

    # Update weights
    w[0] -= learning_rate * grad_w0
    w[1] -= learning_rate * grad_w1
