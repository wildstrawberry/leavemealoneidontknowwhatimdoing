%martin's proposal

n = 256
P.<x> = ZZ["x"]
phi = x^n + 1
f = P.random_element(x=-10, y=11, degree=n-1)
g = P.random_element(x=-10, y=11, degree=n-1)
h = (f * g.change_ring(QQ).inverse_mod(phi)) % phi
q = next_prime(2^64) # pick anything

h_ = h.change_ring(GF(q))
H = matrix([((x^i*h_)%phi).coefficients() for i in range(n)])
B = (identity_matrix(ZZ, n)).augment(H)
B = B.stack(matrix(ZZ, n, n).augment(q*identity_matrix(ZZ, n)))
B = B.LLL()

g_ = P(B[0][:n].list())  # some rotation of g
f_ = (g_*h) % phi        # some rotation of f
print(f)
print(f_)
print()
print(g)
print(g_)
