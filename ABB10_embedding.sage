import random
R = ZZ  # the base ring
F2 = FiniteField(2)

def ABB10(u0, u1, u2, u3):
    """ the Full-rank-difference mapping in ABB10 """
    M = Matrix(ZZ, [[u0, u1, u2, u3], [u3, u0-u3, u1, u2], [u2, u3-u2, u0-u3, u1], [u1, u2-u1, u3-u2, u0-u3]] )
    return M

def diff_ABB10(B):
    """ experiments with ABB10 difference  """
    M1 = ABB10(1, random.randint(-B, B), random.randint(-B, B), 0)
    M2 = ABB10(1, random.randint(-B, B), random.randint(-B, B), 0)
    print M1, M1.det()
    print M2, M2.det()
    print M1-M2, (M1-M2).det()

diff_ABB10(10)
