

""" ---------- 3D Vector ----------- """
class Vec_3:

    def __init__(self, i, j, k):
        self.i = i
        self.j = j
        self.k = k

    def __add__(self, other):
        # add
        i = self.i + other.i
        j = self.j + other.j
        k = self.k + other.k
        return Vec_3(i, j, k)

    def __sub__(self, other):
        # subtract
        i = self.i - other.i
        j = self.j - other.j
        k = self.k - other.k
        return Vec_3(i, j, k)

    def __mul__(self, other):
        # scalar multiplication
        i = self.i * other
        j = self.j * other
        k = self.k * other
        return Vec_3(i, j, k)

    def __truediv__(self, other):
        # scalar multiplication
        i = self.i / other
        j = self.j / other
        k = self.k / other
        return Vec_3(i, j, k)

    def __floordiv__(self, other):
        # scalar multiplication
        i = self.i // other
        j = self.j // other
        k = self.k // other
        return Vec_3(i, j, k)

    def show(self):
        print(self.i, self.j, self.k)
        return
    
def dot(a, b):
    # dot product
    i = a.i * b.i
    j = a.j * b.j
    k = a.k * b.k
    return i + j + k

def cross(a, b):
    # cross product
    i = a.j * b.k - a.k * b.j
    j = a.k * b.i - a.i * b.k
    k = a.i * b.j - a.j * b.i
    return Vec_3(i, j, k)

def mag(a):
    # magnitude
    d = dot(a, a)
    return d**(.5)

def norm(a):
    m = mag(a)
    if m != 0:
        return a / m
    print("Normalized zero vector.")
    return 0/0

def mat_multi(M, v):
    """
    M example:
    M = [
    [a, b, c],
    [d, e, f],
    [g, h, i]
    ]
    """
    if len(M) != 3:
        print("Invalid Matrix.")
        return 0/0
    if len(M[0]) != 3 or len(M[1]) != 3 or len(M[2]) != 3:
        print("Invalid Matrix.")
        return 0/0
    r1 = Vec_3(M[0][0], M[0][1], M[0][2])
    i = dot(r1, v)
    r2 = Vec_3(M[1][0], M[1][1], M[1][2])
    j = dot(r2, v)
    r3 = Vec_3(M[2][0], M[2][1], M[2][2])
    k = dot(r3, v)
    return Vec_3(i, j, k)
