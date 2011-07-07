from numpy import cross, vdot
from sympy import Matrix, pi

#
# computes the i-th reciprocal lattice vector of a 3D lattice,
# where r is the array of the lattice base vectors.
#
def kvec3d(i, r):
    j = (i+1) % 3
    k = (i+2) % 3
    c = Matrix(cross(r[j].T, r[k].T)).T
    part3 = r[i].T * c
    return 2*pi/part3[0] * c
