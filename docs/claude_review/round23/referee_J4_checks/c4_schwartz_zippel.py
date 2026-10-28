# Referee J4, check 4 (exact): J4 step 6.
# (a) The Q-structure of T_5: a class [z] has Q-coordinates w_j = z_j/z_1 = (x_j, y_j) in Q^2 with x^2+y^2 = 1.
#     A closed Z defined over Q by polynomials in (x_j, y_j) pulls back, via u_j = x_j + i y_j and the cone
#     (C^*)^5 -> (C^*)^5/C^*, to a homogeneous F(u_1..u_5) (over Q(i)), and [z] in Z  <=>  F(z_1,..,z_5) = 0
#     with z_j read as complex numbers.  Tested on three Z's against every ordered tuple of distinct points.
# (b) Schwartz-Zippel: #{t in Sigma^5 : F(t) = 0} <= deg F * M^4, Sigma = M points of one circle.
import sys
from itertools import product, permutations
from fractions import Fraction as Fr
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from gz import *

def wcoords(t):  # Q-coordinates of w_j = z_j / z_1
    z1 = t[0]; n1 = norm(z1); out = []
    for z in t:
        c = gmul(z, gconj(z1)); out.append((Fr(c[0], n1), Fr(c[1], n1)))
    return out
# Z defined in Q-coordinates
def Z1(w):  # w2*w3 = w4*w5 in the torus group law
    (x2, y2), (x3, y3), (x4, y4), (x5, y5) = w[1:]
    return x2*x3-y2*y3 == x4*x5-y4*y5 and x2*y3+x3*y2 == x4*y5+x5*y4
def Z2(w): return w[1][0] == w[2][0]                       # Re w2 = Re w3
def Z3(w): return sum(x for x, y in w[1:]) == 0            # x2+x3+x4+x5 = 0
# the homogeneous F in u (Gaussian-integer arithmetic, exact)
def mul(*a):
    r = (1, 0)
    for x in a: r = gmul(r, x)
    return r
def sub(a, b): return (a[0]-b[0], a[1]-b[1])
def add(*a): return (sum(x[0] for x in a), sum(x[1] for x in a))
def F1(u): return sub(mul(u[1], u[2]), mul(u[3], u[4]))                                  # deg 2
def F2(u):  # (u2/u1+u1/u2) = (u3/u1+u1/u3), times u1 u2 u3
    u1, u2, u3 = u[0], u[1], u[2]
    return sub(add(mul(u2, u2, u3), mul(u1, u1, u3)), add(mul(u3, u3, u2), mul(u1, u1, u2)))  # deg 3
def F3(u):  # sum_{j>=2} (u_j/u1 + u1/u_j) = 0, times u1 * u2 u3 u4 u5
    terms = []
    for j in range(1, 5):
        others = [u[l] for l in range(1, 5) if l != j]
        terms.append(mul(u[j], u[j], *others)); terms.append(mul(u[0], u[0], *others))
    return add(*terms)                                                                         # deg 5
tests = [("w2w3=w4w5", Z1, F1, 2), ("Re w2=Re w3", Z2, F2, 3), ("x2+x3+x4+x5=0", Z3, F3, 5)]
for N, M in ((1176852625, 14), (5**2*13**2*17*29, 16), (2*65**2*17, 16)):
    pts = lattice_points(N)
    Sig = pts[:M]
    for name, Zq, F, D in tests:
        agree = True; zeros_all = 0; zeros_distinct = 0
        for t in product(Sig, repeat=5):
            fz = F(t) == (0, 0)
            zeros_all += fz
            if len(set(t)) == 5:
                zeros_distinct += fz
                if Zq(wcoords(t)) != fz: agree = False
        print("N=%d M=%d Z=%-16s [z] in Z <=> F(z)=0 on all distinct ordered tuples: %s;"
              " zeros in Sigma^5: %d (distinct: %d) <= D*M^4 = %d : %s"
              % (N, M, name, agree, zeros_all, zeros_distinct, D*M**4, zeros_all <= D*M**4))
