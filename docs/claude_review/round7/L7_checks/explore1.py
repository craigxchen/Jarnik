import sys
from fractions import Fraction as Fr
import qpoly as Q, cub13q as C, ratrec as RR
tau, s34, s35 = Fr(-2), Fr(3), Fr(5, 2)
f = lambda x: C.E(tau, s34, s35, x)
p, q = RR.reconstruct(f, maxdeg=40)
print('E as function of s45: deg num', Q.deg(p), 'deg den', Q.deg(q))
print('num roots rational:', Q.rational_roots(p))
