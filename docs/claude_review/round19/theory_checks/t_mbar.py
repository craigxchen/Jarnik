from mbar import *
from math import factorial
import time
# tests: Mbar_{0,5}: D^2 = -1 for each boundary, psi^2 = 1; Mbar_{0,6}: psi_1^3 = 1, int psi^a multinomial
L5 = frozenset(range(5))
print("M05 D12^2 =", integrate(L5, [D(L5,{0,1}), D(L5,{0,1})]), " psi0^2 =", integrate(L5, [P(0), P(0)]), " psi0 psi1 =", integrate(L5,[P(0),P(1)]))
L6 = frozenset(range(6))
print("M06 psi0^3 =", integrate(L6,[P(0)]*3), " psi0^2 psi1 =", integrate(L6,[P(0),P(0),P(1)]), " (expect 1,3), psi0psi1psi2=", integrate(L6,[P(0),P(1),P(2)]), "(expect 6)")
L8 = frozenset(range(8))
t = time.time()
print("M08 psi0^2psi1psi2psi3 =", integrate(L8,[P(0),P(0),P(1),P(2),P(3)]), " expect", factorial(5)//2)
print("time", time.time()-t)
# D_T . psi^a vs formula
T = frozenset({0,1,2}); a = [1,0,0,2,1,0,0,0]  # sum over T = 1 = |T|-2, T^c sum 3 = |T^c|-2
cl = [D(L8,T)] + [P(i) for i in range(8) for _ in range(a[i])]
print("D_{012}.psi_0 psi_3^2 psi_4 =", integrate(L8, cl), " formula:", 1 * (factorial(3)//(2)))
