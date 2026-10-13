"""Exhaustive zonotope search for the sdc7 profile, T=1, L=3, target 10 (= rho): is zeta(A) <= 9?"""
import time
from zono import zeta_search
from model8 import profile
t0 = time.time()
z = zeta_search(8, 1, profile('sdc7'), 3, 10)
print("sdc7 T=1 L=3: zonotope with sum L_i >= 10 inside A:", z, "(%.0fs)" % (time.time() - t0))
