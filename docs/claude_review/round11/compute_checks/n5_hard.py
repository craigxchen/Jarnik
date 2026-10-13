"""The one undecided n = 5 instance: profile ((0,2),(0,3),(1,),(2,),(2,4),(3,4),(4,)), T=3, L=3, rho = 13.
Exhaustive search for a zonotope with sum L_i >= 13 (then zeta = rho, unfaithful), else zeta <= 12 (faithful m=13)."""
import time
from zono import zeta_search
t0 = time.time()
prof = ((0, 2), (0, 3), (1,), (2,), (2, 4), (3, 4), (4,))
z = zeta_search(5, 3, prof, 3, 13)
print("n=5", prof, "T=3 L=3 target 13:", z, "(%.0fs)" % (time.time() - t0))
