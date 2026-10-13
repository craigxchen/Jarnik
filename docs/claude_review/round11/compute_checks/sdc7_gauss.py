"""Factor the sdc7 tangent-cone form Phi (and the base-locus polynomial in x = w - 1) over Q(i) (sympy, exact).
If Phi = -q * conj(q) with q a quartic over Q(i), the real zero set of Phi is {Re q = Im q = 0}."""
import json, math, time
import sympy as sp
D = json.load(open("sdc7_T1_L2_m8_basis.json")); A = D['A']; c = D['basis'][0]; d = 7
U = sp.symbols('u0:7')
Phi = sp.expand(sum(ck * sum(k[t] * U[t] for t in range(d)) ** 8 for k, ck in zip(A, c) if ck))
t0 = time.time()
print("Phi has", len(sp.Poly(Phi, *U).terms()), "terms")
f = sp.factor_list(Phi, gaussian=True)
print("factorisation over Q(i) (%.1fs):" % (time.time() - t0))
print("  content:", f[0])
for g, e in f[1]:
    print("  factor of degree", sp.Poly(g, *U).total_degree(), "with", len(sp.Poly(g, *U).terms()), "terms, multiplicity", e)
    if sp.Poly(g, *U).total_degree() <= 4:
        print("   ", g)
