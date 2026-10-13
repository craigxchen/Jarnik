"""Exact ingredients for the failure of ZE-real in the codimension-one model sdc7 (T=1, L=2, m=8, dim V = 1).
(1) the single basis vector c is even (c_{-k} = c_k), so F(u) = sum_k c_k e^{i<k,u>} = sum_k c_k cos<k,u> is real;
(2) Phi(u) = sum_k c_k <k,u>^8 takes both signs: exact integer evaluation at integer points u+, u-, both
    nondegenerate-direction irrelevant (any points suffice);
(3) P_c = (monomial) * G(w - 1) with G an irreducible form of degree 8 over Q (from cone_reduce: one non-binomial
    factor, form of degree 8 = m), so Phi = 8! G has no linear factor over Q (checked again here by factoring Phi).
"""
import json, itertools, math
import numpy as np
import flint
D = json.load(open("sdc7_T1_L2_m8_basis.json")); A = [tuple(k) for k in D['A']]; c = D['basis'][0]; d = 7
idx = {k: i for i, k in enumerate(A)}
even = all(c[idx[tuple(-x for x in k)]] == c[i] for i, k in enumerate(A))
print("c even:", even, " support:", sum(1 for x in c if x), " sup:", max(abs(x) for x in c))
def Phi(u):
    return sum(ck * sum(a * b for a, b in zip(k, u)) ** 8 for k, ck in zip(A, c) if ck)
npos = nneg = nzero = 0
zeros = []
for u in itertools.product(range(-2, 3), repeat=d):
    if not any(u):
        continue
    v = Phi(u)
    if v > 0: npos += 1
    elif v < 0: nneg += 1
    else:
        nzero += 1
        if len(zeros) < 2000: zeros.append(u)
print("exact signs of Phi on the 5^7-1 integer points of [-2,2]^7: positive %d, negative %d, zero %d" % (npos, nneg, nzero))
supp = {}
for u in zeros:
    key = sum(1 for x in u if x)
    supp[key] = supp.get(key, 0) + 1
print("support sizes of the integer zeros:", supp)
print("some zeros with support >= 4:", [u for u in zeros if sum(1 for x in u if x) >= 4][:10])
# Phi as an exact polynomial, factored over Q
ctx = flint.fmpz_mpoly_ctx.get(('u', d), 'lex')
from core import monomials
mons = monomials(d, 8)
coef = {}
for a in mons:
    mult = math.factorial(8)
    for ai in a: mult //= math.factorial(ai)
    s = sum(ck * np.prod([k[t] ** a[t] for t in range(d)]) for k, ck in zip(A, c) if ck)
    if s: coef[a] = int(mult * s)
P = ctx.from_dict(coef)
cont, facs = P.factor()
print("Phi factorisation over Q: content", cont, "; factors (total degree, #terms, multiplicity):",
      [(f.total_degree(), len(f.to_dict()), e) for f, e in facs])
