"""Multiplicities in the reduced points z_j (and in N) at the shared-root places alpha in F_729, best5 at L=32."""
import numpy as np, itertools, collections
import descent_check as dc
from fullgcd import hasse_ord
from fams import fams
K = dc.K; F27 = dc.F27; I = dc.I
forms, signs, g = fams['best5']; L = 32
vs = list(itertools.product(F27, repeat=3)); X = np.array(vs)
IX = K.mul(np.full(X.shape, I), X)
E = np.array([[1 + e*dc.chi27[dc.lam_val(l, x)] for x in vs] for l, e in zip(forms, signs)])
mult = {}
for alpha in range(729):
    aL = K.pow(alpha, L) if alpha else 0
    val = IX[:, 0]
    if alpha:
        for b in (1, 2): val = K.add[val, K.mul(IX[:, b], np.full(len(vs), K.pow(alpha, b)))]
    fib = np.nonzero(val == aL)[0]
    if len(fib) == 0: continue
    ords = np.array([hasse_ord(alpha, list(IX[idx]), L) for idx in fib])
    s = [int((E[j, fib]*ords).sum()) for j in range(5)]
    mult[alpha] = [x - min(s) for x in s]
frob3 = lambda a: K.pow(a, 3) if a else 0
hist = collections.Counter(max(m) for m in mult.values())
print('max multiplicity of (t-alpha) in some z_j, over alpha in F_729 with nonempty fibre: histogram', sorted(hist.items()))
Nm = collections.Counter()
for a in mult:
    b = [x for x in mult if frob3(x) == a]
    for j in range(5):
        Nm[mult[a][j] + (mult[b[0]][j] if b else 0)] += 0
    Nm[mult[a][0] + (mult[b[0]][0] if b else 0)] += 1
print('multiplicity of (t-alpha) in N (via z_1), histogram', sorted(Nm.items()))
