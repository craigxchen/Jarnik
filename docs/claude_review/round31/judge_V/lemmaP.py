# Lemma P sanity check over Z: A subset of {+-1..+-6}, A cap -A empty, size <= 6: when is prod(X+a)-prod(X-a) constant?
import itertools
from sympy import symbols, expand, Poly
X = symbols('X'); found = {}
els = [s * k for k in range(1, 7) for s in (1, -1)]
for n in range(1, 7):
    for A in itertools.combinations(els, n):
        if any(-a in A for a in A): continue
        p, q = 1, 1
        for a in A: p *= (X + a); q *= (X - a)
        if Poly(expand(p - q), X).degree() <= 0: found.setdefault(n, []).append(A)
print({n: len(v) for n, v in found.items()}, "examples size 3:", found.get(3, [])[:3])
