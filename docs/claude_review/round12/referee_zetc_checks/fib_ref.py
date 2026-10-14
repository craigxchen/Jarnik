# Independent recomputation of the fibration bound (zetc.md Lemma 5.1) for Syl(16,L), L=2,3.
import itertools, flint
from fractions import Fraction
def newton_degrees(B, maxdeg=30):
    # e_i = least t such that eval at B[i] of monomials deg<=t not in span of evals at B[0..i-1]
    # computed incrementally: maintain for each t the echelon of previous points' eval vectors.
    mons = [e for t in range(maxdeg+1) for e in itertools.product(range(t+1), repeat=3) if sum(e) == t]
    def ev(b, t):
        return [b[0]**e[0]*b[1]**e[1]*b[2]**e[2] for e in mons if sum(e) <= t]
    es = []
    for i, b in enumerate(B):
        for t in range(maxdeg+1):
            rows = [ev(bb, t) for bb in B[:i]]
            M0 = flint.fmpq_mat(rows) if rows else None
            r0 = M0.rank() if rows else 0
            r1 = flint.fmpq_mat(rows + [ev(b, t)]).rank()
            if r1 > r0: es.append(t); break
        else: raise Exception
    return es
for L in (2, 3):
    B = [b for b in itertools.product(range(-4, 5), repeat=3)
         if abs(b[0]+b[1]) <= L and abs(b[0]+b[2]) <= L and abs(b[1]+b[2]) <= L]
    # fibre degree delta = sum (4 - |b_R|) (slices of 4 coordinates, degree <= 4-|j|)
    B.sort(key=lambda b: (sum(map(abs, b)), b))
    es = newton_degrees(B)
    vals = [e + 12 - sum(map(abs, b)) for e, b in zip(es, B)]
    print("L", L, "|B|", len(B), "max(e_i+delta_i) =", max(vals), "=> rho(Syl16,L) <=", max(vals) + 6, flush=True)
