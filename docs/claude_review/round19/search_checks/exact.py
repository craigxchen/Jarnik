# Exact gain of a unit pencil f = (p/q) A / C, i.e. B = q*C - p*A (same fibre structure as C - lam A).
# This is the referee's formula (round18/referee_theory_checks/ref_gain.py and ref_gain_lambda.py, function gain),
# copied verbatim in substance: Q-factorisation by flint, nu_T of interior factors by exact symbolic substitution
# x_i -> x_{i0} + e*x_i (i in T \ {i0}) and lowest power of e (no specialisation).
import collections
from flint import fmpz_mpoly_ctx
from common import edges, Ts, nu_unit

ctx = fmpz_mpoly_ctx.get(('x', 8), 'lex')  # x0..x6 = labels, x7 = e
G = ctx.gens(); X = G[:7]; E = G[7]
lin = {e: X[e[0]] - X[e[1]] for e in edges}
lin_sig = {}
for e in edges:
    lin_sig[str(lin[e])] = e; lin_sig[str(-lin[e])] = e
subs = {}
for T in Ts:
    i0 = min(T)
    subs[T] = [X[i] if (i not in T or i == i0) else X[i0] + E * X[i] for i in range(7)] + [E]

def nuT(p, T):
    q = p.compose(*subs[T])
    return min(m[7] for m in q.monoms())

def poly(ms):
    p = ctx.from_dict({(0,) * 8: 1})
    for e in ms: p = p * lin[e]
    return p

def exact_gain(A, C, pp, qq, detail=False):
    """gain of the unit pencil with B = qq*C - pp*A (lam = pp/qq).  Returns None if B == 0."""
    pa = poly(A); pc = poly(C); B = qq * pc - pp * pa
    if B.is_zero(): return None
    c, facs = B.factor()
    pairB = collections.Counter(); inter = []
    for f, m in facs:
        if f.total_degree() == 1 and str(f) in lin_sig: pairB[lin_sig[str(f)]] += m
        elif f.total_degree() >= 1: inter.append((f, m))
    nuY = {(i, T): nuT(f, T) for i, (f, m) in enumerate(inter) for T in Ts}
    D = len(A); Lb = 31 * D; Fb = 0; Ecomp = []
    for T in Ts:
        a = nu_unit(A, T); cc = nu_unit(C, T); mu = min(a, cc); Lb -= mu
        b = sum(m for (u, v), m in pairB.items() if u in T and v in T) + \
            sum(m * nuY[(i, T)] for i, (f, m) in enumerate(inter))
        if max(a, b, cc) > mu:
            Fb += 1
            if a == cc: Ecomp.append((sorted(T), b - mu))   # E_T in the fibre over 1, multiplicity b - mu
    for e in edges:
        if e in A or e in C or e in pairB: Fb += 1
    ybeta = [31 * f.total_degree() - sum(nuY[(i, T)] for T in Ts) for i, (f, m) in enumerate(inter)]
    Fb += sum(ybeta)
    g = int(Lb - Fb)
    if not detail: return g
    return dict(gain=g, Lbeta=Lb, Fredbeta=Fb, B=str(B), content=str(c),
                pair_factors={str(e): m for e, m in pairB.items()},
                interior=[dict(Y=str(f), mult=m, deg=f.total_degree(), Ybeta=yb) for (f, m), yb in zip(inter, ybeta)],
                E_in_F1=Ecomp)
