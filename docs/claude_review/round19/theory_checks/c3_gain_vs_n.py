"""c3: consistency test of Theorem A (gain(f) <= 31 - n(f)) on all unit pencils s*A/C, s=+-1, D <= DMAX,
plus special rational lambda (passed on the command line as A-C-lambda triples is not needed: we also
scan lambda in a small set of special values for every (A,C)).  Adapted from round18 ref_gain.py (exact:
Q-factorisation by python-flint, nu_T by exact substitution).  Output: jsonl per A-orbit, resumable.
n(f) = number of boundary divisors in the three fibres.  Reports max of gain + n (theorem: <= 31).
"""
import itertools, sys, time, collections, json, os
from fractions import Fraction
from flint import fmpz_mpoly_ctx, fmpq_mpoly_ctx
D = int(sys.argv[1]); OUT = sys.argv[2]; ONLY = int(sys.argv[3]) if len(sys.argv) > 3 else None
C0 = int(sys.argv[4]) if len(sys.argv) > 4 else 0; C1 = int(sys.argv[5]) if len(sys.argv) > 5 else 10**9
LAMS = [Fraction(1), Fraction(-1), Fraction(2), Fraction(-2), Fraction(1,2), Fraction(-1,2), Fraction(4), Fraction(-4),
        Fraction(1,4), Fraction(-1,4), Fraction(3), Fraction(-3), Fraction(1,3), Fraction(-1,3)]
ctx = fmpz_mpoly_ctx.get(('x', 8), 'lex')
G = ctx.gens(); X = G[:7]; E = G[7]
edges = [(i, j) for i in range(7) for j in range(i+1, 7)]
lin = {e: X[e[0]] - X[e[1]] for e in edges}
Ts = [frozenset(c) for r in range(3, 7) for c in itertools.combinations(range(7), r)]
subs = {}
for T in Ts:
    i0 = min(T); subs[T] = [X[i] if i not in T or i == i0 else X[i0] + E*X[i] for i in range(7)] + [E]
def nuT(p, T):
    q = p.compose(*subs[T]); return min(m[7] for m in q.monoms())
def canon(ms):
    verts = sorted(set(v for e in ms for v in e)); k = len(verts); best = None
    for perm in itertools.permutations(range(k)):
        mp = {verts[i]: perm[i] for i in range(k)}
        t = tuple(sorted(tuple(sorted((mp[a], mp[b]))) for (a, b) in ms))
        if best is None or t < best: best = t
    return best
def poly(ms):
    p = ctx.from_dict({(0,)*8: 1})
    for e in ms: p = p*lin[e]
    return p
def nu_unit(ms, T): return sum(1 for (a, b) in ms if a in T and b in T)
lin_sig = {}
for e in edges: lin_sig[str(lin[e])] = e; lin_sig[str(-lin[e])] = e
done = set()
if os.path.exists(OUT):
    for line in open(OUT): done.add(tuple(map(tuple, json.loads(line)['A'])))
reps = sorted(set(canon(ms) for ms in itertools.combinations_with_replacement(edges, D)))
for ia, A in enumerate(reps):
    if ONLY is not None and ia != ONLY: continue
    if tuple(A) in done and ONLY is None: continue
    Aset = set(A); pa = poly(A); others = [e for e in edges if e not in Aset]
    cnt = 0; worst = None; hist = collections.Counter(); gmax = None
    for ic, C in enumerate(itertools.combinations_with_replacement(others, D)):
        if ic < C0 or ic >= C1: continue
        pc = poly(C)
        for lam in LAMS:
            # B = q*C - p*A  for lambda = p/q  (the fibre over 1 of lambda*A/C)
            B = lam.denominator*pc - lam.numerator*pa
            c, facs = B.factor()
            pairB = collections.Counter(); inter = []
            for f, m in facs:
                if f.total_degree() == 1 and str(f) in lin_sig: pairB[lin_sig[str(f)]] += m
                else: inter.append((f, m))
            if lam not in (1, -1) and not pairB and all(m == 1 for f, m in inter) and len(inter) == 1:
                # possibly generic lambda; still evaluate (cheap enough at D<=3)
                pass
            Lb = 31*D; Fb_bd = 0; Fb_int = 0; cache = {}
            for T in Ts:
                a = nu_unit(A, T); cc = nu_unit(C, T); mu = min(a, cc); Lb -= mu
                b = sum(m for (u, v), m in pairB.items() if u in T and v in T)
                for f, m in inter:
                    if (str(f), T) not in cache: cache[(str(f), T)] = nuT(f, T)
                    b += m*cache[(str(f), T)]
                if max(a, b, cc) > mu: Fb_bd += 1
            for e in edges:
                if e in Aset or e in C or e in pairB: Fb_bd += 1
            for f, m in inter:
                Fb_int += 31*f.total_degree() - sum(cache[(str(f), T)] for T in Ts)
            g = int(Lb - Fb_bd - Fb_int); Fb_bd = int(Fb_bd); cnt += 1
            s = g + Fb_bd; hist[s] += 1
            if worst is None or s > worst[0]: worst = (s, g, Fb_bd, list(C), str(lam))
            if gmax is None or g > gmax: gmax = g
    rec = {'A': [list(e) for e in A], 'chunk': [C0, C1], 'pencils': cnt, 'max_gain_plus_n': worst[0], 'worst': worst, 'max_gain': gmax,
           'hist_gain_plus_n_top': sorted(hist.items(), reverse=True)[:5]}
    with open(OUT, 'a') as fo: fo.write(json.dumps(rec) + '\n')
    print(json.dumps(rec), flush=True)
