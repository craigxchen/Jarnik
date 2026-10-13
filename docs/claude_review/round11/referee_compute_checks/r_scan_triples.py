# Referee scan of admissible slab-box models at n = 3 and n = 4 (all sign profiles), conditions as in compute.md §3:
# (i) slabs cut, (ii') max range of k(S_p) < m, [(ii) 2*range < m], (iii) m <= rho, (v) full rank.
# For each distinct A: exact rho and a zonotope witness with sum L_i = rho (certifies zeta = rho).
import sys, itertools, time
from rzono import *
n = int(sys.argv[1]); Tmax = int(sys.argv[2]); Lmax = int(sys.argv[3])
base = n - 1
pts = list(range(n))
# non-constant patterns up to complement: represent by the side not containing ... (any rep works: |k(S)| = |k(S^c)|)
pats = []
seen = set()
for r_ in range(1, n):
    for S in itertools.combinations(pts, r_):
        S = frozenset(S); Sc = frozenset(pts) - S
        key = min(tuple(sorted(S)), tuple(sorted(Sc)))
        if key in seen: continue
        seen.add(key); pats.append(S)
print("patterns:", [sorted(S) for S in pats])
def full_rank(profile, others):
    import flint
    rows = [[(1 if x in S else 0) - (1 if base in S else 0) for x in others] for S in profile]
    return flint.fmpz_mat(rows).rank() == len(others)
inst = inst2 = 0
distinct = {}
t0 = time.time()
for r_ in range(1, len(pats) + 1):
    for prof in itertools.combinations(pats, r_):
        for T in range(1, Tmax + 1):
            for L in range(1, Lmax + 1):
                others, cons = functionals(n, list(prof), T, L, base)
                if not full_rank(prof, others): continue
                d = len(others)
                A = points(d, cons, T)
                if len(A) == (2 * T + 1) ** d: continue          # (i)
                key = tuple(A)
                if key not in distinct:
                    distinct[key] = {'rho': rho(A, T), 'cons': cons, 'T': T, 'ex': (prof, T, L)}
                rh = distinct[key]['rho']
                rngs = []
                for S in prof:
                    lam = [(1 if x in S else 0) - (1 if base in S else 0) for x in others]
                    vals = [sum(a * b for a, b in zip(lam, k)) for k in A]
                    rngs.append(max(vals) - min(vals))
                ms = [m for m in range(max(rngs) + 1, rh + 1)]
                ms2 = [m for m in range(2 * max(rngs) + 1, rh + 1)]
                if ms:
                    inst += 1; distinct[key]["adm"] = True
                if ms2:
                    inst2 += 1; distinct[key]["adm2"] = True
adm = {k: v for k, v in distinct.items() if v.get('adm')}
print(f"n={n} T<={Tmax} L<={Lmax}: (profile,T,L,m) instances with (i),(ii'),(iii),(v): {inst}; with (ii): {inst2}; "
      f"distinct admissible A: {len(adm)} ({sum(1 for v in adm.values() if v.get('adm2'))} also with (ii)) ({time.time()-t0:.0f}s)", flush=True)
nfound = 0; notfound = []
for key, v in []:
    d = len(key[0])
    w = find_zonotope(v['cons'], d, v['rho'], E=2)
    if w is None:
        w = find_zonotope(v['cons'], d, v['rho'], E=3, node_limit=5_000_000)
    if w is not None:
        nfound += 1
    else:
        notfound.append((v['ex'], v['rho'], len(key)))
print(f"zonotope witness with sum L_i = rho found for {nfound}/{len(adm)} distinct admissible A ({time.time()-t0:.0f}s)")
for x in notfound[:20]:
    print("  NOT FOUND:", [sorted(S) for S in x[0][0]], "T,L=", x[0][1:], "rho", x[1], "|A|", x[2])
