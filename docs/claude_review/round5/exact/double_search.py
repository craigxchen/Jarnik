"""Search for double-Hadamard supports T = core_ref(H) cup core_ref'(P H) (row permutation P,
reference columns ref, ref') satisfying the integral margin condition
      F_1(c) + F_2(c) >= 2M - 2 + ||c||_1        (sigma_c >= ||c||_1/2 at unit weights)
on all ternary c with support <= smax (fast float screen), then exact analysis (inner.analyse) of
the survivors, including rational classes and windings.
usage: python3 double_search.py name ntrials smax seed"""
import sys, random, itertools, time
import numpy as np
from fractions import Fraction
from fixed import sylvester, paley
from core import pattern_of, Group
from inner import analyse, find_hom


def core_ref(H, ref):
    M = len(H)
    out = []
    for a in range(M):
        if a == ref:
            continue
        out.append(pattern_of(M, [x for x in range(M) if H[x][a] * H[x][ref] == -1]))
    return out


def tern(M, smax):
    rows = []
    for h in range(2, smax + 1, 2):
        for supp in itertools.combinations(range(M), h):
            for pos in itertools.combinations(supp, h // 2):
                if supp[0] not in pos:
                    continue
                c = np.zeros(M, dtype=np.int64)
                for x in supp:
                    c[x] = 1 if x in pos else -1
                rows.append(c)
    return np.array(rows)


name = sys.argv[1]; ntr = int(sys.argv[2]); smax = int(sys.argv[3]); seed = int(sys.argv[4])
H = {'S8': lambda: sylvester(3), 'P12': lambda: paley(11), 'S16': lambda: sylvester(4), 'P20': lambda: paley(19)}[name]()
M = len(H)
C = tern(M, smax)
n1 = np.abs(C).sum(axis=1)
rng = random.Random(seed)
c1 = core_ref(H, 0)
best = None
t0 = time.time()
for t in range(ntr):
    perm = list(range(M)); rng.shuffle(perm)
    ref = rng.randrange(M)
    H2 = [H[perm[x]] for x in range(M)]
    c2 = core_ref(H2, ref)
    if set(c1) & set(c2):
        continue
    T = c1 + c2
    P = np.array([[(s >> x) & 1 for x in range(M)] for s in T], dtype=np.int64)
    S = np.abs(C @ P.T).sum(axis=1)            # sum_S |c(S)|
    slack2 = 2 * S - 2 * (M - 1) - n1          # 2 sigma_c - ||c||_1  (unit weights)
    m = int(slack2.min())
    nbad = int((slack2 < 0).sum())
    if best is None or (nbad, -m) < best[0]:
        best = ((nbad, -m), perm, ref)
    if nbad == 0:
        G = Group(M, T)
        print(f'[{time.time()-t0:.0f}s] trial {t}: screen pass (min 2sigma-n = {m}), group {G.invs if G.full_rank else "rank"}', flush=True)
        if G.full_rank and G.order <= 5000:
            A = analyse(M, T, [Fraction(1)] * len(T), 0, 1)
            if A.get('inconclusive'):
                print('   inconclusive'); continue
            eta = find_hom(G, list(A['bad'].keys())) if not A['int_viol'] else None
            print('   exact: int_viol', len(A['int_viol']), 'bad rational classes', len(A['bad']), 'hom', eta,
                  'min slack', A['min_slack'][0] if A['min_slack'] else None, flush=True)
            if not A['int_viol'] and eta is not None:
                print('FOUND Plotkin-tight admissible double core: perm', perm, 'ref', ref, flush=True)
                break
print('best screen', best[0])
