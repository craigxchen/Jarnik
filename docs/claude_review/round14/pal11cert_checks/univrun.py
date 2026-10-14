"""Universal count search for pal11 (L=2, m=8): one exhaustive run per feasible (sigma,t) profile (profiles2.py).

Explicit types (all flats have 1-dimensional complement, so the lists contain EVERY generator of the type):
  E = (2,6) (4,6) (4,9) (6,2) (6,4) (7,5) (12,1).
Hard (parametric) types: (3,9) (4,8) (4,10) (6,6) (7,6) (8,4) (8,5) (8,6) (12,2).
Method per profile p:
  * all types explicit: enumerate all 8-multisets with exactly p's type composition (FREE = 0);
  * p = 7 pairs + (12,2): enumerate 7 pairs + one representative of the (12,2) planes (reps_12_2.pkl; dominance:
    the representative touches a subset of what the original touches, so the replaced multiset is count-feasible);
  * otherwise p has exactly one hard slot: enumerate the other 7 (explicit types) and require that some arbitrary
    nonzero 8th generator fits (FREE = 1; saturated directions must have rank <= 10).
Every count-feasible multiset with profile p is found (or a count-feasible dominated replacement is), so if every
run has 0 leaves, NO product of 8 real linear forms satisfies the width counts, hence there is no product
certificate at all.  usage: python3 univrun.py [profile indices]"""
import sys, pickle, subprocess, time
from math import gcd
from collections import Counter
import numpy as np
from fplib import pal_patterns
exec(open('profiles2.py').read().split("print('feasible profiles:'")[0])
dat = np.load('dirs_pal11.npz'); D = dat['D']; tau = dat['tau']
M_, pats = pal_patterns(11)
core = [P for P in pats if len(P) == 6]
Cv = np.zeros((11, 11), dtype=np.int64)
for c, P in enumerate(core):
    for x in P: Cv[c, x - 1] = 1
EXPL = {(2, 6), (4, 6), (4, 9), (6, 2), (6, 4), (7, 5), (12, 1)}
def norm(z):
    g = 0
    for x in z: g = gcd(g, abs(int(x)))
    z = [int(x) // g for x in z]
    f = next(x for x in z if x)
    return tuple(z) if f > 0 else tuple(-x for x in z)
lists = {ty: set() for ty in EXPL}
for (s, t, R, Z, dim, basis) in pickle.load(open('typeflats_2_3_4_5_6_7_8.pkl', 'rb')) + pickle.load(open('excflats.pkl', 'rb')):
    if (s, t) in EXPL:
        assert dim == 1
        z = [0] * 12
        for x, v in zip(R, basis[0]): z[x] = v
        lists[(s, t)].add(norm(z))
lists = {ty: sorted(v) for ty, v in lists.items()}
lists['reps12_2'] = sorted(set(pickle.load(open('reps_12_2.pkl', 'rb'))))
print({k: len(v) for k, v in lists.items()}, flush=True)
def tc(z): return int(((Cv @ np.array(z[1:])) != 0).sum())
def act(a, z):
    w = [0] * 12; w[0] = z[0]
    for i in range(11): w[1 + (a * i) % 11] = z[1 + i]
    return norm(w)
def write_pool(classes, fn):
    """classes: list of (classname, list of z).  Pool sorted by (t, class, z); orbit reps within class."""
    items = []
    for ci, (name, L) in enumerate(classes):
        for z in L: items.append((tc(z), ci, z))
    items.sort()
    index = {(ci, z): i for i, (t, ci, z) in enumerate(items)}
    rep = [1 if i == min(index[(ci, act(a, z))] for a in (1, 3, 4, 5, 9)) else 0 for i, (t, ci, z) in enumerate(items)]
    P = np.array([z for (t, ci, z) in items], dtype=np.int64)[:, 1:]
    NZ = (P @ D.T) != 0
    with open(fn, 'w') as f:
        f.write(f'{len(items)} {len(D)}\n')
        f.write(' '.join(str(int(t)) for t in tau) + '\n')
        for i, (t, ci, z) in enumerate(items):
            idx = np.nonzero(NZ[i])[0]
            f.write(f'{rep[i]} ' + ' '.join(map(str, z)) + f' {len(idx)} ' + ' '.join(map(str, idx)) + '\n')
    with open(fn + '.cls', 'w') as f:
        f.write(' '.join(str(ci) for (t, ci, z) in items) + '\n')
    return len(items), sum(rep)
sel = [int(a) for a in sys.argv[1:]] if len(sys.argv) > 1 else range(len(full))
summary = open('univrun_summary.txt', 'a')
for pi in sel:
    prof = full[pi]
    cnt = Counter(prof)
    hard = [ty for ty in cnt if ty not in EXPL]
    reserve = None
    if not hard:
        classes = [(str(ty), lists[ty]) for ty in sorted(cnt)]; comp = [cnt[ty] for ty in sorted(cnt)]; m = 8; free = 0
        mode = 'full'
    elif hard == [(12, 2)] and cnt[(12, 2)] == 1:
        rest = [ty for ty in sorted(cnt) if ty != (12, 2)]
        classes = [(str(ty), lists[ty]) for ty in rest] + [('reps12_2', lists['reps12_2'])]
        comp = [cnt[ty] for ty in rest] + [1]; m = 8; free = 0; mode = 'full+reps(12,2)'
    else:
        assert len(hard) == 1 and cnt[hard[0]] == 1, prof
        rest = [ty for ty in sorted(cnt) if ty != hard[0]]
        classes = [(str(ty), lists[ty]) for ty in rest]; comp = [cnt[ty] for ty in rest]; m = 7; free = 1
        mode = 'free slot ' + str(hard[0])
        sf, tf = hard[0]
        reserve = [tf, 0, tf, sf, tf, 0, tf + sf, tf + sf, 2 * tf + sf, 3 * tf + sf]
    fn = f'upool_{pi}.txt'
    npool, nrep = write_pool(classes, fn)
    t0 = time.time()
    procs = []
    for part in range(4):
        procs.append(subprocess.Popen(['./csearch9', fn, 'weights.txt', str(m), f'uleaves_{pi}_{part}.txt', '4', str(part),
                                       'dirs_pal11.txt', fn + '.cls', ','.join(map(str, comp)), str(free)] + ([','.join(map(str, reserve))] if reserve else []),
                                      stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True))
    outs = [p.communicate()[0].strip() for p in procs]
    nleaves = sum(sum(1 for _ in open(f'uleaves_{pi}_{part}.txt')) for part in range(4))
    line = f'profile {pi} {prof} | mode {mode} reserve {reserve} | pool {npool} (reps {nrep}) | written leaves {nleaves} | {time.time()-t0:.0f}s'
    print(line, flush=True); summary.write(line + '\n')
    for o in outs: summary.write('    ' + o + '\n')
    summary.flush()
