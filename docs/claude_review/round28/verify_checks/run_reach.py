"""reach of (SHIFT_e) for a slab box, mod p (own code).
usage: run_reach.py SPEC PRIME SHIFTS [EMIN] [EMAX] [MAXCLASSES]
SHIFTS = pair | cube | paircube (cube = all {-1,0,1}^d classes incl. pair ones)
Output: one JSON line appended to verify_checks/reach_<prime>.jsonl
"""
import sys, json, time
import numpy as np
import flint
sys.path.insert(0, '/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round28/verify_checks')
from vlib import *

spec = sys.argv[1]; p = int(sys.argv[2]); which = sys.argv[3]
emin = int(sys.argv[4]) if len(sys.argv) > 4 else 0
emax_arg = int(sys.argv[5]) if len(sys.argv) > 5 else None
maxcl = int(sys.argv[6]) if len(sys.argv) > 6 else None
t0 = time.time()
n, pats, L, A = spec_box(spec)
d = n - 1
log(spec, 'n', n, 'cuts', len(pats), 'L', L, '|A|', len(A))
# H_A for all degrees (emax = 2d is safe; use rho bound)
by = std_monos(d, 2 * d)
# degree cap for H_A: smallest e with N_{<=e} >= |A|, then grow until sum(H) == |A| (so H is complete)
hg = 0; tot = len(by[0])
while tot < len(A) and hg < 2 * d: hg += 1; tot += len(by[hg])
hg = min(2 * d, hg + 1)
while True:
    blkA = least_blocks(A, by, p, hg)
    if sum(b.shape[0] for b in blkA) == len(A) or hg == 2 * d: break
    hg = min(2 * d, hg + 2)
H = [b.shape[0] for b in blkA]
while H and H[-1] == 0: H.pop()
rho = len(H) - 1
assert sum(H) == len(A)
log('H_A', H, 'rho', rho)
emax = rho - 1 if emax_arg is None else min(emax_arg, rho)
index = {tuple(k.tolist()): i for i, k in enumerate(A)}
classes = lpair_classes(d) if which == 'pair' else cube_classes(d)
sets = []
for l in classes:
    Al = shifted(A, np.array(l), index)
    if len(Al) > emin: sets.append((len(Al), l, Al))
sets.sort(key=lambda t: -t[0])
if maxcl: sets = sets[:maxcl]
log('classes', len(classes), 'used', len(sets))
cache = {}
def blocks_of(i):
    if i not in cache:
        Al = sets[i][2]
        b = least_blocks(Al, by, p, emax)
        cache[i] = [x.astype(np.uint32) if e >= emin else None for e, x in enumerate(b)]
    return cache[i]

res = {'spec': spec, 'p': p, 'shifts': which, 'n': n, 'cuts': len(pats), 'L': L, 'size': len(A), 'H': H,
       'rho': rho, 'emin': emin, 'emax': emax, 'nclasses_used': len(sets), 'dims': {}}
reach = None
for e in range(emin, emax + 1):
    if H[e] == 0:
        res['dims'][e] = [0, 0]; continue
    basis = None; r = 0; used = 0
    for i in range(len(sets)):
        if sets[i][0] <= e: break
        b = blocks_of(i)[e]
        if b.shape[0] == 0: continue
        used += 1
        rows = (basis.tolist() if basis is not None else []) + b.astype(np.int64).tolist()
        M = flint.nmod_mat(rows, p)
        R, r = M.rref()
        basis = np.array([[int(R[a, c]) for c in range(R.ncols())] for a in range(r)], dtype=np.int64)
        if r == H[e]: break
    # containment check: span(basis) inside M_e(A)
    if basis is not None:
        rows = blkA[e].tolist() + basis.tolist()
        rc = flint.nmod_mat(rows, p).rank()
        assert rc == H[e], ('containment failed', e, rc, H[e])
    res['dims'][e] = [r, H[e]]
    log('e', e, 'sum', r, 'H', H[e], 'shifts used', used)
    if r < H[e]:
        reach = e; break
if reach is None:
    reach = emax + 1 if emax == rho - 1 else ('>=%d' % (emax + 1))
    if emax == rho - 1:
        reach = rho  # SHIFT_rho must fail (V(A,rho+1)=0); record reach = rho
res['reach'] = reach
res['time'] = round(time.time() - t0, 1)
log('RESULT', spec, which, 'reach', reach, 'time', res['time'])
with open('/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round28/verify_checks/reach_%d.jsonl' % p, 'a') as f:
    f.write(json.dumps(res) + '\n')
