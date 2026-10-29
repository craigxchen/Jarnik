"""Theory-lane checks on M8, M8', sdc7: base coordinate sigma = sum k' = -k_base; range, top level, rigidity (R),
number of legal pair moves; EXACT (fmpz) membership sigma^e in J(S)_e for S = A and S = A_l (both signs, all L_pair)."""
import sys
import numpy as np
import flint
sys.path.insert(0, '/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round28/verify_checks')
from vlib import *

def in_J_exact(S, phi, e, by):
    """phi^e in J(S)_e  <=>  (phi(k)^e)_k in column span of E_{<e}(S) over Q (exact fmpz ranks)."""
    monos = [a for t in range(e) for a in by[t]]
    E = evalmat(S, monos).astype(np.int64)
    v = (S.astype(np.int64) @ np.array(phi, dtype=np.int64)) ** e
    r1 = flint.fmpz_mat(E.tolist()).rank()
    r2 = flint.fmpz_mat(np.concatenate([E, v[:, None]], axis=1).tolist()).rank()
    return r1 == r2

def top_info(A):
    d = A.shape[1]
    sig = A.astype(int).sum(axis=1)
    h = sig.max()
    idx = {tuple(k.tolist()): i for i, k in enumerate(A)}
    top = A[sig == h]
    moves = 0; rigid_pts = 0
    for k in top:
        cnt = 0
        for x in range(d):
            for y in range(d):
                if x != y:
                    t = k.astype(int).copy(); t[x] += 1; t[y] -= 1
                    if tuple(t.tolist()) in idx: cnt += 1
        moves += cnt
        rigid_pts += (cnt == 0)
    return int(sig.min()), int(h), len(top), moves, rigid_pts, top

for spec in ['syl:3:1', 'syl:3:2', 'syl:3:3', 'sdc:7:1', 'sdc:7:2', 'sdc:7:3']:
    n, pats, L, A = spec_box(spec)
    mn, h, ntop, moves, rp, top = top_info(A)
    print(spec, '|A|', len(A), 'sigma range', [mn, h], 'top points', ntop, 'legal pair moves from top', moves,
          '(R) holds' if moves == 0 else '(R) fails', flush=True)
    if spec == 'syl:3:2':
        # top level = indicators of Fano lines?  lines = complements of the 7 patterns (patterns omit base 7)
        pat_side = [tuple(sorted(set(range(7)) - set(S))) if 7 in S else tuple(S) for S in pats]
        lines = sorted(tuple(sorted(set(range(7)) - set(P))) for P in pat_side)
        tops = sorted(tuple(int(x) for x in np.nonzero(k)[0]) for k in top)
        print('  top level points are 0/1 vectors:', all(set(k.tolist()) <= {0, 1} for k in top),
              ' supports == Fano lines (complements of patterns):', tops == lines, tops)
        d = 7; by = std_monos(d, 8)
        sig = [1] * d
        print('  EXACT sigma^7 in J(A)_7:', in_J_exact(A, sig, 7, by), flush=True)
        print('  EXACT sigma^6 in J(A)_6:', in_J_exact(A, sig, 6, by), flush=True)
        idx = {tuple(k.tolist()): i for i, k in enumerate(A)}
        allin = True; ranges = set()
        for l in lpair_classes(d):
            for s in (1, -1):
                Al = shifted(A, np.array(l) * s, idx)
                sg = Al.astype(int).sum(axis=1); ranges.add((int(sg.min()), int(sg.max())))
                ok = in_J_exact(Al, sig, 6, by)
                allin &= ok
        print('  EXACT sigma^6 in J(A_l)_6 for all 56 signed L_pair shifts:', allin, ' sigma ranges on A_l:', sorted(ranges), flush=True)
