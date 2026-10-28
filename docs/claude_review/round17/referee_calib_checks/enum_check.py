#!/usr/bin/env python3
"""Referee: brute-force enumeration of A_L(q) (independent of calib's enumA.c).
Patterns: P_a = {a} u {x: chi(a-x)=1}, P'_a = {a} u {x: chi(a-x)=-1}.
Writes A to data/A_q{q}_L{L}.npy (full list, int8)."""
import sys, itertools, numpy as np, time, os
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data')
os.makedirs(OUT, exist_ok=True)

def chi(a, q):
    a %= q
    if a == 0: return 0
    return 1 if pow(a, (q-1)//2, q) == 1 else -1

def patterns(q):
    P = []
    for a in range(q):
        P.append([a] + [x for x in range(q) if chi(a-x, q) == 1])
        P.append([a] + [x for x in range(q) if chi(a-x, q) == -1])
    return P

def enum(q, L):
    P = patterns(q)
    Pm = np.zeros((2*q, q), dtype=np.int64)
    for i, S in enumerate(P):
        Pm[i, S] = 1
    # vectorised brute force in chunks over first coordinates
    res = []
    vals = np.array([-1, 0, 1], dtype=np.int8)
    # build all of {-1,0,1}^q in blocks: split into head (h coords) and tail
    t = min(q, 11)
    tail = np.array(list(itertools.product([-1, 0, 1], repeat=t)), dtype=np.int64)
    tailS = tail @ Pm[:, q-t:].T
    h = q - t
    for head in itertools.product([-1, 0, 1], repeat=h):
        hv = np.array(head, dtype=np.int64)
        s = tailS + (Pm[:, :h] @ hv if h else 0)
        ok = np.all(np.abs(s) <= L, axis=1)
        if ok.any():
            sel = tail[ok]
            res.append(np.hstack([np.tile(hv, (sel.shape[0], 1)), sel]) if h else sel)
    A = np.vstack(res).astype(np.int8)
    return A

if __name__ == '__main__':
    q, L = int(sys.argv[1]), int(sys.argv[2])
    t0 = time.time()
    A = enum(q, L)
    np.save(os.path.join(OUT, f'A_q{q}_L{L}.npy'), A)
    supp = (A != 0).sum(1)
    print(f'q={q} L={L}: |A|={len(A)}  max supp={supp.max()}  time={time.time()-t0:.1f}s', flush=True)
    if L == 1:
        # check A_1 = {0} u {+-e_x} u {e_x - e_y}
        S = set(map(tuple, A.tolist()))
        T = {tuple([0]*q)}
        for x in range(q):
            e = [0]*q; e[x] = 1; T.add(tuple(e)); e[x] = -1; T.add(tuple(e))
            for y in range(q):
                if y != x:
                    e = [0]*q; e[x] = 1; e[y] = -1; T.add(tuple(e))
        print('  A_1 == {0,+-e_x,e_x-e_y}:', S == T, ' 1+q(q+1) =', 1+q*(q+1))
