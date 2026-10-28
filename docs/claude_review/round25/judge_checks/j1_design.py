# Judge check J1 (exact): extended double-error-correcting BCH cut designs.
# Columns c_x = (1, x, x^3) in F_2^{1+2t}, x in GF(2^t). Claim (referee Prop R1):
# any 5 columns independent; every 5-subtuple exactly balanced; Z carries 2^{m-1} cuts (incl. trivial).
import itertools, sys
from math import comb

def gf_mul(a, b, t, poly):
    r = 0
    while b:
        if b & 1: r ^= a
        b >>= 1; a <<= 1
        if a >> t: a ^= poly
    return r

POLYS = {3: 0b1011, 4: 0b10011, 5: 0b100101}

def columns(t):
    poly = POLYS[t]
    cols = []
    for x in range(2 ** t):
        x3 = gf_mul(gf_mul(x, x, t, poly), x, t, poly)
        cols.append(1 | (x << 1) | (x3 << (1 + t)))
    return cols, 1 + 2 * t

def rank(vecs):
    basis = []
    for v in vecs:
        for b in basis:
            v = min(v, v ^ b)
        if v: basis.append(v)
    return len(basis)

def parity(x): return bin(x).count('1') & 1

for t in (3, 4, 5):
    cols, m = columns(t)
    M = len(cols)
    bad = 0; unbalanced = 0; n = 0
    us = [u for u in range(1, 2 ** m)]
    for S in itertools.combinations(range(M), 5):
        n += 1
        if rank([cols[x] for x in S]) < 5:
            bad += 1
    # balance check on a full sweep for t<=4, and a deterministic subsample for t=5
    sample = list(itertools.combinations(range(M), 5))
    if t == 5: sample = sample[::997]
    for S in sample:
        counts = {}
        for u in us:
            pat = tuple(parity(u & cols[x]) for x in S)
            counts[pat] = counts.get(pat, 0) + 1
        # each nonzero-pattern count should be 2^{m-5}; zero pattern 2^{m-5}-1 (u != 0)
        ok = all(counts.get(p, 0) == 2 ** (m - 5) for p in itertools.product((0, 1), repeat=5) if any(p)) \
             and counts.get((0,) * 5, 0) == 2 ** (m - 5) - 1
        if not ok: unbalanced += 1
    # cuts of [M] carried by Z: a prime u gives the cut {x : <u,c_x> = 1} up to complement
    cuts = set()
    for u in us:
        pat = tuple(parity(u & c) for c in cols)
        comp = tuple(1 - b for b in pat)
        cuts.add(min(pat, comp))
    print(f"t={t} M={M} m={m} r=2^m-1={2**m-1}: 5-subsets {n}, dependent {bad}; "
          f"balance checked on {len(sample)}, failures {unbalanced}; "
          f"cut classes carried {len(cuts)} (incl. trivial) vs 2^(M-1)-1 = {2**(M-1)-1}")
