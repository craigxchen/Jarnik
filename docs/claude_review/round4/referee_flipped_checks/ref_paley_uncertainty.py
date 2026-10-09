"""Referee check of Lemma 6.2 (integer prime-field uncertainty for Paley) of round4/flipped.md.

Claim: q = 3 mod 4 prime, H(x,a) = -chi(a-x) - [a=x] (x,a in F_q).  For every nonzero n in Z^{F_q}
with s = |supp n| <= (q+1)/2:   #{x in F_q : (H n)_x != 0} >= (q+1)/2 - s.

Tests (numpy, exact int64 arithmetic; entries small enough that no overflow occurs):
  A. q=7: ALL nonzero n in {-3..3}^7 (any support; lemma checked where s <= 4).
  B. q=11: all n with support <= 6 and entries in {-2,-1,1,2}.
  C. q in {19,...,83}: random n, supports 1..(q+1)/2, entries including multiples of q and of q^2
     (the 'divide by a power of q' step), and large entries.
Also records min(#nonzero rows + s) per q (author reports (q+1)/2 + 1).
"""
import itertools, random
import numpy as np

def legendre(a, q):
    a %= q
    if a == 0: return 0
    return 1 if pow(a, (q-1)//2, q) == 1 else -1

def Hfin(q):
    return np.array([[-legendre(a - x, q) - (1 if a == x else 0) for a in range(q)] for x in range(q)], dtype=np.int64)

def check_batch(H, N, q):
    """N: (k, q) integer matrix of nonzero vectors. returns min over rows of (#nonzero(Hn) + s) and asserts lemma."""
    T = N @ H.T                      # (k, q): (H n)_x
    nz = (T != 0).sum(axis=1)
    s = (N != 0).sum(axis=1)
    bound = (q + 1)//2
    ok = s > bound
    assert np.all(ok | (nz >= bound - s)), "LEMMA FAILS"
    m = (nz + s)[s <= bound]
    return int(m.min()) if m.size else None, int((s <= bound).sum())

def test_A():
    q = 7; H = Hfin(q)
    vals = np.array(list(itertools.product(range(-3, 4), repeat=q)), dtype=np.int64)
    vals = vals[(vals != 0).any(axis=1)]
    return q, check_batch(H, vals, q)

def test_B():
    q = 11; H = Hfin(q); best = 10**9; cnt = 0
    for s in range(1, 7):
        ent = np.array(list(itertools.product([-2, -1, 1, 2], repeat=s)), dtype=np.int64)
        for S in itertools.combinations(range(q), s):
            N = np.zeros((len(ent), q), dtype=np.int64); N[:, list(S)] = ent
            m, c = check_batch(H, N, q); best = min(best, m); cnt += c
    return q, (best, cnt)

def test_C(rng):
    out = []
    for q in [19, 23, 31, 43, 47, 59, 67, 71, 79, 83]:
        H = Hfin(q); best = 10**9; cnt = 0
        for trial in range(4000):
            s = rng.randint(1, (q+1)//2)
            S = rng.sample(range(q), s)
            n = np.zeros(q, dtype=np.int64)
            mode = trial % 4
            for a in S:
                if mode == 0: v = rng.choice([-1, 1]) * rng.randint(1, 3)
                elif mode == 1: v = rng.choice([-1, 1]) * rng.randint(1, 10**6)
                elif mode == 2: v = rng.choice([-1, 1]) * q * rng.randint(1, 5) if rng.random() < 0.7 else rng.choice([-1, 1])
                else: v = rng.choice([-1, 1]) * rng.choice([1, q, q*q]) * rng.randint(1, 3)
                n[a] = v
            m, c = check_batch(H, n[None, :], q); best = min(best, m); cnt += c
        out.append((q, best, cnt))
    return out

if __name__ == '__main__':
    rng = random.Random(393)
    q, (m, c) = test_A(); print(f"A q={q}: {c} vectors with s<=(q+1)/2 (all of {{-3..3}}^7); min(#nz+s)={m}; (q+1)/2={(q+1)//2}")
    q, (m, c) = test_B(); print(f"B q={q}: {c} vectors (supp<=6, entries +-1,+-2); min(#nz+s)={m}; (q+1)/2={(q+1)//2}")
    for q, m, c in test_C(rng):
        print(f"C q={q}: {c} random vectors (supp up to (q+1)/2, entries incl. multiples of q, q^2, up to 1e6); min(#nz+s)={m}; (q+1)/2={(q+1)//2}")
    print("PASS")
