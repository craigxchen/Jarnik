"""Referee's independent check of the positional digit-sum covering (X.1 inside X.2).

Exact integer arithmetic only.  Does not import anything of Astra's.

Objects (referee's notation):
  H        maximal word mass, B = 64 H^2, g = 2B labels per block, a_hi = B/(4H) = 16H.
  rows     c_t in {0,1}^b, t = 1..T, detecting every nonzero w in Z^b with ||w||_1 <= H.
  U        = 100 g B^T ;  for a prime r, k_r = least k with r^k >= U, D = r^{k_r}.
  label j  of class cl(j): allowed residues mod D are the integers L in the cells
           [D z / B^T, D (z+1) / B^T],  z = sum_t z_t B^{T-t},  z_t in {a_t, a_t+1},
           a_t = 16H * c_t[cl(j)].   For r = 2 only integers of the prescribed parity.
Claims checked:
  (A) attainment: for every block (g labels, arbitrary classes, arbitrary parities when r = 2)
      and every target S mod r^e (e >= k_r, S of the right parity when r = 2) there are allowed
      L_j with sum L_j = S mod r^e (construction: digits least significant first, then
      endpoints, then a correction by a multiple of D on one slot).
  (B) digit intervals: B^{t-1} L_j / D mod 1 lies in [a_t/B, (a_t+2)/B]  (exact rational test).
  (C) words: for every signed word v with ||v||_1 <= H whose class-sum vector w is nonzero,
      D does not divide sum v_j L_j.
"""
import sys, itertools, random
from math import ceil, log2
from fractions import Fraction


def detecting_rows(b, H, T, rng, tries=200):
    """Random search for T rows in {0,1}^b detecting all nonzero w with ||w||_1 <= H.
    Verification is exhaustive over all such w."""
    ws = all_class_vectors(b, H)
    for _ in range(tries):
        rows = [[rng.randrange(2) for _ in range(b)] for _ in range(T)]
        if all(any(sum(r[k] * w[k] for k in range(b)) != 0 for r in rows) for w in ws):
            return rows
    return None


def all_class_vectors(b, H):
    out = set()
    for h in range(1, H + 1):
        for combo in itertools.product(range(b), repeat=h):
            for signs in itertools.product((1, -1), repeat=h):
                w = [0] * b
                for k, sg in zip(combo, signs):
                    w[k] += sg
                if any(w):
                    out.add(tuple(w))
    return sorted(out)


def build_block(H, rows, classes, parities, r, e, S):
    """Return integers L_j (mod r^e) of a block with sum = S mod r^e. Exact arithmetic."""
    B = 64 * H * H
    g = len(classes)
    T = len(rows)
    BT = B ** T
    U = 100 * 2 * B * BT  # g = 2B in the construction; here U uses the nominal g = 2B
    k = 0
    D = 1
    while D < U:
        D *= r
        k += 1
    assert e >= k
    mod = r ** e
    assert D * 1 == r ** k and mod % D == 0
    assert D >= 100 * g * BT
    step = 2 if r == 2 else 1
    if step == 2:
        assert (S - sum(parities)) % 2 == 0
    Sd = S % D
    target = (Sd * BT) // D - g // 2          # target residue of sum z mod B^T
    # digits, least significant first, carrying
    z = [0] * g
    rem = target % BT
    unit = 1
    for t in reversed(range(T)):              # t = T-1 is least significant (weight B^0)
        low = [16 * H * rows[t][classes[j]] for j in range(g)]
        inc = (rem - sum(low)) % B            # number of slots taking a_t + 1
        assert 0 <= inc <= g
        dig = [low[j] + (1 if j < inc else 0) for j in range(g)]
        rem = (rem - sum(dig)) // B
        z = [z[j] + unit * dig[j] for j in range(g)]
        unit *= B
    Z = sum(z)
    assert (Z - target) % BT == 0
    # endpoints (ceil/floor of D z / B^T), parity adjustment
    lo, hi = [], []
    for j in range(g):
        a = -((-D * z[j]) // BT)
        bnd = (D * (z[j] + 1)) // BT
        if step == 2:
            a += (parities[j] - a) % 2
            bnd -= (bnd - parities[j]) % 2
        assert a <= bnd
        lo.append(a)
        hi.append(bnd)
    # choose the translate S* of Sd (mod D) in [sum lo, sum hi]
    base = sum(lo)
    fl = Fraction(D, BT)
    diff = Fraction(Sd) - Fraction(D * Z, BT)
    shift = (diff - (g // 2) * fl) // D          # translate S* = Sd - shift*D
    Sstar = Sd - shift * D
    diff2 = Fraction(Sstar) - Fraction(D * Z, BT)
    # Astra's location claim: S* - D Z / B^T in [floor(g/2), floor(g/2)+1) * D / B^T
    assert (g // 2) * fl <= diff2 < (g // 2 + 1) * fl
    # and this translate lies in the attainable interval [sum lo, sum hi]
    assert base <= Sstar <= sum(hi), (base, Sstar, sum(hi))
    resid = Sstar - base
    assert resid % step == 0
    L = lo[:]
    for j in range(g):
        use = min(resid, hi[j] - L[j])
        use -= use % step
        L[j] += use
        resid -= use
    assert resid == 0
    assert (sum(L) - Sd) % D == 0
    # lift from D to r^e: correct slot 0 by a multiple of D
    corr = (S - sum(L)) % mod
    assert corr % D == 0
    L[0] += corr
    assert (sum(L) - S) % mod == 0
    return L, D, z


def check_intervals(H, rows, classes, L, D):
    B = 64 * H * H
    cnt = 0
    for j, v in enumerate(L):
        for t in range(len(rows)):
            a_t = 16 * H * rows[t][classes[j]]
            res = (B ** t * v) % D             # B^{t} with t starting at 0 = B^{t-1}, t>=1
            # need a_t/B <= res/D <= (a_t+2)/B
            assert a_t * D <= B * res <= (a_t + 2) * D, (j, t)
            cnt += 1
    return cnt


def words_up_to(labels, H):
    """All signed words (as dicts) of mass 1..H on the given labels (multisets with signs)."""
    for h in range(1, H + 1):
        for combo in itertools.combinations_with_replacement(labels, h):
            for signs in itertools.product((1, -1), repeat=h):
                v = {}
                for j, sg in zip(combo, signs):
                    v[j] = v.get(j, 0) + sg
                v = {j: c for j, c in v.items() if c}
                if v:
                    yield v


def main():
    rng = random.Random(20261120)
    total_blocks = total_int = total_words = total_untested_div = 0
    for H, b, mode in [(2, 3, 'onehot'), (2, 4, 'random'), (3, 3, 'onehot'), (3, 3, 'random')]:
        if mode == 'onehot':
            rows = [[int(i == k) for k in range(b)] for i in range(b)]
        else:
            T = ceil(2 * H * log2(2 * b)) + 2
            rows = detecting_rows(b, H, T, rng)
            assert rows is not None
        # exhaustive verification of detection
        ws = all_class_vectors(b, H)
        assert all(any(sum(r[k] * w[k] for k in range(b)) for r in rows) for w in ws)
        B = 64 * H * H
        g = 2 * B
        nblocks = 2
        for r in (2, 3, 5, 7, 13):
            U = 100 * g * B ** len(rows)
            k = 0
            while r ** k < U:
                k += 1
            for e in (k, k + 2):
                for trial in range(4):
                    blocks = []
                    allL, allcl = [], []
                    for blk in range(nblocks):
                        # classes in a block: random (stress) or all equal
                        if trial % 2 == 0:
                            classes = [rng.randrange(b) for _ in range(g)]
                        else:
                            classes = [trial % b] * g
                        par = [rng.randrange(2) for _ in range(g)]
                        S = rng.randrange(r ** e)
                        if r == 2:
                            S += (sum(par) - S) % 2
                            S %= r ** e
                        L, D, z = build_block(H, rows, classes, par, r, e, S)
                        if r == 2:
                            assert all(L[j] % 2 == par[j] for j in range(g))
                        total_int += check_intervals(H, rows, classes, L, D)
                        total_blocks += 1
                        allL += L
                        allcl += classes
                    # word test: a subset of labels covering both blocks
                    n = len(allL)
                    if H == 2:
                        labs = list(range(0, n, 16)) + [1, 2, n - 1]
                    else:
                        labs = list(range(0, n, 64)) + [1, n - 1]
                    labs = sorted(set(labs))
                    for v in words_up_to(labs, H):
                        w = [0] * b
                        for j, c in v.items():
                            w[allcl[j]] += c
                        lam = sum(c * allL[j] for j, c in v.items())
                        if any(w):
                            assert lam % D != 0, v
                            total_words += 1
                        elif lam % D == 0:
                            total_untested_div += 1
        print(f"H={H} b={b} rows={mode} T={len(rows)} detection verified on {len(ws)} class vectors")
    print("blocks built", total_blocks, "digit-interval checks", total_int,
          "tested-word checks", total_words, "untested words with D | Lambda", total_untested_div)
    print("ALL DIGIT CHECKS PASSED")


if __name__ == '__main__':
    main()
