"""Exhaustive check of the integer prime-field uncertainty lemma (Lemma 5.2 of flipped.md):

For q = 3 mod 4 prime, H = normalised Paley-I matrix of order q+1 (fcommon.paley_I), finite rows x in F_q,
nonconstant columns a in F_q, and every nonzero integer vector n on the nonconstant columns with
|supp n| = s <= (q+1)/2:
      #{x in F_q : (H n)_x != 0}  >=  (q+1)/2 - s.
Also checks the congruence (H n)_x = -P(x) - n_x (mod q), P(X) = sum_a n_a (a-X)^((q-1)/2), used in the proof,
and reports the minimum of #{x in F_q: (Hn)_x != 0} + s over the enumerated vectors (tightness).
"""
import itertools, sys
import numpy as np
from fcommon import paley_I


def run(q, maxsupp, amp):
    H = np.array(paley_I(q), dtype=np.int64)
    Hf = H[1:, 1:]                      # finite rows x (1..q) x nonconstant columns a (1..q)
    d = (q - 1) // 2
    vals = [v for v in range(-amp, amp + 1) if v != 0]
    count = 0
    best = None
    rng = np.random.default_rng(q)
    for s in range(1, maxsupp + 1):
        for supp in itertools.combinations(range(q), s):
            for coeffs in itertools.product(vals, repeat=s):
                n = np.zeros(q, dtype=np.int64)
                n[list(supp)] = coeffs
                Hn = Hf @ n
                nz = int(np.count_nonzero(Hn))
                bound = (q + 1) // 2 - s
                assert nz >= bound, (q, supp, coeffs, nz, bound)
                count += 1
                if best is None or nz + s < best[0]:
                    best = (nz + s, supp, coeffs)
                if count % 997 == 0:          # spot-check the congruence used in the proof
                    for x in range(q):
                        P = sum(int(n[a]) * pow(a - x, d, q) for a in range(q)) % q
                        assert (int(Hn[x]) + P + int(n[x])) % q == 0
    return count, best


def main():
    total = 0
    for q, maxsupp, amp in [(7, 4, 3), (11, 4, 3), (19, 3, 3), (23, 3, 2), (43, 2, 2)]:
        c, best = run(q, maxsupp, amp)
        total += c
        print(f"q={q:3d}: {c:8d} vectors (supp<={maxsupp}, |entries|<={amp}); min (#nonzero rows + s) = {best[0]}"
              f"  [lemma bound (q+1)/2 = {(q+1)//2}]")
    print("TOTAL", total, "PASS")


if __name__ == "__main__":
    main()
