# C1: exhaustive check of Lemma 5.1 (odd encoding) for all N = 0 mod 4, 16 <= N <= 1200, b in {0,1};
# and that it FAILS to exist for N in {4, 8, 12} where needed (counting).
from consD import check_lemma51, lemma51_map
bad = []
for N in range(16, 1201, 4):
    for b in (0, 1):
        if not check_lemma51(N, b):
            bad.append((N, b))
print("Lemma 5.1 exhaustive N=16..1200:", "ALL OK" if not bad else bad)
for N in (4, 8, 12):
    for b in (0, 1):
        try:
            lemma51_map(N, b); print(N, b, "map exists (counting)")
        except AssertionError:
            print(N, b, "no orbit injection (counting fails)")
