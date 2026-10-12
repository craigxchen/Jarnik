"""Referee check R6: the random two-basis examples of lw.md §7.
Regenerates the profiles with lw_checks/dcore.random_permuted_bases (that function DEFINES 'seed s'),
then recomputes s(c) for the reported light characters with independent integer code, and checks
Lemma 7.1's input (flipped.md Lemma 6.2 transposed) on random sparse c for permuted Paley bases.
Usage: python3 r6_random_bases.py
"""
import sys, random
sys.path.insert(0, '../lw_checks')
from dcore import paley_skew, random_permuted_bases
out = []
def log(s): print(s); out.append(s)
def s_of(c, types): return sum(1 for t in types if sum(a*b for a, b in zip(c, t)) != 0)
cases = [(7, 1, [1, 1, -1, -1, -1, 1, 1, -1]), (11, 1, [1, 0, -1, 0, -1, 1, 0, 0, 0, -1, 0, 1]),
         (11, 2, [0, 1, 0, -1, 0, 1, -1, 0, 0, -1, 0, 1]), (11, 2, [1, 0, 1, -1, -1, 0, 0, 0, 0, 1, -1, 0])]
for q, seed, c in cases:
    K = paley_skew(q); types = sum(random_permuted_bases(K, 2, seed), [])
    assert sum(c) == 0
    log(f"rand2 q={q} seed={seed}: c={c}: s(c) = {s_of(c, types)} of {len(types)} types (M/2 = {(q+1)//2})")
# Lemma 7.1 input: |supp H^T c restricted to nonconstant columns| >= (q+1)/2 - |supp c| for the Paley core,
# hence for every row permutation of it.
random.seed(7)
for q in (19, 23, 43, 47):
    K = paley_skew(q); M = q + 1
    cols = [[K[x][a] for x in range(M)] for a in range(1, M)]
    worst = 10**9
    for _ in range(3000):
        k = random.randint(2, M // 4)
        rows = random.sample(range(M), k); c = [0]*M
        for r in rows: c[r] = random.choice([-2, -1, 1, 2])
        c[rows[0]] -= sum(c)
        if not any(c): continue
        sz = sum(1 for x in c if x)
        nz = s_of(c, cols)
        worst = min(worst, nz - ((q + 1)//2 - sz))
        assert nz >= (q + 1)//2 - sz, (q, c)
    log(f"Paley q={q}: |supp H^T c| - ((q+1)/2 - |supp c|) >= {worst} on 3000 random c with |supp c| <= M/4 (Lemma 7.1 input holds)")
open('r6_random_bases_out.txt', 'w').write('\n'.join(out) + '\n')
