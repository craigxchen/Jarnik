"""Two-class ratios for subsets of recurrence sequences (Fibonacci, Pell, Lucas, 2^n, squares,
x^2-type square roots).  For each sequence and each B, the best log q2 / log diam over all
B-subsets of a window of consecutive terms (exact q2)."""
import itertools
from math import log
from tc import q2

def fib(N):
    a, b = 0, 1; out = []
    for _ in range(N): out.append(a); a, b = b, a + b
    return out
def pell(N):
    a, b = 0, 1; out = []
    for _ in range(N): out.append(a); a, b = b, 2 * b + a
    return out
seqs = {
    'Fibonacci': fib(80), 'Pell': pell(60), 'pow2': [2 ** n for n in range(80)],
    'pow3': [3 ** n for n in range(60)], 'squares': [n * n for n in range(1, 400)],
}
for name, s in seqs.items():
    for B in (5, 6, 7):
        best = (0,)
        if name == 'squares':
            windows = [list(range(x, x + 14)) for x in (10, 100, 300)]
        else:
            windows = [list(range(i, i + 13)) for i in (20, 35)]
        for W in windows:
            vals = [s[i] for i in W if i < len(s)]
            for sub in itertools.combinations(vals, B):
                L = max(sub) - min(sub)
                if L < 2: continue
                r = log(q2(list(sub))) / log(L)
                if r > best[0]: best = (r, sub)
        print(f"{name:10s} B={B}: best ratio {best[0]:.4f}  example {best[1] if len(str(best[1]))<160 else '...'}")
