"""Exhaustive check of X.2: b points, positive integer distances <= D, every triple degenerate
(largest = sum of other two)  =>  isometric embedding in Z (for b >= 5)."""
from itertools import combinations
import sys

def degenerate(a, b, c):
    x = sorted((a, b, c))
    return x[2] == x[0] + x[1]

def embeds(b, d):
    # x_0 = 0, x_1 = d01; then each x_i determined by d_0i and d_1i up to the choices
    xs = [0]
    def rec(i):
        if i == b:
            return True
        for s in (1, -1):
            x = s * d[(0, i)]
            if all(abs(x - xs[j]) == d[(j, i)] for j in range(i)):
                xs.append(x)
                if rec(i + 1):
                    return True
                xs.pop()
            if d[(0, i)] == 0:
                break
        return False
    return rec(1)

def run(b, D):
    edges = list(combinations(range(b), 2))
    d = {}
    count = fails = 0
    examples = []
    def rec(e):
        nonlocal count, fails
        if e == len(edges):
            count += 1
            if not embeds(b, d):
                fails += 1
                if len(examples) < 3:
                    examples.append(dict(d))
            return
        i, j = edges[e]
        for val in range(1, D + 1):
            d[(i, j)] = val
            ok = True
            # check all triples (k, i, j) with k < j, k != i, whose three edges are assigned
            for k in range(j):
                if k == i:
                    continue
                e1 = (min(i, k), max(i, k)); e2 = (k, j)
                if e1 in d and e2 in d:
                    if not degenerate(d[e1], d[e2], val):
                        ok = False; break
            if ok:
                rec(e + 1)
        del d[(i, j)]
    rec(0)
    return count, fails, examples

if __name__ == "__main__":
    for b, D in [(4, 6), (5, 10), (6, 8), (7, 6)]:
        c, f, ex = run(b, D)
        print(f"b={b}, D={D}: {c} all-degenerate tables, {f} not embeddable in Z; examples {ex[:1]}")
        sys.stdout.flush()
