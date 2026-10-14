"""For a witness pickle: at each modulus, which per-column conjugation patterns eps make the label genuine-class."""
import sys, pickle, math, itertools
from rlib14 import *
def eps_classes(ps, q, ell):
    m, L = genuine_labels(ps, q); n = m // 4; r = len(ps)
    c = tuple(x % n for x in ell); out = []
    for eps in itertools.product([1, -1], repeat=r):
        if eps[0] == -1: continue           # global conjugation is the unit u = -1
        Le = [(e * x) % n for e, x in zip(eps, L)]
        if n == 1 or any(tuple((u * x) % n for x in Le) == c for u in range(1, n) if math.gcd(u, n) == 1):
            out.append(eps)
    return out
if __name__ == "__main__":
    w = pickle.load(open(sys.argv[1], "rb"))
    for q in sorted(w['labels']):
        print(q, (q - chi4(q)) // 4, eps_classes(w['ps'], q, w['labels'][q]))
