import sys
sys.setrecursionlimit(100000)
from peel import B, L, phi, gens
from fractions import Fraction as F
for M in range(2, 11):
    T = 40 if M <= 6 else 28
    worst = None; maxratioB = F(0); viol = 0; eqL = 0; tot = 0
    for t in range(1, T + 1):
        for n in range(0, t + 1):
            p = t - n
            b = B(M, p, n); ph = phi(M, p, n)
            tot += 1
            if b > ph: viol += 1; worst = (p, n, b, ph)
            if F(b, t) > maxratioB: maxratioB = F(b, t); arg = (p, n)
    k = (M + 1) // 2
    print(f"M={M}: t<= {T}: violations B>phi: {viol}  max B/t = {maxratioB} at {arg}  conj 2-1/ceil(M/2) = {F(2) - F(1, k)}", flush=True)
