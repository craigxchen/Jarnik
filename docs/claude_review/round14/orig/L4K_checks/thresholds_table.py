# Table: exact free threshold vs pigeonhole and Theorem 1 bounds; Q_7 log-mass share of the free range below Q_K.
import re, math
from sympy import primerange
rows = []
for line in open('free_out.txt'):
    m = re.match(r"K=(\d) r=(\d) (\w+) P=(\[.*?\]) Q_K=(\d+) checked q<=(\d+): non-free q = (\[.*?\]);", line)
    if not m: continue
    K, r, tag, P, QK, hi, nf = int(m[1]), int(m[2]), m[3], eval(m[4]), int(m[5]), int(m[6]), eval(m[7])
    qstar = max(nf)
    nstar = (qstar - (1 if qstar % 4 == 1 else -1)) // 4
    pig = 2 * r + 1                       # forced iff n < |S_1| = 2r+1 (h = 1 for K = 2, 3)
    B = r; AK = K * (2 * (K + 1) ** 2 * B + 2 * K + 1)
    th = lambda lo, hi_: sum(math.log(q) for q in primerange(lo, hi_ + 1) if q % 8 == 7)
    share = th(qstar + 1, QK) / th(3, QK)
    rows.append((K, r, tag, qstar, nstar, pig, AK, QK, share))
print("K r set  q*  n*=max non-free n  pigeonhole: forced if n<  Thm1(b): free if n>  Q_K  theta7 share of (q*,Q_K]")
for row in rows:
    print("%d %d %-5s %4d %4d %4d %5d %7d %.4f" % row)
