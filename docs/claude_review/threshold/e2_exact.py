# E2: exact (integer) checks of the combinatorial statements.
#  (a) M=2, s=0, D=2: ball measure (1,-2,1) <-> F = -|z1-z2|^2 has order 2; every shell measure at
#      s=0 has order <= 1 (shells S_{0,0}={0}, S_{0,2}={(1,-1),(-1,1)}).  => no decoupling by shell at fixed s.
#  (b) M=3, s=0, D=4: the measure of F_3^2 lives on the ball A_{0,4}, mixes shells t=0,2,4, has exact
#      order 6; the single shells S_{0,4} has reg <= 4 (mod p => rigorous upper bound), S_{0,2} ord 3,
#      S_{0,0} ord 0.  So at (s,D)=(0,4) the ball order 6 exceeds every shell's order.
#  (c) Vandermonde/Ramana measures: exact order = C(M,2) for M = 2..9, ratio C(M,2)/t_min.
#  (d) Generic ball count: |A_{0,D}| exact formula and leading coefficient C(2n,n)/2^n / n!.
#  (e) monotonicity tau(M+1) >= tau(M): embedding k -> (k,0) keeps the order (checked on examples).
import itertools, math
from fractions import Fraction
from taulib import *

def measure_order(meas, maxdeg):
    pts = list(meas.keys()); c = [meas[k] for k in pts]
    return order_exact(pts, c, maxdeg)

def conv(m1, m2):
    out = {}
    for k1, a in m1.items():
        for k2, b in m2.items():
            k = tuple(x + y for x, y in zip(k1, k2))
            out[k] = out.get(k, 0) + a * b
    return {k: v for k, v in out.items() if v != 0}

def perm_sign(p):
    s = 1; p = list(p)
    for i in range(len(p)):
        while p[i] != i:
            j = p[i]; p[i], p[j] = p[j], p[i]; s = -s
    return s

def vandermonde_measure(v):
    M = len(v); out = {}
    for p in itertools.permutations(range(M)):
        k = tuple(v[p[i]] for i in range(M))
        out[k] = out.get(k, 0) + perm_sign(p)
    return {k: x for k, x in out.items() if x}

# A measure c on Z^M <-> homogeneous F of degree D: F = sum_k c_k (u1 v1)^{(D-|k|)/2} u^{k+} v^{k-}.
# Product of forms <-> convolution of measures (degrees add). Order along Y = order of measure
# (per s-slice; all our measures have a single s).

print("(a) M=2")
ball_meas = {(1, -1): 1, (0, 0): -2, (-1, 1): 1}
print("   ball measure (1,-2,1) on A_{0,2}: order =", measure_order(ball_meas, 5), " (F = u1v2+u2v1-2u1v1 = -|z1-z2|^2 on X_2)")
sh = shell(2, 0, 2)
print("   shell S_{0,2} =", sh, " reg_p =", reg_mod(sh, P1)[0], "(order <= 1); S_{0,0} = {0}: order 0")

print("(b) M=3")
F3 = vandermonde_measure((1, 0, -1))
print("   F_3 measure:", F3, " order =", measure_order(F3, 6), " t =", max(sum(abs(x) for x in k) for k in F3))
F3sq = conv(F3, F3)
shells_used = sorted(set(sum(abs(x) for x in k) for k in F3sq))
print("   F_3^2 measure: support size", len(F3sq), " shells used", shells_used, " order =", measure_order(F3sq, 9))
for t in shells_used:
    part = {k: v for k, v in F3sq.items() if sum(abs(x) for x in k) == t}
    print(f"      shell-{t} part: order {measure_order(part, 9)}")
r4, h4 = reg_mod(shell(3, 0, 4), P1); r4b, _ = reg_mod(shell(3, 0, 4), P2)
print("   shell S_{0,4}: reg_p =", r4, r4b, "(=> every measure on S_{0,4} has order <= 4 over Q)")
print("   ball A_{0,4}: reg_p =", reg_mod(ball(3, 0, 4), P1)[0], "(F_3^2 attains it)")

print("(c) Vandermonde measures")
for M in range(2, 7):  # M>=7: order = C(M,2) by the antisymmetrization argument (sum mu V = M! V(v) != 0)
    h = M // 2
    v = tuple(range(-h, M - h))  # M distinct consecutive integers, minimal l1
    meas = vandermonde_measure(v)
    t = sum(abs(x) for x in v); s = sum(v)
    o = measure_order(meas, math.comb(M, 2) + 1)
    print(f"   M={M} v={v} s={s} t={t} |supp|={len(meas)} order={o} (C(M,2)={math.comb(M,2)}) ratio={Fraction(o, t)} = {o/t:.4f}")

print("(d) generic ball count")
def count_exact(M, D):
    # number of k in Z^M with sum 0, |k|_1 <= D (D even): sum over k=|k^+| of N_{2k}
    tot = 0
    for kk in range(0, D // 2 + 1):
        if kk == 0: tot += 1; continue
        for p in range(1, M):
            for q in range(1, M - p + 1):
                tot += math.factorial(M) // (math.factorial(p) * math.factorial(q) * math.factorial(M - p - q)) * math.comb(kk - 1, p - 1) * math.comb(kk - 1, q - 1)
    return tot
for M in (3, 4, 5):
    for D in (2, 4, 6, 8):
        assert count_exact(M, D) == len(ball(M, 0, D)), (M, D)
print("   exact count formula agrees with enumeration (M=3,4,5; D=2..8)")
for M in range(2, 13):
    n = M - 1
    D = 4000
    lead = Fraction(math.comb(2 * n, n), 2 ** n)
    approx = count_exact(M, D) * math.factorial(n) / D ** n
    print(f"   M={M}: n!|A_0,D|/D^n at D={D}: {approx:.4f}  vs C(2n,n)/2^n = {float(lead):.4f};  generic ratio (C(2n,n)/2^n)^(1/n) = {float(lead)**(1/n):.4f};"
          f"  2M^(-1/(M-1)) = {2*M**(-1/n):.4f}")

print("(e) monotonicity: F_3 measure embedded in Z^4, Z^5 keeps order")
for M in (4, 5):
    emb = {k + (0,) * (M - 3): v for k, v in F3.items()}
    print(f"   M={M}: order = {measure_order(emb, 6)}")
