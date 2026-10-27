# regime arithmetic and constants, interval arithmetic (mpmath iv)
from mpmath import iv, mp
import sympy
iv.dps = 40
def regime(q): return iv.mpf(q) >= 316*iv.log(8*iv.mpf(q))   # d = n-1 = q, r = 2q, 4r = 8q
cands = [q for q in sympy.primerange(3, 4000) if q % 4 == 3]
first = None
for q in cands:
    v = iv.mpf(q) - 316*iv.log(8*iv.mpf(q))
    if v.a > 0: first = q; break
    if v.b >= 0: print("undecided interval at", q)
prev = [q for q in cands if q < first][-1]
print(f"first regime prime q = {first}; previous q=3 mod 4 prime {prev}: q - 316 ln 8q = {iv.mpf(prev) - 316*iv.log(8*iv.mpf(prev))}")
print(f"at 3251: q - 316 ln 8q = {iv.mpf(3251) - 316*iv.log(8*iv.mpf(3251))}")
for q in [3251, 3259, 3271, 5479]:
    w = iv.sqrt(8*iv.mpf(q)*iv.log(8*iv.mpf(q))); mprime = int(mp.ceil((4*w).b)); assert int(mp.ceil((4*w).a)) == mprime
    LBw = iv.mpf(q)/2 - 11 - iv.mpf(1.5)*iv.sqrt(q)
    print(f"q={q}: w={w}, j={int(mp.floor(w.a/2))}, m'={mprime}, m'<q {mprime < q}, 2m'<q {2*mprime < q}, "
          f"q/2-11-1.5sqrt q = {LBw}, /r = {LBw/(2*q)}, (q-3)/2/r = {(q-3)/(4*q):.5f}")
# m' < q for all q >= 3251: 4 sqrt(8 q ln 8q) + 1 < q  <=>  for q >= 3251, check derivative-free: h(q) = q - 1 - 4 sqrt(8 q ln 8q) increasing for q >= 3251
q = iv.mpf(3251); print("q-1-4w at 3251:", q - 1 - 4*iv.sqrt(8*q*iv.log(8*q)))
# monotonicity of c(q) = (q/2 - 11 - 1.5 sqrt q)/(2q) = 1/4 - 5.5/q - 0.75/sqrt q: derivative 5.5/q^2 + 0.375 q^{-3/2} > 0 (by hand)
# tau_D with c0 = 1/6 and with c0 = 0.2349, at q = 3251 (C symbolic): print the C-free parts
import math
q = 3251
for c0 in [1/6, 0.2349]:
    print(f"c0={c0}: 2/c0 = {2/c0:.4f}; (2/c0)*2 ln(8q+5) = {(2/c0)*2*math.log(8*q+5):.2f}; 4 ln(8q+5)+8 = {4*math.log(8*q+5)+8:.2f}; "
          f"2 ln(1000 q (q+1)^2) = {2*math.log(1000*q*(q+1)**2):.2f}")
# Prop 2.5 condition 2 leading order: o * sigma_min > r/2  with sigma_min >= c0 r  <=> o > 1/(2 c0)
for c0 in [1/6, 0.23515, 0.2349]:
    print(f"c0={c0}: 1/(2c0) = {1/(2*c0):.4f} -> minimal integer o = {math.floor(1/(2*c0)) + 1}")
