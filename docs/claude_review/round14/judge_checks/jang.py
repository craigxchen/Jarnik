"""Judge check J1 (round 14): per-unit (PU) goodness and exact per-v angle compatibility (ANG)
of the rigtheory joint witnesses (rigtheory_checks/witness*.pkl).

ANG for one v (proved necessary in README sec 5): with N = Norm Z_v, y = 2<v,phi>, (1.2) of sharp.md for s=0..3 says
   2 m2^2/N - 1 <= cos y <= 1 - 2 m0^2/N,   m3^2/N - 1 <= sin y <= 1 - m1^2/N,
so the rational box B_v must meet the unit circle. B_v (convex) meets it iff B_v nonempty and
min_{B_v} |z|^2 <= 1 <= max_{B_v} |z|^2. All arithmetic in Fraction (exact).
Demand: genuine moduli G at all levels (q^vq(D_s(Z_v),q)), other moduli level 1 (the witness has no level-2
collisions, checked by referee vwit.py; re-checked here).
For each column twist t in {0,1}^r the free labels are shifted by (m/2)t (G stays genuine): these are the
associate-twist relatives (referee_rigtheory Prop. T) restricted to the free moduli.
"""
import sys, glob, os, pickle, math, itertools
from fractions import Fraction as Fr
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from jglib import two_sq, genuine, nu, Zv, Ds, VK, vq

def box_meets_circle(lo1, hi1, lo2, hi2):
    if lo1 > hi1 or lo2 > hi2: return False
    def mn(a, b): return Fr(0) if a <= 0 <= b else min(a*a, b*b)
    def mx(a, b): return max(a*a, b*b)
    return mn(lo1, hi1) + mn(lo2, hi2) <= 1 <= mx(lo1, hi1) + mx(lo2, hi2)

def analyse(w, t):
    K, ps, G, q0, lab, l2 = w['K'], w['ps'], w['G'], w['q0'], w['labels'], w['lev2']
    r = len(ps); pis = [two_sq(p) for p in ps]; V = VK(r, K)
    labt = {}
    for q, ell in lab.items():
        m = q-1 if q % 4 == 1 else q+1
        labt[q] = [(x + (m//2)*tj) % m for x, tj in zip(ell, t)]
        # level-2 collision-freeness of the twisted lift
        M2 = m*q; lift = [(x + (M2//2)*tj) % M2 for x, tj in zip(l2[q], t)]
        for v in V:
            if sum(a*b for a, b in zip(v, lift)) % (M2//4) == 0: return None
    pu_viol = 0; ang_fail = []
    for v in V:
        z = Zv(pis, v); N = z[0]**2 + z[1]**2
        ms = []
        for s in range(4):
            d = Ds(z, s); mm = 1
            for q in G: mm *= q**vq(d, q)
            for q, ell in labt.items():
                m = q-1 if q % 4 == 1 else q+1
                if sum(a*b for a, b in zip(v, ell)) % m == s*(m//4): mm *= q
            ms.append(mm)
            if mm*mm > N*(2 if s % 2 else 1): pu_viol += 1
        m0, m1, m2, m3 = ms
        if not box_meets_circle(Fr(2*m2*m2, N) - 1, 1 - Fr(2*m0*m0, N), Fr(m3*m3, N) - 1, 1 - Fr(m1*m1, N)):
            ang_fail.append((v, tuple(ms), N))
    return pu_viol, ang_fail

def main():
    D = sys.argv[1]
    for f in sorted(glob.glob(os.path.join(D, 'witness*.pkl'))):
        w = pickle.load(open(f, 'rb'))
        r = len(w['ps']); res = []
        for t in itertools.product([0, 1], repeat=r):
            a = analyse(w, t)
            if a is None: res.append((t, 'lev2coll')); continue
            pu, af = a
            res.append((t, pu, len(af), af[:1]))
        t0 = res[0]
        good_ang = [x[0] for x in res if x[1] != 'lev2coll' and x[1] == 0 and x[2] == 0]
        print(os.path.basename(f), f"K={w['K']} ps={w['ps']} G={w['G']} q0={w['q0']}",
              f"| t=0: PUviol={t0[1]} ANGfail={t0[2]} ex={t0[3]}",
              f"| twists PU-good={sum(1 for x in res if x[1]==0)}/16 PU&ANG-good={good_ang}", flush=True)

if __name__ == "__main__": main()
