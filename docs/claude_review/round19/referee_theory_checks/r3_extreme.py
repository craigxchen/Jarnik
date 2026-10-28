"""Referee r3: every claimed extreme case + a family of Belyi-type unit composites, exact gains.
Checks gain, gain + n <= Lambda_{k(f)}.beta (Cor 2.5), Thm 3.3 bound gain <= Lambda_k - N - n_1.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as Fr
from r3_lib import *

def run(name, eA, eC, lam):
    A, C = unit_poly(eA), unit_poly(eC)
    e = {}
    for p, v in eA.items(): e[p] = e.get(p, 0) + v
    for p, v in eC.items(): e[p] = e.get(p, 0) - v
    k = kf_unit(e)
    d = gain_data(A, C, Fr(lam))
    N = d['n0'] + d['ninf']
    bound = LAMBDA[k] - N - d['n1']
    ok = d['gain'] <= bound and d['gain'] + d['n'] <= LAMBDA[k]
    print(f"{name:42s} lam={str(Fr(lam)):>8s} k={k} L.beta={d['L']:4d} gain={d['gain']:5d} n={d['n']:3d} N={N:3d} n1={d['n1']:2d} "
          f"R_int(0,1,inf)={d['R_int']}  ThmBound={bound:4d} gain+n={d['gain']+d['n']:3d}<=Lam_k={LAMBDA[k]} {'OK' if ok else 'VIOLATION'}  B={d['Bfactors']}", flush=True)
    return d

# z = (x0-x1)/(x0-x2), 1-z = (x1-x2)/(x0-x2)
z = (0, 1); w = (1, 2); den = (0, 2)
for D in range(1, 6):
    run(f"z^{D}", {z: D}, {den: D}, 1)
run("4z(1-z)", {z: 1, w: 1}, {den: 2}, 4)
run("r18 referee: -1/4 (x0-x1)^2/((x0-x2)(x1-x2))", {(0, 1): 2}, {(0, 2): 1, (1, 2): 1}, Fr(-1, 4))
# lam z^a (1-z)^b with critical value 1: lam = (a+b)^(a+b)/(a^a b^b)
for a in range(1, 5):
    for b in range(1, 5):
        if a + b <= 5:
            run(f"z^{a}(1-z)^{b}", {z: a, w: b}, {den: a + b}, Fr((a + b) ** (a + b), a ** a * b ** b))
# lam z^a / (1-z)^b , a>b>=1: f = lam (x0-x1)^a (x0-x2)^(b-a) / (x1-x2)^b ; critical point z = a/(a-b)
for a in range(2, 6):
    for b in range(1, a):
        zc = Fr(a, a - b)
        val = zc ** a / (1 - zc) ** b
        eA = {z: a}; eC = {w: b}
        if a - b > 0:
            eC[den] = a - b
        run(f"z^{a}/(1-z)^{b}", eA, eC, 1 / val)
# Ptolemy-type and anchored Ptolemy (round 18 sanity: -32)
run("Ptolemy (x0-x1)(x2-x3)/((x0-x2)(x1-x3))", {(0, 1): 1, (2, 3): 1}, {(0, 2): 1, (1, 3): 1}, 1)
run("Ptolemy lam=-1", {(0, 1): 1, (2, 3): 1}, {(0, 2): 1, (1, 3): 1}, -1)
run("linear (x0-x1)/(x2-x3)", {(0, 1): 1}, {(2, 3): 1}, 1)
run("8-point unit (x0-x1)(x2-x3)(x4-x5)/((x0-x2)(x1-x4)(x3-x6))", {(0, 1): 1, (2, 3): 1, (4, 5): 1}, {(0, 2): 1, (1, 4): 1, (3, 6): 1}, 1)
run("star unit k=8 (min N found by r2)", {(0, 4): 1, (0, 5): 1, (0, 6): 1}, {(0, 1): 1, (0, 2): 1, (0, 3): 1}, 1)
