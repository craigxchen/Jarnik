#!/usr/bin/env python3
"""Referee checks R8-R10.
R8 Corollary 4.5 parameter chain.  For a grid of eta in (0,1/4], C >= 1, beta in (1/3, 1/2-eta], and
   log n equal to (1+1e-9) times the hypothesis threshold (the worst case), with m = ceil(2/eta),
   w = ceil(m/beta), omega = w+1, check condition (*) of Lemma 4.4 in logarithmic form:
     omega*log(omega) + m(m+1)/2*log n + omega(omega-1)/2*log X < beta*m*omega*log n,
   where log X = log C + (3 beta/2 - 1/2) log n.  Also check w <= 6/eta + 4 and t = w - m >= 0.
   (Floating point; margins reported.)
R9 Matveev constant: c(k) = min((1/2)(e k/2)^2 30^(k+3) k^3.5, 2^(6k+20)) >= 7e5 for k = 1..60,
   and the exponent comparison 5.9e6 * sum log p >= (1/2) * a_max * sum log p iff a_max <= 1.18e7.
R10 Prop 4.1 bookkeeping: with L = {(a,c) in Z^M x Z : sum a + 2c = 0}, sigma(a,c) = (-a,-c) and
   span{psi_j} = {(2a, -sum a)} contains 2L with index 2 (so 'span 2L' is a misstatement, harmless).
"""
import math
# R8
worst = None; fails = 0; checked = 0
for eta in [0.25, 0.2, 0.1, 0.05, 0.02, 0.01, 0.005, 0.001]:
    m = math.ceil(2/eta)
    for C in [1, 1.5, 10, 1e3, 1e6]:
        thr = (16/(5*eta))*(math.log(C) + 2*math.log(6/eta + 5))
        logn = thr*(1 + 1e-9)
        for i in range(1, 400):
            beta = 1/3 + (0.5 - eta - 1/3)*i/399
            if beta <= 1/3: continue
            w = math.ceil(m/beta); om = w + 1
            assert w - m >= 0 and w <= 6/eta + 4 + 1e-12
            logX = math.log(C) + (1.5*beta - 0.5)*logn
            assert logX >= 0
            lhs = om*math.log(om) + m*(m+1)/2*logn + om*(om-1)/2*logX
            rhs = beta*m*om*logn
            checked += 1
            marg = (rhs - lhs)/rhs
            if lhs >= rhs: fails += 1
            if worst is None or marg < worst[0]: worst = (marg, eta, C, beta)
print(f"R8 Cor 4.5 chain: {checked} (eta,C,beta) cases at the threshold log n: failures {fails}; "
      f"smallest relative margin {worst[0]:.4f} at eta={worst[1]}, C={worst[2]}, beta={worst[3]:.4f}")

# R9
cs = []
for k in range(1, 61):
    c = min(0.5*(math.e*k/2)**2*30**(k+3)*k**3.5, 2.0**(6*k+20))
    cs.append(c)
print(f"R9 Matveev c(k): min over k=1..60 = {min(cs):.4e} at k={cs.index(min(cs))+1} (claim >= 7e5)")
E = 7e5*4*(1+math.log(2))*(math.pi/2)*0.8
print(f"   exponent constant 7e5*4*(1+log2)*(pi/2)*0.8 = {E:.4e}; beats (1/2) a_max iff a_max <= {2*E:.3e}")
# product >= 0.8 sum for x_j >= 1.6: minimum of prod/sum at x_j = 1.6
print("   min over k of 1.6^(k-1)/k =", min(1.6**(k-1)/k for k in range(1, 30)))

# R10
import itertools
M = 3
L = [(a, (-sum(a))//2) for a in itertools.product(range(-3, 4), repeat=M) if sum(a) % 2 == 0]
span_psi = set()
for coef in itertools.product(range(-3, 4), repeat=M):
    span_psi.add((tuple(2*x for x in coef), -sum(coef)))
twoL = set((tuple(2*x for x in a), 2*c) for a, c in L)
psi1 = ((2, 0, 0), -1)
print(f"R10 psi_1 = {psi1} in 2L? {psi1 in twoL};  2L subset of span(psi) (on a box)? "
      f"{all(v in span_psi for v in twoL if max(abs(t) for t in v[0]) <= 6)}")
