# Referee J4, check 1 (exact rationals): the averaging constants of J4 step 3.
# (a) q(M) := max_{0<=kap<=M} P(random 5-subset meets both sides of a cut of sizes kap, M-kap)
#          = 1 - min_kap (C(kap,5)+C(M-kap,5))/C(M,5)  <=  15/16 + 15/(16(M-4))   for M >= 9 ?
# (b) reverse Markov: X in [0,1], E X <= q(M) => P(X <= 15/16+eta) >= eta - (q(M)-15/16)... and >= eta/2 when M >= 9+2/eta ?
#     exact: P(X <= a) >= 1 - E X / a is Markov on X itself; the judge uses P(Y >= 1/16-eta) >= E Y - (1/16-eta), Y = 1-X.
# (c) the largest admissible eta for given delta: (15/16+eta)(1+2 delta) < 1, i.e. delta < 1/30 iff eta can be > 0.
from fractions import Fraction as Fr
from math import comb
bad = []
worst = None
for M in range(5, 3001):
    mn = min(Fr(comb(k, 5)+comb(M-k, 5), comb(M, 5)) for k in range(0, M+1) if 0 < k < M)
    q = 1 - mn
    bound = Fr(15, 16) + Fr(15, 16*(M-4))
    if M >= 9 and q > bound: bad.append(M)
    slack = bound - q
    if M >= 9 and (worst is None or slack < worst[1]): worst = (M, slack)
    if M in (5, 6, 7, 8, 9, 10, 20, 50, 100, 1000, 3000):
        print("M=%4d  q(M)=%.6f  15/16+15/(16(M-4))=%.6f  q-15/16=%.3e" % (M, q, bound, q - Fr(15, 16)))
print("violations of (a) for 9<=M<=3000:", bad)
print("smallest slack (M, bound-q):", worst[0], float(worst[1]))
# (b) worst case of reverse Markov with Y = 1 - X, E Y >= 1/16 - 15/(16(M-4)):
viol = []
for eta in [Fr(1, d) for d in (17, 20, 32, 64, 100, 1000)]:
    Mmin = 9 + 2/eta
    import math
    M = math.ceil(Mmin)
    for MM in range(M, M+200):
        EY = Fr(1, 16) - Fr(15, 16*(MM-4))
        lb = EY - (Fr(1, 16) - eta)   # P(Y >= 1/16-eta) >= E Y - (1/16-eta) since Y <= 1
        if lb < eta/2: viol.append((eta, MM))
    print("eta=%s: M_min=%d, lower bound at M_min = %.6f >= eta/2 = %.6f" %
          (eta, M, float(Fr(1, 16) - Fr(15, 16*(M-4)) - (Fr(1, 16)-eta)), float(eta/2)))
print("violations of (b):", viol)
# (c)
for d in (Fr(1, 31), Fr(1, 30)):
    print("delta=%s: (15/16)(1+2delta) = %s" % (d, Fr(15, 16)*(1+2*d)))
