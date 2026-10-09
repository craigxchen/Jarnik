"""Referee numerics for sharp.md Theorem 6.4 (P3-all) parameters and moment bound.

 (1) beta, theta = beta/log X, tau_* as functions of L = log M (A = 1, 2; b = 33): the first L at
     which e^beta > 1 (beta > 0) and theta <= 1/2, and tau_* / (log^2 M / loglog M) against 2A.
 (2) the per-kk moment inequality of step 4:
       pi(24 kk) log(1 + (log_3(24X) + 3.4) e^(2 beta)) + 58 kk e^beta (loglog X + 1)
         <= kk [20 beta + 10 loglog X + 10 + 58 e^beta (loglog X + 1)]
     on a grid of (L, kk) with pi exact (kk <= 10^5) -- i.e. the q <= 24 kk part is <= kk(20 beta + 10 loglog X + 10).
 (3) the exponent identity: log(4er) - theta tau_*/2 + bracket = -4 (exact algebra, evaluated).
"""
import math
def primes_upto(n):
    s = bytearray([1])*(n+1); s[0:2]=b"\x00\x00"
    for i in range(2,int(n**.5)+1):
        if s[i]: s[i*i::i]=bytearray(len(s[i*i::i]))
    return s
S = primes_upto(24*10**5+10)
pi_cum = [0]*(len(S))
c=0
for i in range(len(S)):
    c += S[i]; pi_cum[i]=c
out=[]
for A in (1,2):
    first_beta=None; first_theta=None
    rows=[]
    for L in [10,20,50,100,200,500,1000,2000,5000,10**4,10**5,10**6]:
        M=math.exp(L) if L<700 else None
        logr = math.log(33)+L          # r ~ 33 M
        logX = math.log(3)+A*L
        llX = math.log(logX); llM = math.log(L)
        eb = (math.log(4*math.e)+logr)/(58*(llX+1)*llM)
        beta = math.log(eb)
        theta = beta/logX
        br = (math.log(4*math.e)+logr) + 14 + 20*beta + 10*llX + 58*eb*(llX+1)
        tau = 2*logX/beta*br if beta>0 else float('nan')
        rows.append(f"   L={L:>8}: beta={beta:8.3f} theta={theta:.4f} tau_*/(L^2/logL)={tau/(L*L/llM) if beta>0 else float('nan'):.3f} (2A={2*A})")
        if first_beta is None and beta>0: first_beta=L
    out.append(f"(1) A={A}:")
    out+=rows
# finer search for beta>0 threshold, A=1
for A in (1,2):
    L=5.0
    while True:
        logr=math.log(33)+L; logX=math.log(3)+A*L
        eb=(math.log(4*math.e)+logr)/(58*(math.log(logX)+1)*math.log(L))
        if eb>1: break
        L+=1
    out.append(f"    A={A}: beta > 0 first at log M ~ {L:.0f}")
# (2)
worst=-1e9
for A in (1,2):
    for L in (2000,5000,10**4,10**5,10**6):
        logr=math.log(33)+L; logX=math.log(3)+A*L; llX=math.log(logX); llM=math.log(L)
        eb=(math.log(4*math.e)+logr)/(58*(llX+1)*llM); beta=math.log(eb)
        for kk in (1,2,3,5,10,30,100,1000,10**4,10**5):
            pik = pi_cum[24*kk]
            lhs_small = pik*math.log(1+(logX/math.log(3)+math.log(24)/math.log(3)+3.4)*math.exp(2*beta))
            rhs_small = kk*(20*beta+10*llX+10)
            worst=max(worst, lhs_small/rhs_small)
out.append(f"(2) max over grid of [q<=24kk part]/[kk(20 beta + 10 loglog X + 10)] = {worst:.3f} (<= 1 required)")
# also note the trivial bound E[q^(theta a_q)] <= X^theta = e^beta for q <= 24 kk
out.append("    (the trivial bound E[q^(theta a_q)] <= q^(theta A_q) <= X^theta = e^beta gives the q <= 24kk part <= pi(24kk) beta, tighter)")
# (3)
L=10**4; A=1
logr=math.log(33)+L; logX=math.log(3)+A*L; llX=math.log(logX); llM=math.log(L)
eb=(math.log(4*math.e)+logr)/(58*(llX+1)*llM); beta=math.log(eb); theta=beta/logX
br=(math.log(4*math.e)+logr)+14+20*beta+10*llX+58*eb*(llX+1)
tau=2*logX/beta*br
bracket=20*beta+10*llX+10+58*eb*(llX+1)
out.append(f"(3) exponent per kk at L=1e4: log(4er) - theta tau_*/2 + bracket = {(math.log(4*math.e)+logr) - theta*tau/2 + bracket:.6f} (= -4)")
print("\n".join(out))
