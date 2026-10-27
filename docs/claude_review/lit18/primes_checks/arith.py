# Arithmetic for primes.md: what H(eta,delta) of p3lower.md Thm C demands, versus what
# Kubilius/Poisson-model, level-of-distribution and robust-ensemble results can supply.
import math
eta, eps = 0.5, 0.1
lam = math.log(1/(1-eta))          # lambda_Q ~ log(1/(1-eta)) when no prime of N in the window
print("eta=%.2f eps=%.2f lambda_Q~%.4f" % (eta, eps, lam))
print("%8s %6s %6s %4s %12s %14s %14s %10s %10s %12s" % ("M","L","l","k","logE_k","log10Q_k","log10pairs","ratio","u_needed","log10D_need"))
for M in [1e2,1e3,1e4,1e6,1e9,1e12]:
    L = math.log(M); l = math.log(L)
    k = math.floor((2-eps/4)*L/l)
    # E_k ~ binom(M,2) lam^k/k!  (Lemma C.1 upper/lower shape)
    logE = math.log(M*(M-1)/2) + k*math.log(lam) - math.lgamma(k+1)
    log10Q = (1-eta)*k*L/math.log(10)        # product of k window primes >= M^((1-eta)k)
    log10pairs = math.log10(M*(M-1)/2)
    # Kubilius-type TV error ~ u^{-u}; need u^{-u} <= P(event) ~ E_k/binom(M,2)
    target = -(logE - math.log(M*(M-1)/2))   # = -log P(c_xy>=k) >0
    u = 1.0
    while u*math.log(u) < target: u += 0.01
    log10D = u*L/math.log(10)                # level D = y^u with y = M
    print("%8.0e %6.2f %6.2f %4d %12.2f %14.2f %14.2f %10.2f %10.2f %12.2f" % (M,L,l,k,logE,log10Q,log10pairs,log10Q/log10pairs,u,log10D))
print()
# Size of a heavy pair integer: |B_xy| >= M^((1-eta)k) vs L1 bound M^b
for M in [1e3,1e6,1e9]:
    L=math.log(M); l=math.log(L); k=math.floor((2-eps/4)*L/l)
    print("M=%.0e: heavy pair needs |B_xy| >= M^%.2f; L1 allows |B_xy|<=M^b, so H & L1 clash once (1-eta)k > b" % (M,(1-eta)*k))
print()
# Graft visibility for robust inputs: the graft moves 2s<=2M primes; a robust error term c X^gamma at
# height X=M^h hides them iff h*gamma>1.  Regime (1.4): h0 = 2c0-3+o(1).
for c0 in [3.2,4,6]:
    h0 = 2*c0-3
    for name,gamma in [("GRH sqrt",0.5),("Huang-Liu-Rudnick sector count X^(2/5)",0.4),("Harman-Lewis X^0.619 area",0.619)]:
        print("c0=%.1f h0=%.1f %-40s h0*gamma=%.2f -> %s" % (c0,h0,name,h0*gamma,"invisible" if h0*gamma>1 else "visible in size"))
print()
# Merikoski (Compositio 2025): arbitrary B with |B| = X^(1/2-delta), delta<1/10 => |B| >= X^0.4.
for h in [3,5,9]:
    print("primes <= X=M^%d: cluster selection |S|<=M = X^%.3f  vs Merikoski adversarial range |B|>=X^0.4 -> %s" % (h,1/h,"too sparse" if 1/h<0.4 else "in range"))
