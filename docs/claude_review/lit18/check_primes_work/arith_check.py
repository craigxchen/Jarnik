# Independent arithmetic checks for check_primes.md (adversarial review of primes.md).
import math
def lgbinom2(M): return math.log(M*(M-1)/2)
eps = 0.1
print("(a) Recompute primes.md table at eta=1/2, eps=0.1, lambda_Q ~ log(1/(1-eta))")
eta = 0.5; lam = math.log(1/(1-eta))
for M in [1e3,1e6,1e9,1e12]:
    L=math.log(M); l=math.log(L); k=math.floor((2-eps/4)*L/l)
    logE = lgbinom2(M) + k*math.log(lam) - math.lgamma(k+1)
    ratio = (1-eta)*k*L/lgbinom2(M)
    claimed = (1-eta)*L/(2*l)            # primes.md asymptotic formula
    corrected = (1-eta)*(1-eps/8)*L/l    # (1-eta) k/2 with k ~ (2-eps/4) L/l
    target = lgbinom2(M)-logE
    u=1.0
    while u*math.log(u) < target: u+=0.01
    print(f"  M={M:.0e} k={k:2d} logE_k={logE:6.2f} ratio log Q_k/log binom={ratio:5.2f}  "
          f"primes.md formula (1-eta)L/(2l)={claimed:5.2f}  corrected (1-eta)(1-eps/8)L/l={corrected:5.2f}  u={u:5.2f} log10 D={u*L/math.log(10):6.1f}")
print()
print("(b) Joint modulus vs family size for other window parameters eta (Thm C holds for any fixed eta<1)")
for eta in [0.1,0.25,0.5,0.75,0.9]:
    row=[]
    for M in [1e6,1e9,1e12,1e30]:
        L=math.log(M); l=math.log(L); k=math.floor((2-eps/4)*L/l)
        row.append(f"M={M:.0e}: (1-eta)k/2={(1-eta)*k/2:5.2f}")
    print(f"  eta={eta:4.2f}  "+"  ".join(row))
print("  -> ratio >1 means the k-fold joint modulus exceeds the number of pairs binom(M,2).")
print("     For fixed eta the ratio ~ (1-eta) L/l -> infinity; it is <=1 only if eta >= 1-2/k, i.e. eta -> 1.")
print()
print("(c) What counting/pigeonhole/entropy can FORCE: a collision mod Q is forced only if #classes(Q) ~ Q < M.")
print("    A product of j window primes is >= M^(j(1-eta)), so forced joint collisions need j < 1/(1-eta).")
for eta in [0.25,0.5,0.75,0.9]:
    jmax = math.ceil(round(1/(1-eta),9))-1
    for M in [1e6,1e12]:
        L=math.log(M); l=math.log(L); k=math.floor((2-eps/4)*L/l)
        print(f"  eta={eta:4.2f} M={M:.0e}: max forced joint count j_max={jmax}, heavy pair needs k={k}")
print("    D_A (round-6 design) has c_xy <= floor(18/(1-eta)): bounded, as is j_max. Any input whose only")
print("    forcing is counting (pigeonhole, large sieve, Shearer entropy) is satisfied by D_A at k > 18/(1-eta).")
print()
print("(d) Entropy deficit at moduli Q > M is automatic: for ANY M-point set, log #classes(Q) - H(X mod Q) >= log(Q/M).")
for M in [1e6]:
    for Q in [M**1.0, M**2, M**5]:
        print(f"  M={M:.0e} Q=M^{math.log(Q)/math.log(M):.0f}: forced deficit >= {math.log(Q/M):6.2f} nats, same for fakes and actual clusters")
print()
print("(e) Size of the 'global' integers in the cluster analogue of CMPSY: sum_{x<y} log|B_xy| <= M W/8, W=c0 M log M.")
for c0 in [3.0,4.0,6.0]:
    avg = c0/4   # average log|B_xy| / log M ~ (M W/8)/binom(M,2)/log M = c0/4 * M/(M-1)
    print(f"  c0={c0}: average |B_xy| ~ M^{avg:.2f}; window primes per average pair <= {avg/(1-0.5):.2f} at eta=1/2")
print("    CMPSY-type counting saturates at omega(B) <= log|B|/loglog|B|: this is the two_thirds rate, not better.")
print()
print("(f) Merikoski 2302.11331 range: |B| = X^(1/2-delta), delta<1/10 => |B| > X^0.4. Selection |S|<=M, X=M^h.")
for h in [1.2,2.0,2.5,3.0,5.0,9.0]:
    print(f"  h={h:3.1f}: M = X^{1/h:.3f} -> {'in size range (can pad B)' if 1/h>=0.4 else 'too sparse'}")
print("    h0 = 2c0-3 >= 3 since c0 >= 3 (two_thirds), so the top height is always out of range; in the")
print("    low-height regime P_+(N) < M^(2-delta) the size condition can be met, but the theorem is a count")
print("    over all primes with b in B (main term >= X^0.9/log X >> M for h > 1.12), so it cannot see the selection.")
for h in [1.12,1.5,2.0]:
    X=None
    M=1e6
    print(f"  M=1e6 h={h}: X^0.9/log X = {(M**h)**0.9/math.log(M**h):.3e}  vs M={M:.0e}")
