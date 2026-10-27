# Regime numbers for CIRCLE(q), sdc_q, C = 1.  mpmath 50 digits; exact integers via python-flint.
# Outputs: parameters, cap threshold ln T*(tau), and the three rigorous "restricted no-go" lower bounds
# (cube-only, S-grid-only), the slab-reduction mass fraction, and exact Hilbert-function ratios of the cube.
import mpmath as mp, flint, sys
mp.mp.dps = 50
def regime(q):
    n = q + 1; d = q; r = 2 * q
    w = mp.sqrt(8 * d * mp.log(4 * r)); j = int(mp.floor(w / 2)); mp1 = int(mp.ceil(4 * w))
    return n, d, r, w, j, mp1
def lnT(q, tau, mprime, C=1):
    # T* = eps_0/(2 sqrt(m'+1) e^{pi^2/(4q)}),  eps_0 = C d e^{-W/4}/(m'+1),  W = r tau + 1/4 (worst case)
    W = 2 * q * tau + mp.mpf(1) / 4
    return mp.log(C * q) - W / 4 - mp.log(mprime + 1) - mp.log(2) - mp.log(mprime + 1) / 2 - mp.pi**2 / (4 * q)
def tauD(q, c0, C=1):
    n = q + 1; r = 2 * q
    return max((2 / mp.mpf(c0)) * (2 * mp.log(4 * r + 5) + mp.log(1000 / mp.mpf(C))),
               4 * mp.log(4 * r + 5) + 8, 2 * mp.log(1000 * q * n**2 / mp.mpf(C)))
def lncos_sorted_lb(q, i):
    # i-th largest |cos(theta + 2 pi x/q)| over x in F_q is >= cos(pi i/(2q)), for every theta (cap.md Lemma 4.1)
    return mp.log(mp.cos(mp.pi * i / (2 * q)))
out = []
def P(s):
    print(s, flush=True); out.append(s)
for q in [3251, 5479]:
    n, d, r, w, j, m = regime(q)
    P(f"=== q={q}: n={n} d={d} r={r} regime d>=316 ln 4r: {d >= 316*mp.log(4*r)} (316 ln 4r={mp.nstr(316*mp.log(4*r),8)})")
    P(f"w={mp.nstr(w,10)} j={j} m'={m} 2j+1={2*j+1} |P|=(q+1)/2={(q+1)//2}  m'<q: {m<q}  m'/(2j)={mp.nstr(mp.mpf(m)/(2*j),6)}")
    sdP = mp.sqrt(mp.mpf(2) / 3 * (q + 1) / 2)
    P(f"cube sd of k(P) = {mp.nstr(sdP,6)};  j/sd = {mp.nstr(j/sdP,6)}")
    for c0 in [0.184, 0.5, 1.0]:
        t = tauD(q, c0)
        P(f"  tau_D(c0={c0}) = {mp.nstr(t,8)} ; ln T*(tau_D) = {mp.nstr(lnT(q,t,m),10)}")
    tfloor = 2 * mp.log(1000 * q * n**2)
    P(f"  c0-free floor tau >= 2 ln(1000 q n^2) = {mp.nstr(tfloor,8)} ; ln T* there = {mp.nstr(lnT(q,tfloor,m),10)}")
    P(f"  absolute floor tau=4 (Thm 2.3): ln T* = {mp.nstr(lnT(q,4,m),10)}")
    lnfac = mp.loggamma(m + 1)
    # cube-only bound: alpha = 2 on top h coords, 1 on one more if m odd
    h = m // 2
    s = 2 * sum(lncos_sorted_lb(q, i) for i in range(1, h + 1))
    if m % 2: s += lncos_sorted_lb(q, h + 1)
    lnB_T = lnfac - m * mp.log(2) - m * mp.log(d) + s
    P(f"  [cube-only]  ln kappa_V^(T)(u_theta) >= {mp.nstr(lnB_T,10)}   (ln m'!={mp.nstr(lnfac,10)}, sum ln cos={mp.nstr(s,8)})")
    # S-grid-only bound: alpha_a <= 2j; put 2j on top floor(m/2j) coords, remainder on next
    gam = 2 / mp.sqrt(q + 1)
    full, rem = divmod(m, 2 * j)
    s2 = sum((2 * j) * lncos_sorted_lb(q, i) for i in range(1, full + 1)) + (rem * lncos_sorted_lb(q, full + 1) if rem else 0)
    lnB_G = lnfac - m * mp.log(2) + m * mp.log(gam) - m * mp.log(d) + s2
    P(f"  [S-grid-only] ln kappa^(grid)(u_theta) >= {mp.nstr(lnB_G,10)}   (|gamma|=2/sqrt(q+1), alpha={full}x{2*j}+{rem})")
    # slab-reduction mass fraction
    sumc_lb = sum(mp.cos(mp.pi * i / (2 * q)) for i in range(1, q + 1))  # times |gamma|
    lnfrac = (mp.log(mp.binomial(m, 2 * j + 1)) + mp.log(q) + (2 * j + 1) * (0 - mp.log(sumc_lb)))
    P(f"  [slab mass] sum_a |c_a| >= |gamma| * {mp.nstr(sumc_lb,10)};  ln(mass fraction touched by S_a^(2j+1)) <= {mp.nstr(lnfrac,10)}")
    # exact cube Hilbert function H_T(m) = [t^m](1+t+t^2)^q
    poly = flint.fmpz_poly([1, 1, 1]) ** q
    co = poly.coeffs()
    H = lambda e: int(co[e]) if 0 <= e < len(co) else 0
    tot = 3**q
    cum = sum(int(co[e]) for e in range(m))
    P(f"  [Hilbert] ln H_T(m') = {mp.nstr(mp.log(H(m)),10)} ; ln(3^q) = {mp.nstr(q*mp.log(3),10)}")
    P(f"  [Hilbert] ln( sum_(e<m') H_T(e) / 3^q ) = {mp.nstr(mp.log(cum)-q*mp.log(3),10)}  (< ln(1/2) => H_A(m')>0)")
    rat = mp.mpf(2 * q) * H(m - 2 * j - 1) / H(m)
    P(f"  [Hilbert] ln( 2q H_T(m'-2j-1)/H_T(m') ) = {mp.nstr(mp.log(rat),10)}  (relative dim of syntactic slab part of J_gen in degree m')")
    # largest m with cumulative < 3^q/2 (generic tangent degree lower bound)
    acc = 0; mstar = None
    for e in range(len(co)):
        if 2 * (acc + int(co[e])) >= tot: mstar = e; break
        acc += int(co[e])
    P(f"  [Hilbert] H_A(e) > 0 for all e <= {mstar} (cumulative cube count first reaches 3^q/2 at e={mstar}); q={q}")
open('regime_numbers_out.txt', 'w').write('\n'.join(out) + '\n')
