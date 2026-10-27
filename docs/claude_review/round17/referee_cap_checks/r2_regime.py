# Referee pass 2: certified (arb ball) recomputation of all regime numbers in cap.md. Own code.
import sys
from flint import arb, ctx, fmpz_poly, fmpz
from sympy import isprime
ctx.prec = 256
out = []
def P(*a):
    s = " ".join(str(x) for x in a); print(s); out.append(s); sys.stdout.flush()
def mid(x, nd=4):
    return f"{float(x.mid()):.{nd}f} (rad {float(x.rad()):.1e})"
pi = arb.pi()
def regime(q):  # n-1 = d = q, r = 2q, 4r = 8q
    return arb(q) >= 316*arb(8*q).log()
# first regime prime = 3 mod 4
first = None
for q in range(3, 6000):
    if q % 4 == 3 and isprime(q) and regime(q):
        first = q; break
P("first prime q=3 mod 4 with q >= 316 ln(8q):", first)
prev = max(p for p in range(3, first) if p % 4 == 3 and isprime(p))
P("previous candidate", prev, "316 ln 8q =", mid(316*arb(8*prev).log()))
# m' < q for all regime q up to 2e5 (and m'/q decreasing beyond)
bad = []
for q in range(first, 200000, 4):
    pass
cnt = 0; worst = 0
for q in range(first, 200001):
    if q % 4 == 3 and isprime(q):
        w = (8*arb(q)*arb(8*q).log()).sqrt(); mp = int(((4*w).ceil()).unique_fmpz())
        cnt += 1; worst = max(worst, mp/q)
        if mp >= q: bad.append(q)
P("regime primes <= 2e5:", cnt, " any with m' >= q:", bad, " max m'/q =", round(worst, 4))
for q in (3251, 5479):
    n = q+1; d = q; r = 2*q
    w = (8*arb(d)*arb(4*r).log()).sqrt()
    j = int((w/2).floor().unique_fmpz()); mp = int((4*w).ceil().unique_fmpz())
    P(f"\n=== q={q}: prime {isprime(q)}, q mod 4 = {q%4}, regime {regime(q)}; w={mid(w)}, j={j}, m'={mp}, m'<q {mp<q}, 4j={4*j}, m'/2={mp/2}, m'/(2j)={mp/(2*j):.4f}")
    t1 = 2*arb(4*r+5).log() + arb(1000).log()   # times 2/c0
    t2 = 4*arb(4*r+5).log() + 8
    t3 = 2*(arb(1000)*q*n*n).log()
    P("tau_D terms: (2/c0)*", mid(t1), "; ", mid(t2), "; floor", mid(t3))
    def lnT(tau, extra=arb(0)):
        W = r*tau + extra
        eps0 = arb(d).log() - W/4 - arb(mp+1).log()
        return eps0 - arb(2).log() - arb(mp+1).log()/2 - pi*pi/(4*q)
    for name, tau in [("floor", t3), ("c0=0.184", max(2*t1/arb("0.184"), t2, t3)), ("c0=0.5", max(2*t1/arb("0.5"), t2, t3))]:
        P(f"  ln T* [{name}, tau={mid(tau,3)}]: W=r tau: {mid(lnT(tau),3)};  W=r tau+1/4: {mid(lnT(tau, arb(0.25)),3)}")
    lT = lnT(t3, arb(0.25))
    # cosines sorted lower bounds
    cs = [ (pi*i/(2*q)).cos() for i in range(1, q+1)]
    S = sum(cs, arb(0))
    P("  sum_{i<=q} cos(pi i/2q) =", mid(S))
    lnfact = arb(mp+1).lgamma()
    # cube: alpha_x <= 2 on top coords
    rem = mp; lnprod = arb(0); i = 0
    while rem > 0:
        a = min(2, rem); lnprod += a*cs[i].log(); rem -= a; i += 1
    cube = lnfact - mp*arb(2).log() - mp*arb(d).log() + lnprod
    P("  Prop3.2 cube: ln kappa^T >=", mid(cube,3), " coords used", i, " gap to floor ln T*:", mid(cube - lT,1))
    gam = arb(2)/arb(q+1).sqrt()
    rem = mp; lnprod = arb(0); i = 0
    while rem > 0:
        a = min(2*j, rem); lnprod += a*cs[i].log(); rem -= a; i += 1
    sg = lnfact - mp*arb(2).log() - mp*arb(d).log() + mp*gam.log() + lnprod
    P("  Prop3.2 S-grid: ln kappa^Y >=", mid(sg,3), " coords used", i, " gap:", mid(sg - lT,1))
    lb = arb(mp+1).lgamma() - arb(2*j+2).lgamma() - arb(mp-2*j).lgamma()
    L31 = lb + arb(q).log() - (2*j+1)*S.log()
    P("  Lemma3.1 prefactor ln:", mid(L31,3))
    surv = 2*lT + 2*mp*arb(d).log() + mp*arb(4).log()
    P("  survival threshold ln N_0 >", mid(surv,2))
    cubeN0 = 2*lnfact + mp*(arb(2)/3).log()
    P("  cube ln N0^T >= ", mid(cubeN0,2), "; fraction ln", mid(surv-cubeN0,1))
    P("  sec6 value per unit l1 mass ln(2^{m'+1} d^{m'} T*) =", mid((mp+1)*arb(2).log() + mp*arb(d).log() + lT,1))
    P("  Prop6.1 bound floor(4jq/(q+1)) =", (4*j*q)//(q+1))
    # exact trinomials
    T = fmpz_poly([1,1,1])**q
    co = T.coeffs()
    def lnz(z): return arb(fmpz(z)).log()
    r35 = arb(2*q).log() + lnz(co[mp-2*j-1]) - lnz(co[mp])
    P("  Prop3.5 ln(2q H_T(m'-2j-1)/H_T(m')) =", mid(r35,3))
    half = sum(co[:q]); P("  sum_{e<q} H_T(e) < 3^q/2:", 2*half < 3**q)
    # slab tail, |P| = (q+1)/2
    p = (q+1)//2
    tp = (fmpz_poly([1,1,1])**p).coeffs()  # index e <-> k(P) = e - p
    tail = 2*sum(tp[p+j+1:])
    lt = lnz(tail) - p*arb(3).log()
    P("  |P|=", p, " sd=", mid((arb(2*p)/3).sqrt(),3), " j/sd=", mid(j/(arb(2*p)/3).sqrt(),3))
    P("  ln Pr_T(|k(P)|>j) =", mid(lt,3), "; ln(2q tail) =", mid(lt + arb(2*q).log(),3))
    P("  Rem6.2: j^2/(2*|X|/2) at |X|=m'/2:", mid(arb(j*j)/(arb(mp)/2),2), " vs (m'/2)ln3+ln4q:", mid(arb(mp)/2*arb(3).log()+arb(4*q).log(),1))
    P("  Rem6.2 hypercontr ln(3^{-3m'}/4) =", mid(-3*mp*arb(3).log()-arb(4).log(),1))
    P("  ln H' = m'W/40 at floor:", mid(mp*r*t3/40,0), " at c0=.184:", mid(mp*r*max(2*t1/arb('0.184'),t2,t3)/40,0))
open(__file__.replace('.py','_out.txt'),'w').write("\n".join(out)+"\n")
