# Referee's independent regime arithmetic for cap.md (no author code imported).
import mpmath as mp, flint, sys
mp.mp.dps = 50
def isprime(n):
    if n<2: return False
    i=2
    while i*i<=n:
        if n%i==0: return False
        i+=1
    return True
# first regime prime q = 3 mod 4 with q >= 316 ln(8q)
first=None
for q in range(7,20000):
    if q%4==3 and isprime(q) and q >= 316*mp.log(8*q):
        first=q; break
print("first regime prime q=3 mod 4:", first)
prev=[q for q in range(3000,first) if q%4==3 and isprime(q)][-1]
print("previous prime 3 mod 4:", prev, " 316 ln(8q)=", mp.nstr(316*mp.log(8*prev),8))
def lncos_sum(q, idx, mult):
    return sum(mult[k]*mp.log(mp.cos(mp.pi*i/(2*q))) for k,i in enumerate(idx))
def run(q, C=1):
    n=q+1; d=q; r=2*q
    w=mp.sqrt(8*d*mp.log(4*r)); j=int(mp.floor(w/2)); m=int(mp.ceil(4*w))
    print(f"\n=== q={q}: d={d} r={r} 316ln4r={mp.nstr(316*mp.log(4*r),8)} regime={d>=316*mp.log(4*r)} w={mp.nstr(w,10)} j={j} m'={m} m'<q={m<q} 4j={4*j} m'/2={m/2}")
    t1=lambda c0: (2/mp.mpf(c0))*(2*mp.log(4*r+5)+mp.log(1000/mp.mpf(C)))
    t2=4*mp.log(4*r+5)+8
    t3=2*mp.log(1000*q*n**2/mp.mpf(C))
    print(" tau_D terms: t2=",mp.nstr(t2,8)," t3=",mp.nstr(t3,8), " t1(0.184)=",mp.nstr(t1(0.184),8))
    def lnT(tau):
        W=r*tau
        lneps0=mp.log(C*d/mp.mpf(m+1))-W/4
        return lneps0-mp.log(2*mp.sqrt(m+1))-mp.pi**2/(4*q)
    floor=max(t2,t3)
    for lab,tau in [("floor",floor),("c0=0.184",max(t1(0.184),t2,t3)),("c0=0.5",max(t1(0.5),t2,t3))]:
        print(f"  ln T* at {lab} (tau={mp.nstr(tau,8)}): {mp.nstr(lnT(tau),10)}")
    LT=lnT(floor)
    lnd=mp.log(d); lnfact=mp.loggamma(m+1)
    # cube mixed difference
    h=m//2  # number of alpha=2
    s=2*sum(mp.log(mp.cos(mp.pi*i/(2*q))) for i in range(1,h+1))+(mp.log(mp.cos(mp.pi*(h+1)/(2*q))) if m%2 else 0)
    cube=lnfact-m*mp.log(2)-m*lnd+s
    print("  cube-only lower bound ln kappa^(T) >=",mp.nstr(cube,10)," (sum ln cos=",mp.nstr(s,8),")  gap to floor T*:",mp.nstr(cube-LT,8))
    # S-grid
    a2=2*j; full=m//a2; rem=m-full*a2
    lng=mp.log(2/mp.sqrt(q+1))
    s2=sum(a2*mp.log(mp.cos(mp.pi*i/(2*q))) for i in range(1,full+1))+(mp.log(mp.cos(mp.pi*(full+1)/(2*q)))*rem if rem else 0)
    # if rem>0 the remaining exponent rem<=2j placed on coordinate full+1 (author: 1 at q=3251)
    grid=lnfact-m*mp.log(2)-m*lnd+m*lng+s2
    print(f"  S-grid: alpha = {full} x {a2} + {rem}; lower bound ln kappa >= {mp.nstr(grid,10)}  gap to floor T*: {mp.nstr(grid-LT,8)}")
    # slab mass fraction
    sc=sum(mp.cos(mp.pi*i/(2*q)) for i in range(1,q+1))
    frac=mp.log(mp.binomial(m,2*j+1))+mp.log(q)-(2*j+1)*mp.log(sc)
    print("  sum_i cos(pi i/2q)=",mp.nstr(sc,10)," ln slab-mass prefactor <=",mp.nstr(frac,10))
    # ln N_0 threshold for survival
    thrN=m*mp.log(4)+2*m*lnd+2*LT
    print("  survival threshold ln N_0 >",mp.nstr(thrN,10))
    print("  ln((m'!)^2 (2/3)^m') =",mp.nstr(2*lnfact+m*mp.log(mp.mpf(2)/3),10))
    print("  per-unit-mass certificate need ln(2^(m'+1) d^m' T*) =",mp.nstr((m+1)*mp.log(2)+m*lnd+LT,10), " cube Lambda value ln(m'!/2^m') =",mp.nstr(lnfact-m*mp.log(2),10))
    # Hilbert function of cube: exact trinomials
    P=flint.fmpz_poly([1,1,1])**q
    HT=lambda e: int(P[e])
    rel=mp.log(2*q)+mp.log(HT(m-2*j-1))-mp.log(HT(m))
    print("  ln(2q H_T(m'-2j-1)/H_T(m')) =",mp.nstr(rel,10))
    cum=sum(HT(e) for e in range(q))
    print("  sum_{e<q} H_T(e) < 3^q/2:", 2*cum < 3**q)
    # exact slab tail: |P|=(q+1)/2 iid uniform{-1,0,1}
    p=(q+1)//2
    Q=flint.fmpz_poly([1,1,1])**p
    tail=sum(int(Q[e]) for e in range(2*p+1) if abs(e-p)>j)
    lt=mp.log(tail)-p*mp.log(3)
    print("  exact ln Pr_T(|k(P)|>j) =",mp.nstr(lt,8)," union bound ln(2q*tail)=",mp.nstr(lt+mp.log(2*q),8), " sd=",mp.nstr(mp.sqrt(mp.mpf(2)*p/3),8), " j/sd=",mp.nstr(j/mp.sqrt(mp.mpf(2)*p/3),6))
    # Prop 6.1
    print("  zonotope bound m <= floor(4jq/(q+1)) =",(4*j*q)//(q+1))
    # hypercontractive threshold
    print("  ln(3^(-3m')/4) =",mp.nstr(-3*m*mp.log(3)-mp.log(4),8))
    # H' = e^{m'W/40}
    for lab,tau in [("floor",floor),("c0=0.184",max(t1(0.184),t2,t3))]:
        print(f"  ln H' at {lab} = {mp.nstr(m*r*tau/40,8)}")
    sys.stdout.flush()
for q in [3251,5479]:
    run(q)
