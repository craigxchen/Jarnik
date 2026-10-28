# Judge check J1 (round 18): exact combinatorial reproduction of numbers the round relies on.
from fractions import Fraction as F
from itertools import combinations
from math import comb, log, ceil
import sympy

def visible_blocks(k):
    # anchored labels 0..k-1, anchor 0; visible T subset of [m]={1..k-1}, 2<=|T|<=m-1
    m=k-1
    return [frozenset(T) for s in range(2,m) for T in combinations(range(1,k),s)]

def boundary_divisors(k):
    # splits S|S^c of {0..k-1}, 2<=|S|<=k-2, represented by side not containing 0
    return visible_blocks(k)

out=[]
for k in range(5,9):
    m=k-1
    vb=len(visible_blocks(k))
    out.append(f"k={k}: visible blocks {vb}, 2^m-m-2={2**m-m-2}, 2^m-m-1={2**m-m-1}, #D_T=2^(k-1)-k-1={2**(k-1)-k-1}")

# K = sum_s (s(n-s)/(n-1) - 2) D_s ; beta has degree 1 on every boundary divisor
for n in (7,8):
    tot=F(0); detail=[]
    for s in range(2,n//2+1):
        cnt = comb(n,s) if 2*s!=n else comb(n,s)//2
        coef = F(s*(n-s),n-1)-2
        tot += coef*cnt; detail.append((s,cnt,coef))
    out.append(f"n={n}: K.beta = {tot}  (s,#D_s,coef) {detail}")

# Theorem E: min # primes r_k = ceil(#visible/(k-3)); p_r; emptiness threshold; minor bound
for k in range(5,9):
    vb=len(visible_blocks(k)); r=ceil(F(vb,k-3)); p=sympy.prime(r)
    for eta in (F(1,10),F(1,2)):
        thr = log(p)/(1+float(eta))   # need e^{(1+eta)w} >= p_r
        ninT = 2**(k-3)-1             # visible T containing fixed i,j>=1
        minor = ninT*(1-float(eta))*thr/log(10)
        out.append(f"k={k} eta={eta}: r={r}, p_r={p}, empty for e^w<{2.718281828459045**thr:.1f}, minors >= 10^{minor:.1f} ({ninT} blocks per minor)")

# Prop 7.5 threshold: 70-182*eta = 0 at eta=5/13 ; 7K.beta check: -56+126 = 70
out.append(f"70-182*5/13 = {70-182*F(5,13)};  -2*28*(1+eta)+(56+70)*(1-eta) at eta=0: {-56+126}")

# forgetful fibres: D_T vertical over Mbar_{0,q} (labels Q={0..q-1}) iff |T∩Q|>=2 and |Q\T|>=2 (T up to complement)
k=8
for q in range(4,8):
    Q=set(range(q)); cnt=0
    for T in visible_blocks(k):
        S=set(T)
        if len(S&Q)>=2 and len(Q-S)>=2: cnt+=1
    out.append(f"forgetful to Mbar_0,{q}: vertical (bounded) contacts {cnt}, formula {(2**(q-1)-q-1)*2**(8-q)}")

# cross-ratio pencil on Mbar_{0,8} via labels {1,2,3,4}: fibres over 0,1,inf = D_T with T∩Q a fixed 2|2 split
Q=[1,2,3,4]
fib={}
for T in visible_blocks(8):
    S=set(T)&set(Q); Sc=set(Q)-set(T)
    if len(S)==2 and len(Sc)==2:
        key=frozenset([frozenset(S),frozenset(Sc)]); fib[key]=fib.get(key,0)+1
Lbeta=list(fib.values())[0]; Fred=sum(fib.values())
out.append(f"cross-ratio: fibres {sorted(fib.values())} -> L.beta={Lbeta}, F_red.beta={Fred}, gain={Lbeta-Fred}")
# expected dims (Kapranov counts) genus 0: (186d+2)-196d
out.append(f"H.beta = #visible T containing a fixed pair i,j>=1 (incl. {{i,j}}) = {sum(1 for T in visible_blocks(8) if {1,2}<=set(T))}; Kapranov: 6H.beta={6*31}, sum(|T|-1)E_T.beta={sum(len(T)-1 for T in visible_blocks(8) if len(T)>=3)} (|T|>=3 terms), K.beta=-6*31+sum(|T|-2)={-6*31+sum(len(T)-2 for T in visible_blocks(8) if len(T)>=3)}")
open('/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round18/judge_checks/j1_out.txt','w').write("\n".join(out)+"\n")
print("\n".join(out))
