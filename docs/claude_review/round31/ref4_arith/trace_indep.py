# Independent check of the trace construction over F_q, q=2^d, and of the
# violation of (IV.2) in characteristic 2. Own GF(2^d) arithmetic.
import sys
from fractions import Fraction
POLY={3:0b1011,4:0b10011,5:0b100101,6:0b1000011,7:0b10000011}
def run(d):
    q=1<<d; mod=POLY[d]
    def mul(a,b):
        r=0
        while b:
            if b&1: r^=a
            b>>=1; a<<=1
            if a&q: a^=mod
        return r
    def pw(a,e):
        r=1
        e%= (q-1)
        while e:
            if e&1: r=mul(r,a)
            a=mul(a,a); e>>=1
        return r
    def tr(x):
        s=0; y=x
        for _ in range(d): s^=y; y=mul(y,y)
        return s
    assert all(tr(x) in (0,1) for x in range(q))
    def ev(p,x):  # p list of coeffs low->high
        r=0
        for c in reversed(p): r=mul(r,x)^c
        return r
    def deg(p):
        for i in range(len(p)-1,-1,-1):
            if p[i]: return i
        return -1
    f={}
    for a in range(1,q):
        L=[0]*(q//2+1)
        for j in range(d):
            L[1<<j]^=pw(a,(1<<j)-q//2)   # a^{2^j-q/2}
        assert L[q//2]==1
        L[q//2]^=1                       # t^{q/2}-L_a  (char 2: minus=plus)
        assert L[0]==0
        f[a]=L[1:]                       # divide by t
    H={a:[x for x in range(1,q) if tr(mul(a,x))==0] for a in f}
    m=q-1; Lp=q-1   # parent N=t^{q-1}-1, squarefree, L=q-1 layers
    alpha={a:len(H[a]) for a in f}
    assert set(alpha.values())=={q//2-1}
    bad=0; maxdelta=Fraction(-10**9)
    for a in f:
        for b in f:
            if b<=a: continue
            D=[x^y for x,y in zip(f[a],f[b])]
            dg=deg(D); common=set(H[a])&set(H[b])
            assert dg>=0
            if not all(ev(D,x)==0 for x in common): bad+=1
            if dg!=len(common): bad+=1     # quotient by product of common primes is constant
            crit=Fraction(alpha[a]+alpha[b],2)-Fraction(Lp,4)
            maxdelta=max(maxdelta,dg-crit)
    s=[sum(1 for a in f if tr(mul(a,x))==0) for x in range(1,q)]
    lhs=Fraction(3,4)*(m-24)*Lp+sum((Fraction(sh)-Fraction(m,2))**2 for sh in s)
    print(f"q={q} m={m} L={Lp} pairs={m*(m-1)//2} failures={bad} max(deg D - crit)={maxdelta} "
          f"layer sizes={set(s)}  IV.2 LHS={lhs} (RHS=0 since delta=E=0)")
for d in [3,4,5,6,7]: run(d)
