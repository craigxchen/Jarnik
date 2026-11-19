import random, math, itertools
from fractions import Fraction as F
from sympy import Matrix
from check import *

def mu4(a): return a in [(F(1),F(0)),(F(-1),F(0)),(F(0),F(1)),(F(0),F(-1))]

# Example of 2.3: deleted prime p0 carries column e1; outside 7e1,e2,e3,e4
p0,p1,p2,p3,p4=SPLIT[:5]
cols={p0:(1,[1,0,0,0,0]), p1:(7,[7,0,0,0,0]), p2:(1,[0,1,0,0,0]), p3:(1,[0,0,1,0,0]), p4:(1,[0,0,0,1,0])}
z,y,yc,Delta,W,b,gp=analyse(cols)
print("example distinct:",len(set(z))==5, " I_full=",lattice_index([b[p] for p in cols])[0],
      " I_out=",lattice_index([b[p] for p in cols if p!=p0])[0])

rng=random.Random(7)
nB=nC=0
for trial in range(120):
    S0=SPLIT[:2]
    cols={}
    # S0 primes with random columns
    for p in S0:
        if rng.random()<0.8:
            e=rng.randint(1,3)
            while True:
                a=[rng.randint(0,e) for _ in range(5)]
                if min(a)==0 and max(a)==e: break
            cols[p]=(e,a)
    r_out=rng.randint(0,4)
    outs=rng.sample(SPLIT[2:60], rng.randint(0,6))
    # outside columns: random in a random rank-r_out structure (scaled by mult for big index)
    mult=rng.choice([1,1,2])
    basis=[[rng.randint(-2,2) for _ in range(4)] for _ in range(r_out)]
    for p in outs:
        if r_out==0: break
        co=[rng.randint(-1,1) for _ in range(r_out)]
        vec=[mult*sum(co[i]*basis[i][j] for i in range(r_out)) for j in range(4)]
        # realise vec as a column: a_j - a_5 = vec_j, choose a_5 = -min(0,min vec)
        a5=-min(0,min(vec)); a=[t+a5 for t in vec]+[a5]
        e=max(a)
        if e==0: continue
        cols[p]=(e,a)
    if len(cols)==0: continue
    z=build(cols)
    if len(set(z))<5: continue
    z,y,yc,Delta,W,b,gp=analyse(cols)
    full=Matrix([b[p] for p in cols]).T
    if full.rank()<4: continue
    outc=[b[p] for p in cols if p not in S0]
    rk=Matrix(outc).T.rank() if outc else 0
    if rk==4:
        I,Hn=lattice_index(outc)
        v=short_dual(Hn,I)
        den=1
        for t in v: den=den*t.denominator//math.gcd(den,t.denominator)
        m=den; c=[int(m*t) for t in v]
        u=(F(1),F(0)); extra=(F(1),F(0)); twoh=0
        for p in cols:
            s=sum(c[j]*b[p][j] for j in range(4))
            if p in S0:
                q=round(F(s,m)); r=s-m*q
                assert 2*abs(r)<=m
                u=qmul(u,qpow(gp[p],q)); extra=qmul(extra,qpow(gp[p],4*r)); twoh+=abs(q)*math.log(p)
            else:
                assert s%m==0
                u=qmul(u,qpow(gp[p],s//m)); twoh+=abs(s//m)*math.log(p)
        ycc=(F(1),F(0))
        for j in range(4): ycc=qmul(ycc,qpow(y[j],c[j]))
        assert qpow(ycc,4)==qmul(qpow(u,4*m),extra)
        assert not mu4(qpow(ycc,4)) and qpow(ycc,4)!=(F(1),F(0))
        vinf=float(max(abs(t) for t in v))
        assert twoh <= 4*vinf*W + sum(math.log(p) for p in S0)+1e-9
        nB+=1; print("B",I,m,flush=True)
    else:
        # annihilator via minors
        if rk==0: c=[1,0,0,0]
        else:
            M=Matrix(outc).T
            # pick r independent columns
            cols_idx=[]
            for i in range(M.shape[1]):
                if M[:,cols_idx+[i]].rank()>len(cols_idx): cols_idx.append(i)
                if len(cols_idx)==rk: break
            Mr=M[:,cols_idx]
            rows=None
            for R in itertools.combinations(range(4),rk):
                if Mr[list(R),:].det()!=0: rows=list(R);break
            j0=[j for j in range(4) if j not in rows][0]
            J=sorted(rows+[j0])
            c=[0]*4
            for idx,jj in enumerate(J):
                minor=Mr[[t for t in J if t!=jj],:].det()
                c[jj]=int((-1)**idx*minor)
            assert c[j0]!=0
        assert all(sum(c[j]*bb[j] for j in range(4))==0 for bb in outc)
        assert max(abs(t) for t in c) <= math.factorial(max(rk,1))*max(1,(W/math.log(5)))**rk
        ycc=(F(1),F(0))
        for j in range(4): ycc=qmul(ycc,qpow(y[j],c[j]))
        w4=qpow(ycc,4)
        prod=(F(1),F(0))
        for p in S0:
            if p in cols: prod=qmul(prod,qpow(gp[p],4*sum(c[j]*b[p][j] for j in range(4))))
        assert w4==prod and w4!=(F(1),F(0))
        nC+=1; print("C",rk,flush=True)
print("case B (outside rank 4) checked:",nB," case C (outside rank<4) checked:",nC)
