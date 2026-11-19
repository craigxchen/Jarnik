# VII.2(b): exponent-two fraction in the unpadded VII.1 examples; own F_27 = F_3[a]/(a^3 - a - 1)... use a^3 = a + 2 (primitive)
import itertools
p=3
def mul(x,y):  # elements tuples (c0,c1,c2)
    r=[0]*5
    for i in range(3):
        for j in range(3): r[i+j]+=x[i]*y[j]
    # a^3 = a + 2 ; a^4 = a^2 + 2a
    for d in (4,3):
        c=r[d]; r[d]=0; r[d-3]+=2*c; r[d-2]+=c
    return tuple(v%3 for v in r[:3])
F=list(itertools.product(range(3),repeat=3))
one=(1,0,0)
def pw(x,e):
    r=one
    while e:
        if e&1: r=mul(r,x)
        x=mul(x,x); e>>=1
    return r
zero=(0,0,0)
def chi(x):
    if x==zero: return 0
    v=pw(x,13); return 1 if v==one else -1
assert len({pw(x,13) for x in F if x!=zero})==2
def add(x,y): return tuple((a+b)%3 for a,b in zip(x,y))
def sm(c,x): return tuple((c*a)%3 for a in x)
def run(forms,signs):
    labels=list(itertools.product(F,repeat=2))
    n=len(labels); q=27; m=len(forms)
    def ex(v,j):
        a,b=forms[j]; val=add(sm(a,v[0]),sm(b,v[1])); return 1+signs[j]*chi(val)
    g={v:min(ex(v,j) for j in range(m)) for v in labels}
    neg=lambda v:(sm(2,v[0]),sm(2,v[1]))
    G=sum(g.values())
    mixed=0; two=0
    for v in labels:
        es=[ex(v,j)-1 for j in range(m)]
        if 1 in es and -1 in es: mixed+=1
        if g[v]+g[neg(v)]==0: two+=1
    bound=1-2*(q+1)/(m*(q-1))
    print(f"forms={forms} signs={signs} n={n} g={G} mixed={mixed} ({mixed/n:.4f}) exp2-labels={two} exp2-fraction of degN0={two/(n-G):.4f} claimed lower bound={bound:.4f}")
run([(1,0),(0,1),(1,1),(1,2)],[1,1,1,1])
run([(1,0),(0,1)],[1,1])
run([(1,0),(0,1),(1,1)],[1,1,1])
