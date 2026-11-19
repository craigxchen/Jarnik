# Parity lemma check (char != 2): if A cap -A = empty and prod(X+a)-prod(X-a) is constant in X, then |A| is odd.
# (i) integers: all A subset {+-1..+-7} with A cap -A empty, |A|<=6;  (ii) F_3[t]_{<3}: all A with |A|<=4.
import itertools, sympy as sp
X=sp.symbols('X')
cnt={}
reps=range(1,8)
for k in range(1,7):
    for absA in itertools.combinations(reps,k):
        for sg in itertools.product((1,-1),repeat=k):
            A=[a*s for a,s in zip(absA,sg)]
            if A[0]<0: continue  # A and -A symmetric
            d=sp.Poly(sp.expand(sp.prod([X+a for a in A])-sp.prod([X-a for a in A])),X)
            if d.degree()<=0:
                cnt.setdefault(k,[]).append(A)
print("integer A with constant difference, by size:",{k:len(v) for k,v in cnt.items()}, "examples:",{k:v[:3] for k,v in cnt.items()})
# F_3[t]: elements as tuples of 3 coefficients; polynomial in X over F_3[t]: use sympy with modulus via t symbol
t=sp.symbols('t')
els=[c for c in itertools.product(range(3),repeat=3) if any(c)]
def poly(c): return sum(ci*t**i for i,ci in enumerate(c))
pairs=[];seen=set()
for c in els:
    if c in seen: continue
    n=tuple((-x)%3 for x in c); seen|={c,n}; pairs.append((c,n))
found={}
for k in range(1,5):
    for P in itertools.combinations(pairs,k):
        for sg in itertools.product((0,1),repeat=k):
            if sg[0]==1: continue
            A=[p[s] for p,s in zip(P,sg)]
            e=sp.expand(sp.prod([X+poly(a) for a in A])-sp.prod([X-poly(a) for a in A]))
            d=sp.Poly(e,X,t,modulus=3)
            degX=max([m[0] for m in d.monoms()],default=-1) if not d.is_zero else -1
            if degX<=0: found.setdefault(k,0); found[k]+=1
print("F_3[t]_{<3}: number of A (up to sign) with X-constant difference, by size:",found)
