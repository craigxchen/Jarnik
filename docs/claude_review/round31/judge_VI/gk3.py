# formal g for k=3, q=27 (own F_27 = F_3[y]/(y^3-y-1)); elements as 3-tuples
import itertools
els=list(itertools.product(range(3),repeat=3))
def mul(a,b):
    r=[0]*5
    for i in range(3):
        for j in range(3): r[i+j]+=a[i]*b[j]
    for d in (4,3):  # y^3=y+1
        c=r[d]; r[d]=0; r[d-3]+=c; r[d-2]+=c
    return tuple(x%3 for x in r[:3])
sq=set(mul(a,a) for a in els if a!=(0,0,0))
assert len(sq)==13
def chi(a): return 0 if a==(0,0,0) else (1 if a in sq else -1)
def lam(f,x): return tuple(sum(f[b]*x[b][c] for b in range(3))%3 for c in range(3))
def g(forms,signs):
    tot=0
    for x in itertools.product(els,repeat=3):
        e=[]
        for f,s in zip(forms,signs):
            l=lam(f,x)
            if s==-1: l=tuple((-c)%3 for c in l)
            e.append(1+chi(l))
        tot+=min(e)
    return tot
A=[(1,0,0),(0,1,0),(0,0,1),(1,1,0),(1,0,1)]
print("Astra 5 forms x,y,z,x+y,x+z: g=",g(A,[1]*5)," h=",729)
R=[(0,0,1),(0,1,0),(0,1,1),(1,0,0),(1,1,2)]   # x2,x1,x1+x2,x0,x0+x1-x2 (coords x0,x1,x2)
print("referee forms, signs (+,+,-,+,+): g=",g(R,[1,1,-1,1,1]))
