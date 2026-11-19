"""Bezoutian octic test over F_p(i), including all geometric pivot values.

The gauge fixes one of the first star's three outer roots at i.  A witness
anchor tests all independent matchings of the last star's roots and barred
roots. The resultant gcd is in F_{p^2}[b], not just a scan of pivot values.
"""
from itertools import permutations, combinations

class Field:
    def __init__(self, p):
        assert p % 4 == 3
        self.p, self.q = p, p*p
        q=self.q
        self.plus=[[0]*q for _ in range(q)]
        self.minus=[[0]*q for _ in range(q)]
        self.times=[[0]*q for _ in range(q)]
        self.bar=[0]*q
        self.inv=[0]*q
        for z in range(q):
            a,b=z%p,z//p
            self.bar[z]=a+p*(-b%p)
            for w in range(q):
                x,y=w%p,w//p
                self.plus[z][w]=(a+x)%p+p*((b+y)%p)
                self.minus[z][w]=(a-x)%p+p*((b-y)%p)
                c=(a*x-b*y)%p+p*((a*y+b*x)%p)
                self.times[z][w]=c
                if c==1:self.inv[z]=w
    def powers(self,z):
        x=1; m=[]
        for _ in range(5):
            x=self.times[x][z];m.append(x//self.p)
        return m

def shave(f):
    while f and f[-1]==0:f.pop()
    return f

def mul(F,f,g):
    if not f or not g:return []
    a=F.plus;m=F.times
    h=[0]*(len(f)+len(g)-1)
    for i,x in enumerate(f):
        for j,y in enumerate(g):h[i+j]=a[h[i+j]][m[x][y]]
    return shave(h)

def evaluate(F,f,x):
    h=0
    for v in reversed(f):h=F.plus[F.times[h][x]][v]
    return h

def gcd(F,f,g):
    a=F.minus;m=F.times
    while g:
        if not f: f,g=g,f;continue
        t=f[:]
        while len(t)>=len(g):
            k=len(t)-len(g);v=m[t[-1]][F.inv[g[-1]]]
            for j,z in enumerate(g):t[k+j]=a[t[k+j]][m[v][z]]
            shave(t)
        f,g=g,t
    return [m[z][F.inv[f[-1]]] for z in f] if f else []

def cubic(F,roots):
    f=[1]
    for w in roots:f=mul(F,f,[F.minus[0][w],1])
    return f

def interpolate(F,values):
    bases=[1]; f=[]
    for t,y in enumerate(values):
        k=F.times[F.minus[y][evaluate(F,f,t)]][F.inv[evaluate(F,bases,t)]]
        f.extend([0]*(len(bases)-len(f)))
        for j,x in enumerate(bases):f[j]=F.plus[f[j]][F.times[k][x]]
        bases=mul(F,bases,[F.minus[0][t],1])
    return shave(f)

def octics(F,u,v,U,V):
    """u,v = root stars, U,V their independent barred-root lists."""
    plus,minus,times=F.plus,F.minus,F.times
    pencils=[]
    for j in range(3):
        L=cubic(F,[V[j],u[(j+1)%3],u[(j+2)%3]])
        R=cubic(F,[U[j],v[(j+1)%3],v[(j+2)%3]])
        pencils.append((L,R))
    table=[[] for _ in range(3)]
    for b in range(9):
        quads=[]
        for L,R in pencils:
            aa=evaluate(F,R,b);bb=evaluate(F,L,b)
            c2=minus[aa][bb]
            c1=plus[times[b][c2]][minus[times[aa][L[2]]][times[bb][R[2]]]]
            c0=plus[times[b][c1]][minus[times[aa][L[1]]][times[bb][R[1]]]]
            quads.append((c2,c1,c0))
        for j in range(3):
            (a,b1,c),(A,B,C)=quads[j],quads[(j+1)%3]
            t=minus[times[a][C]][times[A][c]]
            z=minus[times[a][B]][times[A][b1]]
            w=minus[times[b1][C]][times[B][c]]
            table[j].append(minus[times[t][t]][times[z][w]])
    return [interpolate(F,values) for values in table]

def candidates(F):
    p=F.p; power=[F.powers(x) for x in range(F.q)]
    def sig(star):
        return tuple(map(lambda v:sum(v)%p,zip(*[power[x] for x in star])))
    def net(a):return len(set(a)|{F.bar[x] for x in a})==2*len(a)
    good={}
    for a,b in combinations(range(p,F.q),2):
        u=(p,a,b)
        if net(u):good.setdefault(sig(u),[]).append(u)
    for v in combinations(range(p,F.q),3):
        for u in good.get(sig(v),[]):
            if net(u+v):yield u,v

def audit(F,u,v):
    U=[F.bar[x] for x in u]
    count={}
    degs=set()
    for vv in permutations(v):
        for VV in permutations([F.bar[x] for x in v]):
            e=octics(F,u,vv,U,VV)
            g=gcd(F,e[0],e[1]); d1=len(g)-1
            g=gcd(F,g,e[2]); d2=len(g)-1
            inf=all(len(f)<9 for f in e)
            ds=(d1,d2,inf)
            count[ds]=count.get(ds,0)+1
            degs.add(tuple(len(f)-1 for f in e))
    return count,degs

if __name__=='__main__':
    for p in (7,11,19):
        F=Field(p)
        for j,(u,v) in enumerate(candidates(F)):
            ans,degs=audit(F,u,v)
            if all(d2==0 and not inf for (d1,d2,inf) in ans):
                print('p',p,'anchor ordinal',j,'u',u,'v',v,'gcd degrees of (E0,E1) and all, infinity:',ans,'octic degrees',sorted(degs),flush=True)
                break
            if j==25:print('no witness in first 26');break
