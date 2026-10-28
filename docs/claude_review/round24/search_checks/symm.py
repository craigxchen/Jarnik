# enumerate symmetric ansaetze and their expected dimension
import itertools, numpy as np, sys
FULL=frozenset(range(1,6))
PL=[frozenset(c) for k in range(1,5) for c in itertools.combinations(range(1,6),k)]  # 30 places
PI={A:i for i,A in enumerate(PL)}
def comp(A): return FULL-A
MAPS={'id':(lambda z:z), '-s':(lambda z:-z), '1/s':(lambda z:1/z), '-1/s':(lambda z:-1/z)}
def compose_maps(m1,a1,m2,a2):
    # returns (function, anti) for m1 after m2 ; anti maps include conj applied first
    f1=MAPS_F[m1] if isinstance(m1,str) else m1
    return None
class Elt:
    def __init__(s,perm,g,anti):  # perm: tuple index->index ; g in MAPS key (real Mobius), anti: conj applied
        s.perm=perm; s.g=g; s.anti=anti
    def f(s,z):
        z=np.conj(z) if s.anti else z
        return MAPS[s.g](z)
    def key(s): return (s.perm,s.g,s.anti)
GMUL={('id','id'):'id'}
names=['id','-s','1/s','-1/s']
def gmul(a,b):  # a after b, as maps
    z=0.37+0.21j; w=MAPS[a](MAPS[b](z))
    for n in names:
        if abs(MAPS[n](z)-w)<1e-12: return n
def mul(x,y):  # x after y
    perm=tuple(x.perm[y.perm[i]] for i in range(30))
    return Elt(perm,gmul(x.g,y.g),x.anti^y.anti)   # real maps commute with conj
def gen_group(gens):
    idp=tuple(range(30)); G={ (idp,'id',False): Elt(idp,'id',False)}
    frontier=list(G.values())
    while frontier:
        new=[]
        for a in frontier:
            for b in gens:
                c=mul(b,a)
                if c.key() not in G: G[c.key()]=c; new.append(c)
        frontier=new
    return list(G.values())
def sym_elt(g,pi,eps):
    perm=[]
    for A in PL:
        B=frozenset(pi[x-1] for x in A)
        if eps<0: B=comp(B)
        perm.append(PI[B])
    return Elt(tuple(perm),g,False)
KAPPA=Elt(tuple(PI[comp(A)] for A in PL),'id',True)
TS=[A for A in PL if 5 not in A]   # 15 rho_T
MEM=[[PI[A] for A in TS if j in A] for j in range(1,5)]
def eqs(rho):  # rho: array of 30 complex values indexed by PL
    out=[]
    for j in range(4):
        c=np.poly(rho[MEM[j]]); out+=list(c[[1,2,3,5,6,7,8]].imag)
    return np.array(out)
def analyse(gens,rng):
    G=gen_group(gens)
    if len(G)>64: return None
    # orbits
    seen=set(); orbs=[]
    for i in range(30):
        if i in seen: continue
        orb={}
        for h in G: orb.setdefault(h.perm[i],h)
        seen|=set(orb); stab=[h for h in G if h.perm[i]==i]
        orbs.append((i,orb,stab))
    # parametrise representatives
    def kind(stab):
        hol=[h for h in stab if not h.anti and h.g!='id']; ant=[h for h in stab if h.anti]
        if any(h.g in ('-s','1/s') for h in hol): return 'bad'
        if any(h.g=='id' and h.anti for h in ant): return 'bad'   # z real
        if hol: # -1/s fixes +-i
            if any(h.anti for h in ant): return 'bad'
            return 'pt'
        if ant:
            gs={h.g for h in ant}
            if len(gs)>1: return 'bad'
            g=gs.pop()
            if g=='-1/s': return 'bad'
            return 'line:'+g   # -s: imaginary axis ; 1/s: unit circle
        return 'free'
    kinds=[kind(st) for _,_,st in orbs]
    if 'bad' in kinds: return ('bad',len(G))
    npar=sum({'free':2,'pt':0}.get(k,1) for k in kinds)
    def build(p):
        rho=np.zeros(30,dtype=complex); q=0
        for (i,orb,st),k in zip(orbs,kinds):
            if k=='free': z=p[q]+1j*p[q+1]; q+=2
            elif k=='pt': z=1j
            elif k=='line:-s': z=1j*np.exp(p[q]); q+=1
            else: z=np.exp(1j*p[q]); q+=1
            for t,h in orb.items(): rho[t]=h.f(z)
        return rho
    # consistency check & jacobian rank at random points
    ranks=[]
    for trial in range(3):
        p=rng.normal(size=npar)
        rho=build(p)
        # check pairing rho[comp]=conj
        ok=all(abs(rho[PI[comp(A)]]-np.conj(rho[PI[A]]))<1e-9 for A in PL)
        if not ok: return ('inconsistent',len(G))
        h=1e-7; J=np.array([(eqs(build(p+h*e))-eqs(build(p-h*e)))/(2*h) for e in np.eye(npar)]).T if npar>0 else np.zeros((28,0))
        sv=np.linalg.svd(J,compute_uv=False) if npar>0 else []
        ranks.append(int(np.sum(np.array(sv)>1e-6*max(1,max(sv) if len(sv) else 1))))
    gauge = 1 if all(h.g in ('id','-s') for h in G) else 0
    return dict(order=len(G),npar=npar,rank=max(ranks),gauge=gauge,expdim=npar-max(ranks)-gauge,kinds=kinds,build=build,orbs=orbs)
if __name__=='__main__':
    rng=np.random.default_rng(0)
    invs=[p for p in itertools.permutations(range(1,6)) if all(p[p[i]-1]==i+1 for i in range(5))]
    # conjugacy class reps of involutions in S5: id, (12), (12)(34)
    reps={'id':(1,2,3,4,5),'(12)':(2,1,3,4,5),'(15)':(5,2,3,4,1),'(12)(34)':(2,1,4,3,5),'(12)(35)':(2,1,5,4,3)}
    res=[]
    for g in ['-s','1/s','-1/s']:
        for pn,pi in reps.items():
            for eps in (1,-1):
                r=analyse([sym_elt(g,pi,eps),KAPPA],rng)
                if isinstance(r,dict): print(g,pn,eps,'order',r['order'],'npar',r['npar'],'rank',r['rank'],'gauge',r['gauge'],'EXPDIM',r['expdim'],flush=True)
                else: print(g,pn,eps,r,flush=True)
