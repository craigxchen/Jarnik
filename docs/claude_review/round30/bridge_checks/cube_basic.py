# Cube cluster (ff Thm C): exact structure, pairwise degrees, gamma = coeff of t^4 in z_i - z_j,
# 5-subsets: gcd degree, contact, alpha = (-K.C)/mu0, cut profile, the 15 quadric cones of round 27.
import sys, itertools
sys.path.insert(0,'/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round30/bridge_checks')
from glib import *

def cube(y):
    assert sum(y)==0 and y[0]*y[1]*y[2]*y[3]==1
    eps=[1,-1,-1,1]
    r=[G(Fr(e)/v, v) for e,v in zip(eps,y)]
    roots=[]                       # coordinate order: +r0,-r0,+r1,-r1,...
    for x in r: roots+=[x,-x]
    sets={'0':[], '1':list(range(8)), 'A':[0,1,2,3],'Ac':[4,5,6,7],'B':[0,1,4,5],'Bc':[2,3,6,7],
          'P':[0,2,4,6],'Pc':[1,3,5,7]}
    pts={}
    for nm,F in sets.items():
        rr=[roots[k].conj() if k in F else roots[k] for k in range(8)]
        pts[nm]=(prod_lin(rr), [1 if k in F else 0 for k in range(8)])
    return roots,pts

def analyse(y,verbose=True):
    roots,pts=cube(y)
    names=list(pts)
    N=pmul(pts['0'][0],pconj(pts['0'][0]))
    for nm in names:
        z=pts[nm][0]; assert psub(pmul(z,pconj(z)),N)==[G(0)] or deg(psub(pmul(z,pconj(z)),N))<0
    D=8; gam={}; degs={}
    for a,b in itertools.combinations(names,2):
        d=psub(pts[a][0],pts[b][0]); degs[(a,b)]=deg(d)
        gam[(a,b)]=d[4] if len(d)>4 else G(0)
    if verbose:
        print("y=",y); print("pair degrees:",sorted(set(degs.values())), {k:v for k,v in degs.items() if v<4})
        print("max |gamma|^2:",max(g.norm() for g in gam.values()))
    # 5-subsets
    res=[]
    for S in itertools.combinations(names,5):
        vecs=[pts[s][1] for s in S]
        g=sum(1 for k in range(8) if len(set(v[k] for v in vecs))==1)
        mu0=D-max(degs[tuple(sorted((a,b),key=names.index))] for a,b in itertools.combinations(S,2))
        BC=2*(D-g); alpha=Fr(BC,mu0)
        # cut profile: each coordinate k gives the cut {j: v_j[k]=1}; count weight 2 (layer of degree 2)
        cuts={}
        for k in range(8):
            A=frozenset(j for j in range(5) if vecs[j][k]==1)
            Am=min(A,frozenset(range(5))-A,key=lambda s:(len(s),sorted(s)))
            if 0<len(A)<5: cuts[Am]=cuts.get(Am,0)+2
        # quadric cones: (z_b-z_a)(z_c-z_a)=(z_d-z_e..) for base a and pairings
        Z=[pts[s][0] for s in S]; cone_hits=[]
        for a in range(5):
            others=[j for j in range(5) if j!=a]
            for (b,c) in [(others[0],others[1]),(others[0],others[2]),(others[0],others[3])]:
                d,e=[j for j in others if j not in (b,c)]
                lhs=pmul(psub(Z[b],Z[a]),psub(Z[c],Z[a])); rhs=pmul(psub(Z[d],Z[a]),psub(Z[e],Z[a]))
                if deg(psub(lhs,rhs))<0: cone_hits.append((S[a],S[b],S[c],S[d],S[e]))
        # rank of the 4 difference vectors
        import sympy as sp
        rk=sp.Matrix([[vecs[j][k]-vecs[0][k] for k in range(8)] for j in range(1,5)]).rank()
        res.append((S,g,mu0,BC,alpha,len(cuts),sorted(cuts.values()),rk,cone_hits))
    return res,gam,degs

if __name__=="__main__":
    y=[Fr(-9,4),Fr(-2,3),Fr(1,4),Fr(8,3)]
    res,gam,degs=analyse(y)
    from collections import Counter
    print("5-subsets:",len(res))
    print("gcd degrees:",Counter(r[1] for r in res)," mu0:",Counter(r[2] for r in res)," alpha:",Counter(str(r[4]) for r in res))
    print("#cuts used:",Counter(r[5] for r in res)," ranks:",Counter(r[7] for r in res))
    print("cone hits total:",sum(len(r[8]) for r in res))
    for r in res[:3]: print(r)
