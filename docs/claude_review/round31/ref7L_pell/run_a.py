import sys; sys.path.insert(0,'.')
from pell_core import *
A=[[P(1,0,1),P(0,1)],[P(0,1),P(1)]]
beta=P(G(0,2),1); betab=pconj(beta)
def family(n,seeds,edges,use_beta=True):
    An=matpow(A,n); Gs=[block(An,s) for s in seeds]; Gb=[pconj(x) for x in Gs]; r=len(seeds)
    z0=pmul(beta if use_beta else P(1),prod(Gb))
    marks=[z0]
    for (j,k) in edges:
        f=[]
        for l in range(r):
            e=(l==j)+(l==k)
            f+= [Gs[l]]*e + [Gb[l]]*(1-e) if e<=1 else [Gs[l]]*2
        # loop jj: G_j^2 and no conj(G_j) -> degree changes; handle generally
        marks.append(pmul(betab if use_beta else P(1),prod(f)))
    return marks,Gs
if __name__=="__main__":
    seedsets=[[Fr(1,2),Fr(1,3),Fr(2,7),Fr(-5,3)],[Fr(3),Fr(-1,5),Fr(7,4),Fr(2,9)],[Fr(1),Fr(2),Fr(3),Fr(5)]]
    print("five marks r=4 edges 12,13,23,14")
    for n in range(1,6):
        for S in seedsets:
            m,_=family(n,S,[(0,1),(0,2),(1,2),(0,3)])
            a=analyse(m); print(n,[str(x) for x in S],"D",a['D'],"8n+1",8*n+1,"h",a['h'],"4n",4*n,"gcd",a['gcddeg'],"eqnorm",a['eqnorm'],"distinct",a['distinct'],"pairc",sorted(set(a['pairc'].values())))
    print("four marks r=3 edges 12,13,23")
    for n in range(1,6):
        for S in seedsets:
            m,_=family(n,S[:3],[(0,1),(0,2),(1,2)])
            a=analyse(m); print(n,[str(x) for x in S[:3]],"D",a['D'],"6n+1",6*n+1,"h",a['h'],"gcd",a['gcddeg'],a['eqnorm'],a['distinct'])
    print("three marks G_j=t^d-a_j+i")
    for d in [1,2,3,4,6,9]:
        for avals in [[0,1,2],[Fr(1,2),-3,Fr(7,5)]]:
            Gs=[padd(padd([Z0]*d+[ONE],P(-G(a))),P(I)) for a in avals]; Gb=[pconj(x) for x in Gs]
            m=[pmul(Gs[j],prod([Gb[l] for l in range(3) if l!=j])) for j in range(3)]
            a=analyse(m); print(d,avals,"D",a['D'],"3d",3*d,"h",a['h'],"2d",2*d,"gcd",a['gcddeg'],a['eqnorm'],a['distinct'])
