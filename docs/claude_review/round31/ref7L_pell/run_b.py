import sys; sys.path.insert(0,'.')
from pell_core import *
import random
def raw_contact(marks):
    D=deg(marks[0]); return D-max(deg(psub(a,b)) for a,b in itertools.combinations(marks,2))
def build(A,beta,n,seeds,edges,anchor_beta_pow=1):
    An=matpow(A,n); Gs=[block(An,s) for s in seeds]; Gb=[pconj(x) for x in Gs]; r=len(seeds)
    betab=pconj(beta)
    marks=[pmul(beta,prod(Gb))]
    for (j,k) in edges:
        f=[Gs[l] if l in (j,k) else Gb[l] for l in range(r)]
        marks.append(pmul(betab,prod(f)))
    return marks,Gs,Gb
A=[[P(1,0,1),P(0,1)],[P(0,1),P(1)]]; beta=P(G(0,2),1)
print("exact identities")
for n in [1,2,3,4]:
    An=matpow(A,n)
    for s in [Fr(1,2),Fr(-3),Fr(5,7)]:
        Gj=block(An,s); Gb=pconj(Gj)
        d=psub(pmul(pmul(Gj,Gj),pconj(beta)),pmul(pmul(Gb,Gb),beta))
        print(" n",n,"s",s,"deg(G^2 bb - Gb^2 b)=",deg(d),"lead",d[-1], "4is=",G(0,4*s), "degG",deg(Gj))
    s1,s2=Fr(1,2),Fr(-3)
    d=psub(pmul(block(An,s1),pconj(block(An,s2))),pmul(pconj(block(An,s1)),block(An,s2)))
    print(" G1 Gb2 - Gb1 G2 =",d, " expected 2i(s1-s2) up to sign:",G(0,2*(s1-s2)))
print("raw contacts (no leading-coefficient normalisation)")
for n in [1,2,3]:
    m,_,_=build(A,beta,n,[Fr(1,2),Fr(1,3),Fr(2,7),Fr(-5,3)],[(0,1),(0,2),(1,2),(0,3)])
    print(" n",n,"raw h",raw_contact(m),"lcs",[x[-1] for x in m])
print("special seeds s_j+s_k=0")
for n in [1,2,3]:
    for S in [[Fr(1),Fr(-1),Fr(2),Fr(3)],[Fr(1),Fr(-1),Fr(-1,1)+Fr(1,2),Fr(-1,2)],[Fr(1),Fr(-1),Fr(5),Fr(-1)]]:
        try:
            m,_,_=build(A,beta,n,S,[(0,1),(0,2),(1,2),(0,3)])
            a=analyse(m); print(" n",n,S,"D",a['D'],"h",a['h'],"gcd",a['gcddeg'],"distinct",a['distinct'],"pairc",a['pairc'])
        except AssertionError as e: print(" degree mismatch",S)
