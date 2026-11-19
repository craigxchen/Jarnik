import sys; sys.path.insert(0,'.')
from pell_core import *
A=[[P(1,0,1),P(0,1)],[P(0,1),P(1)]]; beta=P(G(0,2),1); bb=pconj(beta)
ok=True
for n in range(1,7):
    An=matpow(A,n)
    for s in [Fr(1,2),Fr(-3),Fr(0),Fr(7,3)]:
        for u in [Fr(2),Fr(-1,5),Fr(1,2)]:
            Gs,Gu=block(An,s),block(An,u); Gsb,Gub=pconj(Gs),pconj(Gu)
            bil=psub(pmul(pmul(Gs,Gu),bb),pmul(pmul(Gsb,Gub),beta))
            exp1=P(G(0,4)*G(-1+s*u),G(0,4)*G((s+u)/2))
            alt=psub(pmul(Gs,Gub),pmul(Gsb,Gu)); exp2=P(G(0,2)*G(s-u))
            ok&= (bil==exp1) and (alt==exp2)
print("bilinear invariants  G_sG_u bbar - Gbar_s Gbar_u beta = 4i((s+u)t/2 - 1 + su)  and  G_s Gbar_u - Gbar_s G_u = 2i(s-u):",ok)
