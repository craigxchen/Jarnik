from fractions import Fraction as Fr
from itertools import combinations_with_replacement
def esym(ks):
    e=[1]+[0]*6
    for k in ks:
        for j in range(6,0,-1): e[j]+=e[j-1]*k
    return e
cnt=0; ex=[]
for K in range(1,31):
    for rest in combinations_with_replacement(range(K+1),4):
        ks=(0,)+rest+(K,); e=esym(ks)
        a=Fr(e[1],2); b=Fr(e[2],2)-Fr(e[1]**2,8); c=Fr(e[3],2)-Fr(e[1]*e[2],4)+Fr(e[1]**3,16)
        r2=e[4]-(b*b+2*a*c)
        if r2==0:
            cnt+=1; r1=e[5]-2*b*c
            if len(ex)<8: ex.append((ks,r1,c))
print(cnt, ex)
