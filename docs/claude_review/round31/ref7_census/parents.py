import sys
from witness import pmul,pconj
def dec(c):
    p=[]
    while c: c,d=divmod(c,9); p.append((d%3,d//3))
    return p
S=set(); E=0
for line in open(sys.argv[1]):
    _,n,*bl=line.split(); N=[(1,0)]
    for b in bl: x=dec(int(b)); N=pmul(N,pmul(x,pconj(x)))
    t=tuple(c[0] for c in N); S.add(t)
    if all(t[j]==0 for j in range(1,len(t),2)): E+=1
print(sys.argv[1],"distinct parents",len(S), "even parents", sum(1 for t in S if all(t[j]==0 for j in range(1,len(t),2))))
