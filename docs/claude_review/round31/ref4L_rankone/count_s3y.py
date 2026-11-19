import itertools, sys
R=int(sys.argv[1]); n=0; ex=[]
def esym(v):
    e=[1,0,0,0,0,0]
    for a in v:
        for k in range(5,0,-1): e[k]+=e[k-1]*a
    return e
for c in itertools.combinations(range(-R,R+1),4):
    y5=-sum(c)
    if y5 in c or 0 in c or y5==0 or y5<c[3]: continue
    e=esym(c+(y5,))
    if e[2]*e[5]+e[3]*e[4]==0:
        n+=1; ex.append((c+(y5,),-e[4]))
print(n); print(ex[:20])
