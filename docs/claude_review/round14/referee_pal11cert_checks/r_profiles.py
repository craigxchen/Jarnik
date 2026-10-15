from collections import defaultdict
import itertools
att=defaultdict(set)
for ln in open('r_flats_out.txt'):
    a=ln.split(); att[int(a[2])].add(int(a[3]))
types=sorted((s,t) for s in att for t in att[s])
res=[]
def rec(start,k,ss,st,cur):
    if k==8:
        res.append(tuple(cur)); return
    for i in range(start,len(types)):
        s,t=types[i]
        # remaining 8-k-1 gens need sigma>=2, t>=1
        if ss+s+2*(7-k)>26 or st+t+1*(7-k)>44: continue
        cur.append(types[i]); rec(i,k+1,ss+s,st+t,cur); cur.pop()
rec(0,0,0,0,[])
print('profiles with sum sigma<=26, sum t<=44:',len(res))
expl={(2,6),(4,6),(4,9),(6,2),(6,4),(7,5),(12,1)}
for p in res:
    np_=sum(1 for x in p if x==(2,6)); par=[x for x in p if x not in expl]
    print(p, 'pairs',np_,'param',par, 'sum s',sum(x[0] for x in p),'sum t',sum(x[1] for x in p))
