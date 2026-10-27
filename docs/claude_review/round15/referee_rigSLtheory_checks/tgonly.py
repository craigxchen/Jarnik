import sys,json,rlib,onemod
P=json.loads(sys.argv[1]); K=int(sys.argv[2]); F=json.load(open(sys.argv[3]))['F']
bad=[]
for q0 in F:
    pis,V,data=onemod.setup(P,K,F,q0)
    P0=set(v for v in V if data[v][2][1]==abs(data[v][0][1]) and data[v][2][3]==abs(data[v][0][3]))
    T=set()
    for v in P0:
        for w in P0:
            u=tuple(a+b for a,b in zip(v,w))
            if any(u) and u in P0: T.update([v,w,u])
    idx=onemod.hnf_index(list(T),len(P))
    if idx!=1: bad.append((q0,idx,len(P0)))
print(sys.argv[1],K,'n_forced',len(F),'index!=1 at:',bad)
