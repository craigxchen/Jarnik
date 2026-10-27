# For a given set,K,q0: take SL-at-phi0 labels outside every eps-class, lift to level 2 w/o level-2 collisions,
# and re-verify (1.2) exactly at phi0 for all v in V_K, all units, with all-level demand.
import sys, json, random, rlib, onemod
P=json.loads(sys.argv[1]); K=int(sys.argv[2]); F=json.load(open(sys.argv[3]))['F']; q0=int(sys.argv[4]); nw=int(sys.argv[5])
res=onemod.run(P,K,F,q0,cap=3000,detail=False)
labs=res['SLlabels'][:nw]
pis,V,data=onemod.setup(P,K,F,q0); r=len(P)
m,L=rlib.labels(pis,q0); n=m//4
m2=q0*m
rng=random.Random(5)
ok_all=True
for c in labs:
    # level-1 direct check
    ok=True; ncol=0
    for v in V:
        D,N,g=data[v]; s=sum(a*b for a,b in zip(v,c))%m
        for u in range(4):
            mu=g[u]*(q0 if s==u*n else 1)
            if mu>abs(D[u]): ok=False
        if s%n==0: ncol+=1
    # class check (independent recomputation)
    from math import gcd
    incl=any(all((c[j]-uu*e[j]*L[j])%n==0 for j in range(r)) for uu in range(1,n+1) if gcd(uu,n)==1 for e in __import__('itertools').product((1,-1),repeat=r))
    # lift: c'' = c~ + m t, with c~ = c or -c (convention), test both conventions
    lifted=None
    for sign in (1,-1):
        ct=[(sign*x)%m for x in c]
        for trial in range(20000):
            t=[rng.randrange(q0) for _ in range(r)]
            cc=[ct[j]+m*t[j] for j in range(r)]
            if all(sum(a*b for a,b in zip(v,cc))%(m2//4)!=0 for v in V):
                lifted=(sign,cc); break
        if lifted is None: break
    print(json.dumps({'q0':q0,'label':c,'SLphi0_exact':ok,'collisions':ncol,'in_eps_class':incl,'lift_found_both_conventions':lifted is not None,'lift':lifted}))
    ok_all &= ok and (not incl) and lifted is not None
print('ALL_OK',ok_all,'n_labels',len(labs),'SL counts',res['SL0'])
