import pickle, itertools, mpmath as mp, numpy as np, collections
mp.mp.dps=60
ci,dl,pts,u,fx1,fy1=pickle.load(open('equi5_cert.pkl','rb'))
u=[mp.mpf(a) for a in u]
X=[mp.mpf(0),mp.mpf(fx1.numerator)/fx1.denominator]+u[:6]; Y=[mp.mpf(1),mp.mpf(fy1.numerator)/fy1.denominator]+u[6:]
rho=[mp.mpc(a,b) for a,b in zip(X,Y)]
D=8
mom=[[mp.im(r**m) for m in (1,2,3)] for r in rho]
allp=list(itertools.product([1,-1],repeat=D))
mv={s:tuple(sum(s[k]*mom[k][m] for k in range(D)) for m in range(3)) for s in allp}
groups=[]
for s in allp:
    for g in groups:
        if max(abs(mv[s][m]-mv[g[0]][m]) for m in range(3))<mp.mpf(10)**-45: g.append(s); break
    else: groups.append([s])
sz=collections.Counter(len(g) for g in groups); print('fibre size histogram',sorted(sz.items(),reverse=True))
reps=[mv[g[0]] for g in groups]
print('min gap between distinct fibres',mp.nstr(min(max(abs(a[m]-b[m]) for m in range(3)) for a,b in itertools.combinations(reps,2)),5))
for g in groups:
    if len(g)>=5:
        A=np.array(g); print(len(g),'rank',np.linalg.matrix_rank(A[1:]-A[0])); print(A)
big=[np.array(g) for g in groups if len(g)>=5]
base5=big[0]; V5=(base5[1:]-base5[0])//2
for g in big[1:]:
    Dg=(g[1:]-g[0])//2
    print(len(g),'rank of span(5-cluster diffs + these diffs):',np.linalg.matrix_rank(np.vstack([V5,Dg])))
# real kernel dim of the moment matrix
Mm=np.array([[float(mom[k][m]) for k in range(D)] for m in range(3)]); print('moment matrix rank',np.linalg.matrix_rank(Mm,1e-12))
