import sys, itertools, numpy as np, sympy as sp, pickle
fn = sys.argv[1]; midgen = 'mid' in fn
GENT = [('q',(i,j,k,l)) for i,j,k,l in itertools.combinations(range(8),4)]
if midgen: GENT += [('m',(i,j,k)) for i,j,k in itertools.combinations(range(8),3)]
def gvec(t):
    v=[0]*8
    if t[0]=='q': i,j,k,l=t[1]; v[i]+=1; v[l]+=1; v[j]-=1; v[k]-=1
    else: i,j,k=t[1]; v[i]+=1; v[k]+=1; v[j]-=2
    return v
# all relations: pairings {a,b}|{c,d} and midpoints (j;{i,k})
REL=[]
for q in itertools.combinations(range(8),4):
    a,b,c,d=q
    for p in (((a,b),(c,d)),((a,c),(b,d)),((a,d),(b,c))): REL.append(('q',frozenset(map(frozenset,p))))
for j in range(8):
    for i,k in itertools.combinations([x for x in range(8) if x!=j],2): REL.append(('m',(j,frozenset((i,k)))))
RIDX={r:n for n,r in enumerate(REL)}
def rvec(r):
    v=[0]*8
    if r[0]=='q':
        p1,p2=list(r[1]);
        for x in p1: v[x]+=1
        for x in p2: v[x]-=1
    else:
        j,ik=r[1]; v[j]-=2
        for x in ik: v[x]+=1
    return v
RV=[sp.Matrix([rvec(r)]) for r in REL]
def apply(perm,r):
    if r[0]=='q': return ('q',frozenset(frozenset(perm[x] for x in p) for p in r[1]))
    j,ik=r[1]; return ('m',(perm[j],frozenset(perm[x] for x in ik)))
PP=10007
def reduce_mod(M,B):
    M=M.copy()%PP
    for piv,row in B: M=(M-np.outer(M[:,piv],row))%PP
    return M
def rref_mod(rows):
    B=[]
    for v in rows:
        v=reduce_mod(np.array([v],dtype=np.int64),B)[0]; nz=np.nonzero(v)[0]
        if len(nz)==0: continue
        p=nz[0]; v=v*pow(int(v[p]),PP-2,PP)%PP
        B=[(q,(r-r[p]*v)%PP) for q,r in B]+[(p,v)]
    return B
RVN=np.array([rvec(r) for r in REL],dtype=np.int64)
perms=list(itertools.permutations(range(8)))
PA=np.array([[RIDX[apply(p,r)] for r in REL] for p in perms],dtype=np.int16)
recs=[]
for line in open(fn):
    if not line.startswith('REC'): continue
    t=line.split(); rank=int(t[2]); ids=list(map(int,t[t.index('rows')+1:]))
    B=rref_mod([gvec(GENT[i]) for i in ids]); assert len(B)==rank
    Rm=reduce_mod(RVN,B); rel=[int(n) for n in np.nonzero(~Rm.any(1))[0]]
    nq=sum(1 for n in rel if REL[n][0]=='q'); nm=len(rel)-nq
    S=np.array(rel)
    imgs=np.sort(PA[:,S],axis=1)
    o=np.lexsort(imgs.T[::-1]); canon=tuple(imgs[o[0]].tolist())
    recs.append(dict(rank=rank,ids=ids,ngen=int(t[4]),nq=nq,nm=nm,canon=canon,rel=rel))
classes={}
for r in recs: classes.setdefault((r['rank'],r['canon']),[]).append(r)
print(len(recs),'records',len(classes),'classes')
out=[]
for (rank,c),L in sorted(classes.items(),key=lambda kv:(kv[0][0],-len(kv[0][1]))):
    r=L[0]; feas=[x for x in L if x['ngen']==x['nq']+(x['nm'] if midgen else 0)]
    print('rank',rank,'nq',r['nq'],'nmid',r['nm'],'records',len(L),'ordered-consistent records',len(feas),'rep',(feas or L)[0]['ids'])
    out.append(dict(rank=rank,nq=r['nq'],nm=r['nm'],nrec=len(L),nfeas=len(feas),rep=(feas or L)[0]['ids'],rel=[REL[n] for n in r['rel']]))
pickle.dump((GENT,out),open(fn+'.classes.pkl','wb'))
