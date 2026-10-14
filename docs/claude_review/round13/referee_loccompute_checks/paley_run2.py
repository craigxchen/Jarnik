from paley_fam import *
V=np.load("paley_V.npy"); w=np.load("paley_w.npy")
acc=Acc(); acc.B=pickle.load(open(sys.argv[1],"rb")); acc.r=acc.B.shape[0]
print("start",acc.r)
def sgn(v):
    v=tuple(int(x) for x in v); f=next(x for x in v if x)
    return v if f>0 else tuple(-x for x in v)
C=[]
C+= [sgn(COLS[i]) for i in range(11)]
C+= [sgn(np.ones(11,dtype=np.int64))]
for j in range(11):
    e=np.zeros(11,dtype=np.int64); e[j]=1; C.append(sgn(e))
for a in range(11):
    for b in range(a+1,11):
        for s in (1,-1):
            e=np.zeros(11,dtype=np.int64); e[a]=1; e[b]=s; C.append(sgn(e))
for v in V[w==7]:
    if np.abs(v).sum()>3: C.append(sgn(v))
C=list(dict.fromkeys(C))
print("candidates",len(C))
def rot(v,t): return tuple(v[(j-t)%11] for j in range(11))
seen=set(); n=0; t0=time.time(); ng=0
for i in range(len(C)):
    for j in range(i+1,len(C)):
        key=min(tuple(sorted((sgn(rot(C[i],t)),sgn(rot(C[j],t))))) for t in range(11))
        if key in seen: continue
        seen.add(key); n+=1
        rows,sz=family_rows([np.array(C[i]),np.array(C[j])])
        if rows is None: continue
        ng+=1
        k=acc.add(rows)
        if k: print("pair",C[i],C[j],"|img|",sz,"new",k,"total",acc.r,flush=True)
print("pairs tried",n,"nongrid",ng,"final",acc.r,"time",time.time()-t0)
pickle.dump(acc.B,open(sys.argv[2],"wb"))
