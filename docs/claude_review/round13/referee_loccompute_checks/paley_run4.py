from paley_fam import *
acc=Acc(); acc.B=pickle.load(open(sys.argv[1],"rb")); acc.r=acc.B.shape[0]
c=int(sys.argv[3]); extra=sys.argv[4]
print("start",acc.r,"sigma + %d coords"%c, "extra",extra)
one=np.ones(11,dtype=np.int64)
def e(j):
    z=np.zeros(11,dtype=np.int64); z[j]=1; return z
seen=set(); t0=time.time()
for S in itertools.combinations(range(11),c):
    key=min(tuple(sorted((i+t)%11 for i in S)) for t in range(11))
    if key in seen: continue
    seen.add(key)
    fams=[[one]+[e(i) for i in key]]
    if extra=='col':
        fams=[[one]+[e(i) for i in key]+[COLS[p]] for p in range(11)]
    for U in fams:
        if np.linalg.matrix_rank(np.array(U))<len(U): continue
        rows,sz=family_rows(U)
        if rows is None: continue
        k=acc.add(rows)
        if k: print("fam sigma+coords",key,"extra" if len(U)>c+1 else "","|img|",sz,"new",k,"total",acc.r,flush=True)
print("final",acc.r,"time",time.time()-t0)
pickle.dump(acc.B,open(sys.argv[2],"wb"))
