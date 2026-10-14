from paley_fam import *
acc=Acc(); acc.B=pickle.load(open(sys.argv[1],"rb")); acc.r=acc.B.shape[0]
which=sys.argv[3]; b=int(sys.argv[4])
print("start",acc.r,"set",which,"b",b)
def sgn(v):
    v=tuple(int(x) for x in v); f=next(x for x in v if x)
    return v if f>0 else tuple(-x for x in v)
one=np.ones(11,dtype=np.int64)
def e(j):
    z=np.zeros(11,dtype=np.int64); z[j]=1; return z
C=[]
C+=[sgn(COLS[i]) for i in range(11)]; C+=[sgn(one)]
C+=[sgn(one-e(j)) for j in range(11)]
if 'c' in which: C+=[sgn(e(j)) for j in range(11)]
if 'k' in which: C+=[sgn(one-COLS[i]) for i in range(11)]
if 'p' in which:
    for a in range(11):
        for bb in range(a+1,11):
            C+=[sgn(e(a)+e(bb)), sgn(e(a)-e(bb))]
if 'q' in which:   # sigma - x_a - x_b  (= -(k_inf + k_a + k_b))
    for a in range(11):
        for bb in range(a+1,11): C+=[sgn(one-e(a)-e(bb))]
C=list(dict.fromkeys(C)); print("candidates",len(C))
def rot(v,t): return tuple(v[(j-t)%11] for j in range(11))
seen=set(); n=0; ng=0; t0=time.time()
for F in itertools.combinations(range(len(C)),b):
    key=min(tuple(sorted(sgn(rot(C[i],t)) for i in F)) for t in range(11))
    if key in seen: continue
    seen.add(key); n+=1
    U=[np.array(C[i]) for i in F]
    if np.linalg.matrix_rank(np.array(U))<b: continue
    rows,sz=family_rows(U)
    if rows is None: continue
    ng+=1
    k=acc.add(rows)
    if k: print("fam",[C[i] for i in F],"|img|",sz,"new",k,"total",acc.r,flush=True)
print("families",n,"nongrid",ng,"final",acc.r,"time",time.time()-t0)
pickle.dump(acc.B,open(sys.argv[2],"wb"))
