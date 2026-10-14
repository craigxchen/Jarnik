from paley_fam import *
V=np.load("paley_V.npy"); w=np.load("paley_w.npy")
acc=Acc()
# width ideal J_w(R=1): all {-1,0,1} dirs with width<=7, orbit reps only (P0 invariance)
def canon(v):
    best=None
    for t in range(11):
        for s in (1,-1):
            c=tuple(int(x) for x in s*np.roll(v,t))
            if best is None or c>best: best=c
    return best
reps={}
for v in V[w<=7]:
    reps.setdefault(canon(v),v)
print("width<=7 orbit reps",len(reps))
for v in reps.values():
    acc.add(width_rows(v))
print("J_w(R=1) P0 dim",acc.r)
pickle.dump(acc.B,open("paley_Jw_B.pkl","wb"))
# columns 5-families (all), as author: should reach 82
cols=[COLS[i] for i in range(11)]
seen=set(); t0=time.time()
for F in itertools.combinations(range(11),5):
    c=min(tuple(sorted((i+t)%11 for i in F)) for t in range(11))
    if c in seen: continue
    seen.add(c)
    rows,sz=family_rows([cols[i] for i in c])
    if rows is None: continue
    n=acc.add(rows)
    print("5-col",c,"|img|",sz,"new",n,"total",acc.r,flush=True)
print("L_col(5) P0 dim",acc.r, time.time()-t0)
pickle.dump(acc.B,open("paley_Lcol5_B.pkl","wb"))
