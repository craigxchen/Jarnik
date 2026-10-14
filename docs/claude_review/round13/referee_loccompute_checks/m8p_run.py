import sys, math
L=int(sys.argv[1]); Dm=int(sys.argv[2])
sys.argv=[sys.argv[0],str(L),str(Dm)]
exec(open('m8.py').read())
V=np.array(list(itertools.product(range(-2,3),repeat=7)),dtype=np.int64)
W=[]
for v in V:
    if not v.any() or math.gcd(*[int(abs(x)) for x in v])!=1 or next(x for x in v if x)<0: continue
    s=len(np.unique(A@v))
    if s<=Dm: W.append((E.pullback([1],[(s,)],[v]),s))
print("width gens",len(W))
for Dt in range(Dm+1):
    if JA[Dt]==0: continue
    acc=Acc(Dt); acc.add(rows_deg(W,Dt)); print("D",Dt,"J_w(R=2)",acc.r,"J(A) ub",JA[Dt],"EQUAL" if acc.r==JA[Dt] else "DIFF",flush=True)
