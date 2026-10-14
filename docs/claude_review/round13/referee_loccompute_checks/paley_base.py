import itertools, numpy as np
D=[0,1,3,4,5,9]
d=11
COLS=np.zeros((11,11),dtype=np.int64)
for i in range(11):
    for x in D: COLS[i,(i+x)%11]=1
def build(L=2):
    K=np.array(list(itertools.product([-1,0,1],repeat=11)),dtype=np.int64)
    Y=K@COLS.T
    ok=np.all(np.abs(Y)<=L,axis=1)
    return K[ok]
if __name__=="__main__":
    A=build(2)
    print("|A|",len(A))
    s=A.sum(1); print("sigma range",s.min(),s.max())
    # widths of all {-1,0,1} primitive directions
    V=np.array(list(itertools.product([-1,0,1],repeat=11)),dtype=np.int64)
    V=V[[next((x for x in v if x),0)>0 for v in V]]
    print("#dirs",len(V))
    widths=np.zeros(len(V),dtype=np.int64)
    for st in range(0,len(V),4000):
        Z=A@V[st:st+4000].T
        widths[st:st+4000]=Z.max(0)-Z.min(0)+1
    np.save("paley_V.npy",V); np.save("paley_w.npy",widths)
    for w in range(1,10):
        sel=V[widths==w]
        supp=np.abs(sel).sum(1)
        print("width",w,"count",len(sel),"support hist",np.bincount(supp,minlength=12).tolist())
