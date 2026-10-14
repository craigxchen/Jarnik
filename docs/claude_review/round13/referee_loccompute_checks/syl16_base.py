import numpy as np, itertools
pts=[p for p in itertools.product([0,1],repeat=4) if any(p)]
d=15
COLS=np.array([[1 if sum(a*x for a,x in zip(av,p))%2 else 0 for p in pts] for av in pts],dtype=np.int8)
def build(L=2):
    out=[]
    # enumerate first 8 coords x last 7 coords in blocks
    K1=np.array(list(itertools.product([-1,0,1],repeat=8)),dtype=np.int8)
    K2=np.array(list(itertools.product([-1,0,1],repeat=7)),dtype=np.int8)
    Y2=(K2.astype(np.int16)@COLS[:,8:].T.astype(np.int16))
    for k1 in K1:
        y1=COLS[:,:8].astype(np.int16)@k1.astype(np.int16)
        ok=np.all(np.abs(Y2+y1)<=L,axis=1)
        if ok.any():
            sel=K2[ok]
            out.append(np.hstack([np.repeat(k1[None,:],len(sel),0),sel]))
    return np.vstack(out)
if __name__=="__main__":
    A=build(2); print("|A|",len(A))
    np.save("syl16_A.npy",A)
    Y=A.astype(np.int16)@COLS.T.astype(np.int16)
    sig=A.sum(1); print("sigma range",sig.min(),sig.max())
    # a line: a,b,a+b
    ia=pts.index((0,0,0,1)); ib=pts.index((0,0,1,0)); ic=pts.index((0,0,1,1))
    print("line image size",len({tuple(r) for r in Y[:,[ia,ib,ic]].tolist()}))
    V=[pts.index(p) for p in pts if p[0]==0]  # a plane: 7 nonzero a with a_0=0
    print("plane cols",V,"plane image size",len({tuple(r) for r in Y[:,V].tolist()}))
