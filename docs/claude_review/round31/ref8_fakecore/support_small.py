# Minimum support of Kv mod 3,4,5,7 (and mod 2) over balanced primitive v with ||v||_1 in {2,4},
# exhaustive up to the affine symmetry x -> a x + b of F_q (which maps (Kv) to +-(Kv) permuted).
import numpy as np, sys, math
def K_of(q):
    sq=set((i*i)%q for i in range(1,q)); c=np.array([0]+[1 if i in sq else -1 for i in range(1,q)])
    idx=np.arange(q); return c[(idx[None,:]-idx[:,None])%q]
MODS=[2,3,4,5,7]
def supp(M,mod):  # M: (batch,q) integer
    return (M%mod!=0).sum(1)
def run(q):
    K=K_of(q); cols=K.T  # cols[j] = K[:,j]
    res={}
    # length 2: v=e0-e1
    kv=(cols[0]-cols[1])[None,:]
    res[2]={m:int(supp(kv,m).min()) for m in MODS}
    best={m:10**9 for m in MODS}; witness={}
    rest=list(range(2,q))
    # shape (+1,+1,-1,-1), case 0,1 positive
    for c in rest:
        ds=np.array([d for d in rest if d>c])
        if len(ds)==0: continue
        kv=(cols[0]+cols[1]-cols[c])[None,:]-cols[ds]
        for m in MODS:
            s=supp(kv,m); k=s.argmin()
            if s[k]<best[m]: best[m]=int(s[k]); witness[m]=('0+1-%d-%d'%(c,ds[k]))
    # case 0 positive, 1 negative, c positive, d negative (c != d)
    for c in rest:
        ds=np.array([d for d in rest if d!=c])
        kv=(cols[0]-cols[1]+cols[c])[None,:]-cols[ds]
        for m in MODS:
            s=supp(kv,m); k=s.argmin()
            if s[k]<best[m]: best[m]=int(s[k]); witness[m]=('0-1+%d-%d'%(c,ds[k]))
    # shape (2,-1,-1): 2e0-e1-e_c ; and (1,1,-2) is the negative
    ds=np.array(rest); kv=(2*cols[0]-cols[1])[None,:]-cols[ds]
    for m in MODS:
        s=supp(kv,m); k=s.argmin()
        if s[k]<best[m]: best[m]=int(s[k]); witness[m]=('2*0-1-%d'%ds[k])
    res[4]=best
    return res,witness
for q in [int(x) for x in sys.argv[1:]]:
    r,w=run(q)
    print("q=%d sqrt(q)=%.2f  bound h=2: %.2f  h=4: %.2f"%(q,math.sqrt(q),math.sqrt(q)-3,math.sqrt(q)-5))
    print("  h=2 min supp by modulus",r[2]); print("  h=4 min supp by modulus",r[4]); print("  witnesses",w)
