# Minimum support of nonzero reductions mod a prime above l of integral eigenvectors of the Paley
# skew conference matrix C (order q+1, C^2=-qI).  Exhaustive enumeration of the reduced code.
# l=2: y = (C-1)b/2 + b*omega, omega=(1+g)/2, b in Z^{q+1} with Cb = b mod 2 (derived in the report);
#      reduction: omega -> r, r in F_2 (q=7 mod 8, split) or r in F_4 root of x^2+x+1 (q=3 mod 8, inert).
# l odd split: code = ker(C - r) over F_l with r^2 = -q.
import numpy as np, itertools, math, sys
def C_of(q):
    sq=set((i*i)%q for i in range(1,q)); c=np.array([0]+[1 if i in sq else -1 for i in range(1,q)])
    idx=np.arange(q); K=c[(idx[None,:]-idx[:,None])%q]
    C=np.zeros((q+1,q+1),dtype=np.int64); C[:q,:q]=K; C[:q,q]=1; C[q,:q]=-1; return C
def rref_mod(M,p):
    M=M.copy()%p; r=0; piv=[]
    for c in range(M.shape[1]):
        k=next((i for i in range(r,M.shape[0]) if M[i,c]),None)
        if k is None: continue
        M[[r,k]]=M[[k,r]]; M[r]=M[r]*pow(int(M[r,c]),-1,p)%p
        for i in range(M.shape[0]):
            if i!=r and M[i,c]: M[i]=(M[i]-M[i,c]*M[r])%p
        piv.append(c); r+=1
        if r==M.shape[0]: break
    return M[:r],piv
def kernel_mod(A,p):
    R,piv=rref_mod(A,p); n=A.shape[1]; free=[c for c in range(n) if c not in piv]; basis=[]
    for f in free:
        v=np.zeros(n,dtype=np.int64); v[f]=1
        for i,c in enumerate(piv): v[c]=(-R[i,f])%p
        basis.append(v)
    return np.array(basis)
def span_basis(G,p):
    R,piv=rref_mod(G,p); return R
def minweight_Fp(G,p):
    k=G.shape[0]; best=10**9
    for coeffs in itertools.product(range(p),repeat=k):
        if not any(coeffs): continue
        w=int((np.array(coeffs)@G%p!=0).sum()); best=min(best,w)
    return best
def code_l2(q):
    C=C_of(q); n=q+1
    ker=kernel_mod((C-np.eye(n,dtype=np.int64))%2,2)   # b mod 2 with Cb=b mod 2
    gens=list(ker)+[2*np.eye(n,dtype=np.int64)[i] for i in range(n)]
    # y(b) = a + b*omega, a=(C-1)b/2 ; represent reduction over F_2[omega]/(minpoly mod 2)
    A=[];B=[]
    for b in gens:
        a=(C-np.eye(n,dtype=np.int64))@b; assert (a%2==0).all(); A.append((a//2)%2); B.append(b%2)
    return np.array(A),np.array(B)
def mw_l2(q):
    A,B=code_l2(q); n=q+1
    if q%8==7:
        out=[]
        for r in (0,1):
            G=span_basis((A+r*B)%2,2); out.append((G.shape[0],minweight_Fp(G,2)))
        return out
    else:  # F_4: element a+b*w, w^2=w+1 ; F_2-span of (A,B) pairs
        G=np.concatenate([A,B],axis=1)%2; G=span_basis(G,2); k=G.shape[0]; best=10**9
        for coeffs in itertools.product((0,1),repeat=k):
            if not any(coeffs): continue
            x=np.array(coeffs)@G%2; w=int(((x[:n]|x[n:])!=0).sum()); best=min(best,w)
        return [(k,best)]
def mw_odd_split(q,l):
    C=C_of(q); n=q+1
    rs=[r for r in range(l) if (r*r+q)%l==0]
    out=[]
    for r in rs:
        G=kernel_mod((C-r*np.eye(n,dtype=np.int64))%l,l); out.append((G.shape[0],minweight_Fp(G,l)))
    return out
for q in [7,11,19,23,31]:
    print("q=%d l=2 (%s): (dim,minweight) ="%(q,'split' if q%8==7 else 'inert F4'),mw_l2(q),
          " sqrt(q)=%.2f, need d^2-d+1>=q"%math.sqrt(q))
for q,l in [(11,3),(23,3),(11,5),(19,5),(19,7),(7,3)]:
    if any((r*r+q)%l==0 for r in range(l)):
        print("q=%d l=%d split: (dim,minweight) ="%(q,l),mw_odd_split(q,l)," sqrt(q)=%.2f"%math.sqrt(q))
    else: print("q=%d l=%d inert (skipped)"%(q,l))
