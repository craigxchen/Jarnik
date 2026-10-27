# Library for H_light on sdc_q. Exact integer arithmetic (numpy int64).
import numpy as np
def chi_vec(q):
    c=np.zeros(q,dtype=np.int64)
    for a in range(1,q):
        c[a]=1 if pow(a,(q-1)//2,q)==1 else -1
    return c
def patterns(q):
    # independent of sdc_sym.py but same definition: P_a={a} u {x: chi(a-x)=1}, P'_a={a} u {x: chi(a-x)=-1}
    ch=chi_vec(q); P=[];Pp=[]
    for a in range(q):
        P.append([a]+[x for x in range(q) if ch[(a-x)%q]==1])
        Pp.append([a]+[x for x in range(q) if ch[(a-x)%q]==-1])
    return P,Pp
def sigma_direct(lam,q,P=None,Pp=None):
    # lam: length q+1, index q = inf (in no pattern)
    if P is None: P,Pp=patterns(q)
    tot=sum(abs(int(lam[S].sum())) for S in P)+sum(abs(int(lam[S].sum())) for S in Pp)
    return tot-q  # sigma = sum(|lam(S)|-1/2) over 2q patterns
def conv_chi(f,q,ch=None):
    # c(a) = sum_t f(t) chi(a-t), exact via sparse loop
    if ch is None: ch=chi_vec(q)
    c=np.zeros(q,dtype=np.int64)
    idx=np.arange(q)
    for t in np.nonzero(f)[0]:
        c+=f[t]*ch[(idx-t)%q]
    return c
def N_formula(f,q,ch=None):
    F=int(f.sum()); h=F+f; c=conv_chi(f,q,ch)
    return int(np.maximum(np.abs(h),np.abs(c)).sum())
def sigma_formula(f,q,ch=None):
    return N_formula(f,q,ch)-q
