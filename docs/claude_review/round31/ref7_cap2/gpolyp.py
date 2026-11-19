# Polynomials over F_9 = Z[i]/(3) as complex numpy arrays with integer parts in {0,1,2} (index = degree)
import numpy as np
P=3
MAXERR=[0.0]
def norm(a):
    re=np.mod(np.rint(a.real),P); im=np.mod(np.rint(a.imag),P)
    r=re+1j*im
    nz=np.nonzero(r)[0]
    return r[:nz[-1]+1] if len(nz) else np.zeros(1,complex)
def mul(a,b):
    la,lb=len(a),len(b); n=la+lb-1
    if min(la,lb)<=64:
        out=np.zeros(n,complex)
        small,big=(a,b) if la<=lb else (b,a)
        for i,c in enumerate(small):
            if c!=0: out[i:i+len(big)]+=c*big
        return norm(out)
    N=1<<(n-1).bit_length()
    r=np.fft.ifft(np.fft.fft(a,N)*np.fft.fft(b,N))[:n]
    err=max(np.abs(r.real-np.rint(r.real)).max(),np.abs(r.imag-np.rint(r.imag)).max())
    MAXERR[0]=max(MAXERR[0],err); assert err<0.2, err
    return norm(r)
def add(a,b):
    n=max(len(a),len(b)); out=np.zeros(n,complex); out[:len(a)]+=a; out[:len(b)]+=b; return norm(out)
def sub(a,b): return add(a,-b)
def deg(a):
    a=norm(a); return -1 if (len(a)==1 and a[0]==0) else len(a)-1
def conj(a): return norm(np.conj(a))
def prod(lst):
    lst=[norm(x) for x in lst]
    if not lst: return np.ones(1,complex)
    while len(lst)>1:
        lst.sort(key=len)
        nxt=[]
        for j in range(0,len(lst)-1,2): nxt.append(mul(lst[j],lst[j+1]))
        if len(lst)%2: nxt.append(lst[-1])
        lst=nxt
    return lst[0]
def powr(a,e):
    r=np.ones(1,complex)
    for _ in range(e): r=mul(r,a)
    return r
def poly(coeffs): return norm(np.array(coeffs,complex))
def tpow(d): z=np.zeros(d+1,complex); z[d]=1; return z
def deriv(a):
    a=norm(a)
    if len(a)==1: return np.zeros(1,complex)
    return norm(a[1:]*np.arange(1,len(a)))
def frob3(a):   # a(t)^3 for a in F_3[t] (real): substitute t^3, coefficients fixed
    out=np.zeros(3*(len(a)-1)+1,complex); out[::3]=a; return out
