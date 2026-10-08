"""Exact audits across every layer selection of small two-prime profiles."""
from itertools import product
from math import gcd


def add(z,w): return z[0]+w[0], z[1]+w[1]
def sub(z,w): return z[0]-w[0], z[1]-w[1]
def mul(z,w): return z[0]*w[0]-z[1]*w[1], z[0]*w[1]+z[1]*w[0]
def conj(z): return z[0],-z[1]
def scale(c,z): return c*z[0],c*z[1]
def norm(z): return z[0]**2+z[1]**2


def div(z,w):
    a,b=mul(z,conj(w)); d=norm(w)
    assert a%d==b%d==0
    return a//d,b//d


def power(z,n):
    out=(1,0)
    for _ in range(n): out=mul(out,z)
    return out


def valuation(z,p):
    out=0
    while z!=(0,0):
        a,b=mul(z,conj(p)); d=norm(p)
        if a%d or b%d: break
        z=(a//d,b//d); out+=1
    return out


def ggcd(z,w):
    while w!=(0,0):
        a,b=mul(z,conj(w)); d=norm(w)
        q=((2*a+d)//(2*d),(2*b+d)//(2*d))
        z,w=w,sub(z,mul(q,w))
    return z


def halfangle(z,anchor):
    if z==anchor: return 0,1,scale(2,anchor)
    if z==scale(-1,anchor): return 1,0,mul((0,1),scale(2,anchor))
    plus=add(z,anchor); minus=sub(z,anchor)
    numerator=mul(minus,conj(plus))
    assert numerator[0]==0
    g,h=numerator[1],norm(plus)
    c=gcd(abs(g),h); g//=c; h//=c
    lam=div(plus,(h,0))
    assert minus==mul((0,g),lam)
    return g,h,lam


primes=[(2,1),(3,2)]
units=[(1,0),(0,1),(-1,0),(0,-1)]
profile_count=row_count=0
for exponents in product(range(1,4),repeat=2):
    levels=list(product(*(range(e+1) for e in exponents)))
    N=5**exponents[0]*13**exponents[1]
    rows=[]
    for j,level in enumerate(levels):
        z=units[j%4]
        for p,e,t in zip(primes,exponents,level):
            z=mul(z,mul(power(p,t),power(conj(p),e-t)))
        rows.append(z)
    for flags in product((0,1),repeat=sum(exponents)):
        selected=[]; offset=0
        for e in exponents:
            selected.append(flags[offset:offset+e]); offset+=e
        fs=[sum(s) for s in selected]
        Q=5**fs[0]*13**fs[1]
        betas=[]
        for level in levels:
            beta=(1,0)
            for p,s,f,t in zip(primes,selected,fs,level):
                bt=sum(s[:t])
                beta=mul(beta,mul(power(p,bt),power(conj(p),f-bt)))
            betas.append(beta)
        for anchor_index,anchor in enumerate(rows):
            anchor_levels=levels[anchor_index]
            dQ=ggcd((Q,0),anchor)
            expected_E=1
            for p,e,f,t0 in zip(primes,exponents,fs,anchor_levels):
                expected_E*=norm(p)**min(f,e-f,t0,e-t0)
            assert norm(dQ)==Q*expected_E
            global_kappa=(0,0); min_ramified=100
            for z,beta,level in zip(rows,betas,levels):
                w=div(z,beta)
                g,h,lam=halfangle(z,anchor)
                gamma=(h,-g)
                kappa=mul(lam,conj(beta))
                assert mul(gamma,kappa)==mul(scale(2,anchor),conj(beta))
                assert mul(conj(gamma),kappa)==scale(2*Q,w)
                assert norm(kappa)*norm(gamma)==4*Q*N
                reduced=div(kappa,dQ)
                assert norm(reduced)*expected_E==norm(lam)
                global_kappa=ggcd(global_kappa,kappa)
                min_ramified=min(min_ramified,valuation(lam,(1,1)))
                for p,e,s,f,t0,t in zip(primes,exponents,selected,fs,anchor_levels,level):
                    bt=sum(s[:t])
                    assert valuation(kappa,p)==min(t0,t)+f-bt
                    assert valuation(kappa,conj(p))==e-max(t0,t)+bt
                if g:
                    c=gcd(abs(lam[0]),abs(lam[1])); sign=1 if g>0 else -1
                    v=div(mul((0,sign),lam),(c,0))
                    d=abs(g)*c; t=-sign*h*c
                    assert gcd(abs(v[0]),abs(v[1]))==1
                    assert sub(z,anchor)==scale(d,v)
                    assert scale(2,anchor)==mul(v,(-d,t))
                row_count+=1
            assert min_ramified in (1,2)
            assert norm(global_kappa)==norm(dQ)*2**min_ramified
            profile_count+=1
print(f'Checked {profile_count} anchored selected-layer profiles and {row_count} exact rows.')
