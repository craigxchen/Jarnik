"""Exact finite tests for Part X's all-height bounded-word residue model.
Run with python3 x_local_check.py. Requires SymPy and python-flint.
The H=2 case allows g=8: two same-sign tokens are every unbalanced word of mass 2.
"""
from math import gcd, isqrt
from itertools import combinations, combinations_with_replacement
from flint import fmpz_mat
from sympy import factorint, primerange
from sympy.ntheory.modular import crt
from ix_factor_check import ONE, I, mul, norm, power, conj, combined

G = 8


def factor(p):
    for a in range(1,isqrt(p)+1):
        b=isqrt(p-a*a)
        if a*a+b*b==p:
            return (a,b)
    raise ValueError(p)


def group(ell, a):
    Q=ell**a
    order=ell**(a-1)*(ell-1 if ell%4==1 else ell+1)
    for x in range(1,ell+3):
        for y in range(1,ell+3):
            if gcd(x*x+y*y,Q)!=1:
                continue
            g=mul((x,y),power((x,-y),-1,Q),Q)
            if all(power(g,order//l,Q)!=ONE for l in factorint(order)):
                break
        else:
            continue
        break
    else:
        raise ValueError((ell,a))
    log={}
    z=ONE
    for i in range(order):
        log[z]=i
        z=mul(z,g,Q)
    assert z==ONE and len(log)==order
    return g,order,log


def split(base_logs, total, order, tie=False):
    # For odd r >= 5 choose representatives 1..(r-1)/2; at 3 use 9,
    # and at 2 use 16 with prescribed parity. Thus every singleton and
    # pair sum has gcd dividing K=3*8 with the local order. 8 slots allow
    # the first two to be set equal wherever parities permit.
    components=[]
    moduli=[]
    for r,t in factorint(order).items():
        full=r**t
        threshold = 16 if r==2 else 9 if r==3 else r
        D=min(threshold,full)
        choices=[]
        for log in base_logs:
            if D==threshold:
                vals=range(1, (D-1)//2+1)
            else:
                vals=range(D)
            choices.append([z for z in vals if r!=2 or z%2==log%2])
        vals=[bag[0] for bag in choices]
        assert not tie or vals[0]==vals[1]
        active=2 if tie else 0
        # g=8 -> 6 free slots even with a tie, all in arithmetic progressions.
        step=2 if r==2 else 1
        excess=(total-sum(vals))%D
        assert excess%step==0
        for j in range(active,len(vals)):
            inc=min(excess,choices[j][-1]-vals[j])
            assert inc%step==0
            vals[j]+=inc
            excess-=inc
        assert excess==0 and sum(vals)%D==total%D
        vals[-1]=(vals[-1]+(total-sum(vals))%full)%full
        components.append(vals)
        moduli.append(full)
    new=[int(crt(moduli,[col[j] for col in components])[0]) for j in range(len(base_logs))]
    assert all((x-y)%2==0 for x,y in zip(new,base_logs))
    assert sum(new)%order==total%order and (not tie or new[0]==new[1])
    return new


def local_check():
    q=11
    labels=list(range(1,q))
    squares={x*x%q for x in labels}
    K=lambda i,j: 0 if i==j else 1 if (j-i)%q in squares else -1
    marks=[i for i in labels if K(i,0)==1]
    block=[j for j in labels for _ in range(G)]
    primes=[p for p in primerange(101,2000) if p%4==1][:len(block)]
    pi=[factor(p) for p in primes]
    A=[[1+K(i,j) for j in block] for i in marks]
    n=len(block)
    D=fmpz_mat([[y-x for x,y in zip(A[0],row)] for row in A[1:]])
    H=fmpz_mat([[int(i==j)-int(marks[0]==j) for j in block] for i in marks[1:]])
    _,trans=D.transpose().hnf(transform=True)
    U=fmpz_mat([[trans[i,j] for i in range(len(marks)-1)] for j in range(n)])
    assert (D*U).is_one()
    F=fmpz_mat([[int(j==k) for k in range(n)] for j in range(n)])+U*(H-D)
    assert D*F==H
    negative=next(ell for ell in primerange(1500,4500) if
                  pow(primes[0],(ell-1)//2,ell)==pow(primes[1],(ell-1)//2,ell))
    exterior=pairs=words=interior=towers=0
    max_gcd=0
    for ell,a in [(3,3),(5,3),(7,2),(11,2),(13,2),(17,2),(19,1),(41,1),(negative,1)]:
        Q=ell**a
        g,order,log=group(ell,a)
        gamma=[mul(z,power(conj(z),-1,Q),Q) for z in pi]
        physical=[log[u] for u in gamma]
        old=[sum(int(F[j,k])*physical[k] for k in range(n))%order for j in range(n)]
        new=[]
        for j in range(0,n,G):
            L=old[j:j+G]
            new.extend(split(L,sum(L)%order,order,tie=(ell==negative and j==0)))
        actual=[]
        for j,p in enumerate(primes):
            diff=(new[j]-physical[j])%order
            assert diff%2==0
            z=mul(pi[j],power(g,diff//2,Q),Q)
            assert norm(z)%Q==p%Q
            assert mul(z,power(conj(z),-1,Q),Q)==power(g,new[j],Q)
            actual.append(z)
        # Verify all exterior row pair identities as Gaussian residues.
        points=[]
        for row in A:
            mark=ONE
            for j,aa in enumerate(row):
                mark=mul(mark,mul(power(actual[j],aa,Q),power(conj(actual[j]),2-aa,Q),Q),Q)
            points.append(mark)
        owners={i:combined(gamma,[int(j==i) for j in block],Q) for i in marks}
        for u,v in combinations(range(len(marks)),2):
            assert mul(points[u],owners[marks[v]],Q)==mul(points[v],owners[marks[u]],Q)
            pairs+=1
        for b in range(1,a+1):
            low=ell**b
            ord_low=order//ell**(a-b)
            assert power(g,ord_low,low)==ONE
            for j in range(n):
                dd=gcd(new[j],ord_low)
                assert 24%dd==0
                max_gcd=max(max_gcd,dd)
                words+=1
            for j,k in combinations_with_replacement(range(n),2):
                dd=gcd(new[j]+new[k],ord_low)
                assert 24%dd==0
                max_gcd=max(max_gcd,dd)
                words+=1
            towers+=1
        if ell==negative:
            assert new[0]==new[1] and ell*ell>primes[0]*primes[1]
        exterior+=1
    # Same-level tests at all 80 physical primes by evaluating the old
    # factorization with the missing physical gamma omitted. The new design
    # keeps that local factorization, including at higher internal powers.
    for skip,p in enumerate(primes):
        for Q in (p,p*p):
            gamma=[mul(z,power(conj(z),-1,Q),Q) if j!=skip else ONE for j,z in enumerate(pi)]
            old=[combined(gamma,[int(F[j,k]) for k in range(n)],Q) if j!=skip else ONE
                 for j in range(n)]
            owners={i:combined(gamma,[int(j==i) for j in block],Q) for i in marks}
            for u,v in combinations(range(len(marks)),2):
                if A[u][skip]==A[v][skip]:
                    ratio=combined(old,[x-y for x,y in zip(A[u],A[v])],Q)
                    assert mul(ratio,owners[marks[v]],Q)==owners[marks[u]]
                    interior+=1
    print(f'{len(marks)} marks, {n} physical prime labels. Exterior designs {exterior}, reduced levels {towers}.')
    print(f'{pairs} exterior pair ratios and {interior} internal same-level ratios pass.')
    print(f'{words} singleton/pair-sum local orders pass, maximum index {max_gcd}, all divide 24.')
    print(f'Negative control: signed difference of labels {primes[0]},{primes[1]} hits +1 at {negative}, exceeding geometric-mean height.')


if __name__=='__main__':
    local_check()
