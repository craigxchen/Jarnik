"""Exact five- and six-mark integer controls with all pair-cut divisibilities."""
from fractions import Fraction
from itertools import combinations, permutations
from math import gcd, prod

from iii_quartet import split_primes
from iii_quartet_exhaust import small_relation
from iii_horizontal_stars import root_of_minus_one


def balanced_flip(u,k,g,parent):
    # A polynomial row in the congruence lattice for a cycle of five pair cuts.
    aa=(3566328745136157,-212247537740,-232553177)
    bb=(-68417652775626,-13737560909,189877)
    def ev(f,w):
        return f[0]+f[1]*w+f[2]*w*w
    syl=(aa+(0,), (0,)+aa, bb+(0,), (0,)+bb)
    det=sum((-1)**sum(ix[i]>ix[j] for i,j in combinations(range(4),2))*prod(row[col] for row,col in zip(syl,ix)) for ix in permutations(range(4)))
    assert det!=0
    cycle=[(0,1),(1,2),(2,3),(3,4),(0,4)]
    r,d=0,1
    for i,j in cycle:
        w=-u[i]*u[j]
        assert ev(aa,w)+(u[i]+u[j])*ev(bb,w)==0
        mod=g[i,j]
        r+=d*((-u[i]-u[j]-r)*pow(d,-1,mod)%mod)
        d*=mod
    a,b=ev(aa,k*k),ev(bb,k*k)
    assert (a-r*b)%d==0
    alpha,beta=(a-r*b)//d,b
    content=gcd(alpha,beta)
    alpha//=content
    beta//=content
    v=(a//content,b//content)
    # alpha*delta-beta*gamma=1 in the lattice basis (d,0),(r,1).
    delta=pow(alpha,-1,abs(beta))
    gamma=(alpha*delta-1)//beta
    w=(d*gamma+r*delta,delta)
    vv=v[0]**2+v[1]**2
    n=(2*(w[0]*v[0]+w[1]*v[1])+vv)//(2*vv)
    w=(w[0]-n*v[0],w[1]-n*v[1])
    shift=(k**7-w[1])//v[1]
    old=[(s,s*s-k*k) for s in u]
    factors=[prod(g[e] for e in cycle if i in e) for i in range(5)]
    weights=[prod(g[e] for e in g if ((i in e)!=(e in cycle))) for i in range(5)]
    for j in range(3):
        row=(w[0]+(shift+j)*v[0],w[1]+(shift+j)*v[1])
        new=[]
        for (x,y),factor in zip(old,factors):
            xx,yy=row[0]*x+row[1]*y,v[0]*x+v[1]*y
            assert xx%factor==yy%factor==0
            new.append((xx//factor,yy//factor))
        if sum(t*t for t in new[0])!=weights[0]:
            break
    else:
        raise AssertionError('Three equal lengths from a nonconstant linear coordinate')
    norms=[x*x+y*y for x,y in new]
    for i,(x,y) in enumerate(new):
        assert gcd(x,y)==1 and 10**12*abs(norms[i]-weights[i])<weights[i]
        assert 10**12*abs(norms[i]**2-parent)<parent
    for i,j in combinations(range(5),2):
        x,y=new[i]
        a,b=new[j]
        shared=g[i,j]*d//(factors[i]*factors[j])
        assert abs(x*b-y*a)==shared*(u[j]-u[i]) and x*a+y*b>0
        assert 8**4*parent*(x*b-y*a)**4<(norms[i]*norms[j])**2
    print('Balanced five: five marked blocks per star, all ten full divisibilities,',
          'norm-squared and radial fourth-power relative errors below 1e-12.',flush=True)


def run(m):
    s=(3,7,11,19,23,31)[:m]
    u=[v*v for v in s]
    edges=list(combinations(range(m),2))
    c={e:u[e[0]]*u[e[1]] for e in edges}
    period=2*prod(abs(a-b) for a,b in combinations(c.values(),2))
    for p in s:
        while period%p==0:
            period//=p
    # K is a multiple of period and is one modulo the selected inert primes.
    k=0
    for p in s:
        k+=period*((1-k)*pow(period,-1,p)%p)
        period*=p
    forced=split_primes(max(c.values())+1,len(edges))
    for e,p in zip(edges,forced):
        residue=root_of_minus_one(p)*s[e[0]]*s[e[1]]%p
        if (residue**2+c[e])%(p*p)==0:
            residue+=p
        k+=period*((residue-k)*pow(period,-1,p*p)%(p*p))
        period*=p*p
    k+=period
    g={e:k*k+c[e] for e in edges}
    assert all(gcd(a,b)==1 for a,b in combinations(g.values(),2))
    assert all(g[e]%p==0 and g[e]%(p*p) for e,p in zip(edges,forced))
    parent=prod(g.values())
    base=[(v+k**(m-3)*(v*v-k*k),v*v-k*k) for v in u]
    norms=[x*x+y*y for x,y in base]
    nonisotropic=0
    for i,(x,y) in enumerate(base):
        assert gcd(x,y)==1
        intended=prod(g[e] for e in edges if i in e)
        assert 10**12*abs(norms[i]-intended)<intended
    for e,p in zip(edges,forced):
        i,j=e
        x,y=base[i]
        a,b=base[j]
        det=x*b-y*a
        assert det==g[e]*(u[j]-u[i]) and x*a+y*b>0
        assert 8**4*parent*det**4<(norms[i]*norms[j])**2
        nonisotropic+=bool(norms[i]%p)
    assert nonisotropic>0
    if m==5:
        balanced_flip(u,k,g,parent)
    cap=(max(u)-min(u))**2
    def block(i,j):
        return g[tuple(sorted((i,j)))]
    def inferred(i,j,k,l):
        a=block(i,j)*block(k,l)
        b=block(i,k)*block(j,l)
        c=block(i,l)*block(j,k)
        ray=small_relation(a,b,c,cap)
        assert ray is not None
        x,y,z=ray # convention a*x+b*y+c*z=0, middle product carries a minus
        return Fraction(-b*y,c*z)
    # Cache all quartet cross ratios needed for every choice of the anchor pair.
    ratios={v:inferred(*v) for v in permutations(range(m),4)}
    for i,j,k,l,h in permutations(range(m),5):
        assert ratios[i,j,k,l]*ratios[i,j,l,h]==ratios[i,j,k,h]
    print('marks, pair atoms, distinguished primes with nonisotropic lines, ordered five-marks:',
          m,len(edges),nonisotropic,len(list(permutations(range(m),5))),flush=True)


if __name__=='__main__':
    for marks in (5,6):
        run(marks)
