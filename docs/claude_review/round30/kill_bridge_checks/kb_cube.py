# Independent kill check: cube cluster, exhaustive clique count, pair degrees, gamma, squarefreeness.
import itertools, sympy as sp, networkx as nx
from fractions import Fraction as Fr
t=sp.symbols('t')
I=sp.I
def cube_roots(y,eps=(1,-1,-1,1)):
    r=[sp.Rational(e)/v + I*v for e,v in zip(eps,y)]
    roots=[]
    for x in r: roots+=[x,-x]
    return roots
def build(y,eps=(1,-1,-1,1)):
    roots=cube_roots(y,eps)
    allr=roots+[sp.conjugate(x) for x in roots]
    assert len(set(sp.nsimplify(x) for x in allr))==16, "coinciding roots"
    assert all(sp.im(x)!=0 for x in roots)
    N=sp.expand(sp.prod([(t-x)*(t-sp.conjugate(x)) for x in roots]))
    Np=sp.Poly(N,t)
    sqf=sp.degree(sp.gcd(Np,Np.diff(t)))==0
    pts={}
    for mask in range(256):
        rr=[sp.conjugate(roots[k]) if (mask>>k)&1 else roots[k] for k in range(8)]
        pts[mask]=sp.Poly(sp.expand(sp.prod([t-x for x in rr])),t,domain='QQ_I')
    return roots,N,sqf,pts
if __name__=="__main__":
    y=[sp.Rational(-9,4),sp.Rational(-2,3),sp.Rational(1,4),sp.Rational(8,3)]
    assert sum(y)==0 and sp.prod(y)==1
    roots,N,sqf,pts=build(y)
    print("roots",roots); print("N squarefree:",sqf, "real coeffs:", all(sp.im(c)==0 for c in sp.Poly(N,t).all_coeffs()))
    G=nx.Graph(); G.add_nodes_from(range(256))
    for a,b in itertools.combinations(range(256),2):
        d=(pts[a]-pts[b])
        if d.is_zero or d.degree()<=4: G.add_edge(a,b)
    cl=list(nx.find_cliques(G)); m=max(len(c) for c in cl)
    best=[c for c in cl if len(c)==m]
    print("max clique size",m,"number of max cliques",len(best))
    from collections import Counter
    print("clique size distribution", Counter(len(c) for c in cl).most_common(6))
    c=sorted(best[0]); print("one max clique masks",[format(x,'08b')[::-1] for x in c])
    for a,b in itertools.combinations(c,2):
        d=pts[a]-pts[b]; co=d.all_coeffs()[::-1]; g=co[4] if len(co)>4 else 0
        print(format(a,'08b')[::-1],format(b,'08b')[::-1],"deg",d.degree(),"gamma",g, "|gamma|=",sp.N(sp.Abs(g),12))
