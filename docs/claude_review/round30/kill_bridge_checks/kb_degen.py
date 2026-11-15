# The family y = {a,-a,1/a,-1/a}: check every ordering against eps=(1,-1,-1,1) for coinciding roots (up to sign/conj).
import itertools, sympy as sp
a=sp.symbols('a',positive=True)
base=[a,-a,1/a,-1/a]; eps=(1,-1,-1,1); bad=0; tot=0
for o in set(itertools.permutations(range(4))):
    y=[base[i] for i in o]; tot+=1
    r=[e/v+sp.I*v for e,v in zip(eps,y)]
    allr=[]
    for x in r: allr+=[x,-x,sp.conjugate(x),-sp.conjugate(x)]
    coinc=any(sp.simplify(allr[i]-allr[j])==0 for i,j in itertools.combinations(range(16),2))
    bad+=coinc
print("orderings",tot,"degenerate for all a:",bad)
