# Exact: extra zero-sum sets on the cube locus (Thm C of round 29).
# y=(a,b,c,d), d=-(a+b+c), abcd=1, eps=(1,-1,-1,1), r_k=eps_k/y_k+i y_k; roots of Q are +-r_k.
# Set F <-> (a_k,b_k) in {0,1}^2 per k (r_k in F, -r_k in F); s=a-b, t=a+b.
# m=1: sum s_k y_k=0; m=2: sum t_k eps_k=0 (Im r^2=2 eps); m=3: sum s_k(3/y_k-y_k^3)=0.
import sympy as sp, itertools
a,b,c=sp.symbols('a b c')
eps=[1,-1,-1,1]
y=[a,b,c,-(a+b+c)]
results=[]
for s in itertools.product([-1,0,1],repeat=4):
    if all(v==0 for v in s): continue
    if s in [(1,1,1,1),(-1,-1,-1,-1)]: continue
    # need some t: t_k=1 if s_k!=0 else t_k in {0,2}, with sum t eps=0
    ok_t=False
    for tt in itertools.product(*[[1] if s[k] else [0,2] for k in range(4)]):
        if sum(tt[k]*eps[k] for k in range(4))==0: ok_t=True
    if not ok_t: continue
    e1=sp.expand(sum(s[k]*y[k] for k in range(4)))
    e3=sp.together(sum(s[k]*(3/y[k]-y[k]**3) for k in range(4)))
    n3=sp.numer(e3)
    g=sp.expand(a*b*c*(-(a+b+c))-1)
    results.append((s,e1,sp.factor(n3),g))
import pickle
for s,e1,n3,g in results: print(s, "| m1:",e1,"| m3 numer:",n3)
pickle.dump([(s,str(e1),str(n3)) for s,e1,n3,g in results],open("cube_locus_eqs.pkl","wb"))
