import ast, flint, re
def load(fn):
    s=open(fn).read().strip().rstrip(':')
    d=ast.literal_eval(s)
    if d==[-1] or (isinstance(d,list) and d[:1]==[-1]): return None
    if d[0]!=0: return ('inf',d)
    P=d[1]; p=P[0]; vars_=P[3]; lin=P[4]
    Q=P[5][1]
    w=flint.nmod_poly(Q[0][1],p)
    den=flint.nmod_poly(Q[1][1],p)
    nums=[flint.nmod_poly(x[0][1],p) for x in Q[2]]
    return dict(p=p,vars=vars_,lin=lin,w=w,den=den,nums=nums,deg=P[2])
def coords(par,sign=-1,usederiv=True):
    p=par['p']; w=par['w']
    dw=w.derivative() if usederiv else flint.nmod_poly([1],p)
    # inverse of dw mod w (assume squarefree w)
    g,a,b=w.xgcd(dw)  # a*w + b*dw = g
    inv=b* flint.nmod_poly([int(g.coeffs()[0])],p).__pow__(1) if False else None
    ginv=pow(int(g.coeffs()[0]),p-2,p)
    inv=(b*ginv)%w
    # variables: those with lin coefficient nonzero form theta
    vs=par['vars']; out={}
    others=[v for v,c in zip(vs,par['lin'])]  # nums correspond to all vars except the last? determine
    return inv
