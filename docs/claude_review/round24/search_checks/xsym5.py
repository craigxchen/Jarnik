import sympy as sp, pickle, numpy as np
import xstruct as X
b1,b2=sp.symbols('b1 b2')
disc,cub,den=pickle.load(open('xsym4.pkl','rb'))
fl=sp.factor_list(disc)[1]
P=np.load('xfam_points.npy')
pts=[]
for p in P[::40]:
    rho,r4,sig,beta,q=X.unpack(p); pts.append((rho[0].imag,rho[1].imag))
for i,(f,e) in enumerate(fl):
    F=sp.lambdify((b1,b2),f,'numpy')
    vals=[abs(F(u,v)) for u,v in pts]
    deg=sp.Poly(f,b1,b2).total_degree()
    print(i,'deg',deg,'mult',e,'max |f| on family:',f"{max(vals):.2e}", ' sample f:', str(f)[:80] if deg<=4 else '')
