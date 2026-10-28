import sympy as sp, pickle
m,W,y=sp.symbols('m W y')
U,V=pickle.load(open('uv_branch1.pkl','rb'))
Q=sp.together(W**4+V*W**2-V**2*W+U*V**2)
num,den=sp.fraction(Q)
print('factor over Q(m):',sp.factor_list(sp.expand(num),W))
# substitute W = V*Z to simplify
Z=sp.symbols('Z')
Q2=sp.factor(sp.expand(sp.together((W**4+V*W**2-V**2*W+U*V**2).subs(W,V*Z)/V**2)))
print('Q(V Z)/V^2 =',Q2)
