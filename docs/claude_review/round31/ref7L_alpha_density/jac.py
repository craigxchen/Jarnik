import sympy as sp
from sympy import I, Rational as Q, Matrix, expand, simplify, cancel
from fam import *
pts = [((Q(1,2),Q(-2),Q(7,3),Q(4)),Q(1)), ((Q(-1,3),Q(5,2),Q(6),Q(-7)),Q(3)), ((Q(0),Q(1),Q(3,5),Q(-11,4)),Q(-2,7)), ((Q(1),Q(2),Q(3),Q(5)),Q(10))]
for n in (1,2):
    G,Gb,zs = marks(n, S)
    # identity: d*conj(c) - c*conj(d) where G = c + s d
    for g in G[:1]:
        c = g.subs(S[0],0); d = expand(g - c).coeff(S[0])
        print('n',n,'d conj c - c conj d =', expand(d*sp.conjugate(c)-c*sp.conjugate(d)))
    q = [zs[i]/zs[0] for i in range(1,5)]
    vars_ = list(S)+[t]
    # log-Jacobian: (dq/dv)/q
    J = Matrix(4,5, lambda i,j: sp.diff(q[i],vars_[j])/q[i])
    for sd,tv in pts:
        sub = dict(zip(S,sd)); sub[t]=tv
        M = J.subs(sub).applyfunc(lambda e: sp.nsimplify(sp.simplify(e)))
        Ms = M[:, :4]
        det = sp.simplify(Ms.det())
        Gv = [g.subs(sub) for g in G]
        pred = sp.simplify(-2*sp.prod([2*I/expand(gv*sp.conjugate(gv)) for gv in Gv]))
        print(' pt',sd,tv,'rank(4x5)=',M.rank(),'det seed block=',det,' predicted -2*prod(2i/|G|^2)=',pred, 'equal', sp.simplify(det-pred)==0)
