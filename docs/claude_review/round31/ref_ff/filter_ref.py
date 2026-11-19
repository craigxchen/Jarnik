# Exact, complete Lemma-1.1 (saturated balance, L=3: min(a+,a-)>=2) test for the rank-5 six-row classes.
import pickle, itertools, sympy as sp, numpy as np
from fractions import Fraction as Fr
from scipy.optimize import linprog
reps=pickle.load(open('reps6.pkl','rb'))
def bits(r): return [(r>>j)&1 for j in range(8)]
def kernel(S):
    M=sp.Matrix([[a-b for a,b in zip(bits(r),bits(S[0]))] for r in S[1:]])
    return M.rank(), M, M.nullspace()
def find_witness(Hm):
    """Hm: 8x3 sympy (rows h_j). Return integer a!=0 with sum a_j h_j=0 and min(a+,a-)<=1, or None. Complete."""
    H=np.array(Hm.tolist(),dtype=float)
    # case one-signed: a>=0, sum a =1, H^T a=0
    res=linprog(np.zeros(8),A_eq=np.vstack([H.T,np.ones(8)]),b_eq=np.r_[0,0,0,1],bounds=[(0,None)]*8,method='highs')
    if res.status==0:
        return ('onesided',res.x)
    # find u with H u >= 1 (exists by Gordan); rationalise and verify exactly
    res=linprog(np.zeros(3),A_ub=-H,b_ub=-np.ones(8),bounds=[(None,None)]*3,method='highs')
    assert res.status==0
    u=[Fr(x).limit_denominator(1000) for x in res.x]
    y=[sum(Fr(int(sp.numer(Hm[j,i])),int(sp.denom(Hm[j,i])))*u[i] for i in range(3)) for j in range(8)]
    assert all(v>0 for v in y)
    hh=[[Fr(int(sp.numer(Hm[j,i])),int(sp.denom(Hm[j,i]))) for i in range(3)] for j in range(8)]
    # a = b - e_k, b>=0 integer, b_k=0, sum b_j h_j = h_k ; then sum b_j y_j = y_k bounds b
    for k in range(8):
        others=[j for j in range(8) if j!=k]
        def dfs(idx,b,rem,vec):
            if idx==len(others):
                if all(vec[i]==hh[k][i] for i in range(3)) and any(b): return list(b)
                return None
            j=others[idx]; m=0
            while m*y[j]<=rem:
                b.append(m)
                r=dfs(idx+1,b,rem-m*y[j],[vec[i]+m*hh[j][i] for i in range(3)])
                b.pop()
                if r: return r
                m+=1
            return None
        r=dfs(0,[],y[k],[Fr(0)]*3)
        if r:
            a=[0]*8
            for j,m in zip(others,r): a[j]=m
            a[k]-=1
            return ('mass1',a)
    return None
out=[]
nrank5=0
for S in reps:
    rk,M,ns=kernel(S)
    if rk!=5: continue
    nrank5+=1
    Hm=sp.Matrix.hstack(*ns)
    w=find_witness(Hm)
    if w is None: out.append(S)
    else:
        if w[0]=='mass1':
            a=sp.Matrix(w[1]); assert (Hm.T*a).is_zero_matrix  # exact check a in V
            ap=sum(x for x in w[1] if x>0); an=-sum(x for x in w[1] if x<0); assert min(ap,an)<=1
print('rank5',nrank5,'survivors',len(out))
pickle.dump(out,open('survivors.pkl','wb'))
# compare with judge's pickle
J=pickle.load(open('/home/user/Jarnik/docs/claude_review/round30/judge_checks/classes_D8_c0.pkl','rb'))
exec(open('canon_helper.py').read())
jk=set(canon_signs(v) for v in J.values())
mk=set(out)
print('judge',len(jk),'mine',len(mk),'equal',jk==mk)
