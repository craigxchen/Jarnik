# Exact checks of counterexamples to (S_min): min_a sum_{x!=a} h(u_x/u_a) <= (1/3-delta) n h(ell_u) + c(n)
import sys; sys.path.insert(0,'.')
from gauss import *
from math import log, sqrt, atan2, pi
def report(name,u,extra=''):
    assert nonreal_pairs(u)
    n=len(u)-1
    hl,G=h_line(u); s,a=minanchor(u); hu=h_proj(u)
    print(f"{name}: n={n} h(u)={hu:.4f} h(ell)={hl:.4f} G={G} minanchor={s:.4f}(a={a}) "
          f"ratio={s/hl:.4f} n/3={n/3:.4f} ratio/n={s/hl/n:.4f} excess S-(n/3)h(ell)={s-n/3*hl:.4f} {extra}")
print("== 1. rank-one (1:g:g^2), g=A+i ==")
for A in [10,10**3,10**6,10**12,10**24]:
    g=(A,1); report(f"A={A}",[(1,0),g,pw(g,2)], f"L=log(A^2+1)={log(A*A+1):.4f}")
print("== 2. rank-one (1:g:g^2:g^3) ==")
for A in [10,10**3,10**6,10**12,10**24]:
    g=(A,1); report(f"A={A}",[pw(g,k) for k in range(4)])
print("== 3. genuine short-arc 3-clusters z=((b-i)^2, b^2+1, (b+i)^2), b even ==")
for b in [10,100,1000,10**4,10**6,10**9]:
    z=[pw((b,-1),2),(b*b+1,0),pw((b,1),2)]
    N=norm(z[1]); assert all(norm(x)==N for x in z)
    g=ggcd(ggcd(z[0],z[1]),z[2]); assert norm(g)==1
    assert mul(z[0],z[2])==mul(z[1],z[1])  # multiplicative dependence z1 z3 = z2^2
    R=sqrt(N)  # = b^2+1
    arc=R*abs(atan2(z[2][1],z[2][0])-atan2(z[0][1],z[0][0]))
    report(f"b={b}",z,f"R={R:.4g} arc/sqrtR={arc/sqrt(R):.6f} log R={log(R):.4f}")
