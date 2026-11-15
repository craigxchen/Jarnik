# Check (S3) on the cube cluster: Blaschke version prod (n_v/conj n_v)^{v_i} = 1+O(t^-(L+1)) vs the N_v version.
import sympy as sp
t=sp.symbols('t'); I=sp.I
y=[sp.Rational(-9,4),sp.Rational(-2,3),sp.Rational(1,4),sp.Rational(8,3)]; eps=[1,-1,-1,1]
r=[eps[k]/y[k]+I*y[k] for k in range(4)]
roots=[r[0],-r[0],r[1],-r[1],r[2],-r[2],r[3],-r[3]]   # base point Q roots
L=3
F={0,1,2,3}  # {+-r0,+-r1}
G={0,2,4,6}  # {r0,r1,r2,r3}
for S in (F,G):
    num=1;den=1;Nn=1;Nd=1
    for k in S:
        num*= (t-roots[k]); den*=(t-sp.conjugate(roots[k]))
    # Blaschke: G = prod (t-conj r)/(t-r) over S (flip)
    B=sp.series(sp.simplify(den/num).subs(t,1/sp.Symbol('s')),sp.Symbol('s'),0,6).removeO()
    print('Blaschke ratio series:',sp.nsimplify(sp.expand(B)))
# N_v version for the pair Q vs z_G: residue vector (single basis vector) v=1 on G (delta=+1 for roots chosen... ) 
# delta has signs relative to UHP: roots in G with Im>0 get +1 or -1
s=sp.Symbol('s')
pos=[k for k in G if sp.im(roots[k])>0]; neg=[k for k in G if sp.im(roots[k])<0]
Np=sp.prod([(t-roots[k])*(t-sp.conjugate(roots[k])) for k in pos]); Nm=sp.prod([(t-roots[k])*(t-sp.conjugate(roots[k])) for k in neg])
print('N+/N- series:',sp.expand(sp.series(sp.expand(Np).subs(t,1/s)/sp.expand(Nm).subs(t,1/s),s,0,5).removeO()))
