# Class of the Thm 4.1 (gauss.md) k=4 curves on Bl_e LM_4: toric intersections, multiplicity at e, B.C,
# Vandermonde-null check. Exact over Q(i) (sympy).
import sympy as sp, itertools
s=sp.symbols('s'); I=sp.I
def fourth(r3):
    # given three roots, the fourth making e1,e2 real (gauss.md §4.1)
    s1=sum(r3); s2=r3[0]*r3[1]+r3[0]*r3[2]+r3[1]*r3[2]
    a,b=sp.re(s1),sp.im(s1); 
    return sp.nsimplify((a*b-sp.im(s2))/b) - I*b
def family(X,m):
    xa=sp.Rational(X); xb=xa*(1-xa)/(m**2+xa); xc=1-xa-xb; u=m*xb
    r={(0,1):u+I*(xa-1),(0,2):u+I*(xb-1),(1,2):u+I*(xc-1),(0,1,2):I}
    r[(0,)]=fourth([r[(0,1)],r[(0,2)],r[(0,1,2)]]); r[(1,)]=fourth([r[(0,1)],r[(1,2)],r[(0,1,2)]]); r[(2,)]=fourth([r[(0,2)],r[(1,2)],r[(0,1,2)]])
    return {T:sp.nsimplify(sp.expand(v)) for T,v in r.items()}
def analyse(r,label):
    Ts=list(r)
    P=[sp.expand(sp.prod([s-r[T] for T in Ts if x in T])) for x in range(3)]
    imP=[sp.simplify(sp.im(sp.expand(p.subs(s,sp.Symbol('t',real=True))))) for p in P]
    consts=[sp.simplify(q) for q in imP]
    # toric intersection numbers over Q(i): at s=r_T the valuation vector of (u_0=1,u_x=P_x/conj P_x) is 1_T (in coords 0..3, index0 = point z_0)
    hits={}
    for T in Ts:
        v=tuple([0]+[1 if x in T else 0 for x in range(3)]); hits[v]=hits.get(v,0)+1
        vc=tuple([0]+[-1 if x in T else 0 for x in range(3)]); hits[vc]=hits.get(vc,0)+1
    # normalise vectors mod 1 to 0/1 form: -1_T = 1_{T^c}
    def norm(v): mn=min(v); return tuple(a-mn for a in v)
    rays={}
    for v,c in hits.items(): rays[norm(v)]=rays.get(norm(v),0)+c
    BC=sum(rays.values())
    # multiplicity at e: u_x - 1 = 2i c_x / conj(P_x): order deg P_x = 4 in t=1/s for each x with c_x != 0
    mult=min(sp.degree(p,s) for p in P)
    # Vandermonde (k=4) divisor: a_A + a_{A^c} = |A||A^c|; D_V.C = sum over hit rays a_A ; ord_e V = 6
    DV=sum(c*0 for c in rays.values())
    DV=sum(rays[v]*sum(1 for a in range(4) for b in range(4) if v[a]==1 and v[b]==0) for v in rays)  # a_A = #(a in A, b notin A)? 
    print(label,"Im P_x =",consts," #rays hit:",len(rays)," all mult 1:",set(rays.values())=={1}," B.C =",BC," mult_e =",mult,
          " B.C/mult =",sp.Rational(BC,mult))
    return rays
r0={(0,1):-6-4*I,(0,2):-6-3*I,(1,2):-6+5*I,(0,):-7+6*I,(1,):-3-2*I,(2,):-4-3*I,(0,1,2):I}
analyse(r0,"instance X=-3,m=3:")
for X,m in [(2,5),(sp.Rational(1,2),sp.Rational(3,2)),(-5,1),(7,2)]:
    analyse(family(X,m),f"family X={X},m={m}:")
# Vandermonde-null: with a_A for V normalised so that a_A + a_{A^c} = |A||A^c|, sum over all 14 rays (each hit once)
tot=sum(len(A)*(4-len(A)) for r in range(1,4) for A in itertools.combinations(range(4),r))
print("D_V.C = sum_A a_A = sum_cuts |A||A^c| =",tot//2*1, "(each ray once; sum over 14 rays of a_A = sum over 7 cuts of width) ; ord_e V * mult =",6*4)
