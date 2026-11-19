"""Independent check of the covering inequality (VII.2) mechanics on Astra's five-mark Pell family.
Log-derivative determinant Delta(t) of (d/dt, d/ds_a, d/ds_b, d/ds_c) on the ratio map.
Checks: Delta != 0, ord_inf Delta >= 4h+1, poles simple, located at charged places, #charged places <= W=2D."""
import sympy as sp, sys, itertools
t=sp.Symbol('t'); I=sp.I
def mats(n):
    A=sp.Matrix([[t**2+1,t],[t,1]]); An=A**n
    return sp.expand(An[0,0]),sp.expand(An[1,0]),sp.expand(An[0,1]),sp.expand(An[1,1])
def run(n,seeds,choose):
    a,c,b,d=mats(n)  # (X,Y)=A^n(1,s): X=a+b s, Y=c+d s
    P=a+I*c; Q=b+I*d; Pb=a-I*c; Qb=b-I*d
    G=[P+s*Q for s in seeds]; Gb=[Pb+s*Qb for s in seeds]
    beta=t+2*I; betab=t-2*I
    edges=[(0,1),(0,2),(1,2),(0,3)]
    # log g_j derivative wrt t and s_j
    dlt=[sp.diff(G[j],t)/G[j]-sp.diff(Gb[j],t)/Gb[j] for j in range(4)]
    dls=[Q/G[j]-Qb/Gb[j] for j in range(4)]
    dl0=1/beta-1/betab  # d/dt log g0, g0=beta/betab
    cols=[]
    cols.append([dlt[j]+dlt[k]-dl0 for (j,k) in edges])
    for m in choose:
        cols.append([(dls[m] if m in (j,k) else 0) for (j,k) in edges])
    M=sp.Matrix(4,4,lambda i,k: cols[k][i])
    Dl=sp.cancel(sp.together(M.det()))
    num,den=sp.fraction(Dl)
    num=sp.Poly(sp.expand(num),t,extension=I) if num.has(I) else sp.Poly(num,t)
    den=sp.Poly(sp.expand(den),t,extension=I) if den.has(I) else sp.Poly(den,t)
    if num.is_zero: return None
    ordinf=den.degree()-num.degree()
    # charged places: roots of z_0 * prod z_k ratio structure = roots of beta,betab,G_j,Gb_j
    charged=sp.Poly(sp.expand(beta*betab*sp.prod(G)*sp.prod(Gb)),t,extension=I)
    sqf=sp.sqf_part(charged)
    # poles simple and at charged places: den | sqf  (den squarefree and divides)
    q,r=sp.div(sqf,den)
    return dict(n=n,choose=choose,ordinf=ordinf,numdeg=num.degree(),dendeg=den.degree(),
                den_divides_charged=r.is_zero, den_squarefree=sp.degree(sp.gcd(den,den.diff(t)),t)==0,
                ncharged=sqf.degree(), h=4*n, W=2*(8*n+1))
for n in [1,2]:
    seeds=[sp.Rational(1),sp.Rational(3),sp.Rational(5),sp.Rational(8)]
    for choose in itertools.combinations(range(4),3):
        print(run(n,seeds,choose)); sys.stdout.flush()
