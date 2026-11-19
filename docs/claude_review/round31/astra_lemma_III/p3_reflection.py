"""Exact saturated exclusion of the double-transposition linear template.

Polynomial coefficients over ZZ or GF(p) use python-flint. The fixed Q(sqrt5)
fibre uses SymPy's exact algebraic field. Triangle order and cofactors are fixed
explicitly, making both projections reproducible. Run without arguments.
"""
from itertools import combinations
from flint import fmpz_mpoly_ctx, fmpz_mod_mpoly_ctx
from p3_two_height import determinant


def geometry(ctx):
    k, v = ctx.gens()
    one = ctx.constant(1)
    labels = (one, -one, k, -k, 0*one)
    edges = list(combinations(range(5), 2))
    real = {(i,j):(labels[i]+labels[j])*(1-v*labels[i]*labels[j]) for i,j in edges}
    X, Y = real[0,2], real[0,3]
    # All imaginary parts are this coefficient pair dotted with (a,b).
    A,B = (k,0*k),(0*k,k)
    C,D = (-k*X,-k*Y),(-X,Y)
    def add(*terms):
        return tuple(map(sum,zip(*terms)))
    def scale(c,pair):
        return c*pair[0], c*pair[1]
    minus_AB = scale(-1,add(A,B))
    E,F = add(minus_AB,C,scale(2,D)),add(minus_AB,scale(2,C),D)
    amps = [[E,A,B,C],[B,A,C],[F,D],[D]]
    imag = {e:amps[e[0]][e[1]-e[0]-1] for e in edges}
    return labels,real,imag


def real_triangle_rows(ctx, labels, real, imag):
    rows = []
    def edge(i,j):
        return tuple(sorted((i,j)))
    def symmetric_product(u,v):
        return u[0]*v[0], u[0]*v[1]+u[1]*v[0], u[1]*v[1]
    for tri in combinations(range(5),3):
        l,m = sorted(set(range(5))-set(tri))
        pair = [[ctx.constant(0) for _ in range(4)] for _ in range(2)]
        for i,j,s in (tri,tri[1:]+tri[:1],tri[2:]+tri[:2]):
            e,f,g=edge(i,j),edge(s,l),edge(s,m)
            c=labels[j]-labels[i]
            x,y,z=real[e],real[f],real[g]
            ab=symmetric_product(imag[e],imag[f])
            ac=symmetric_product(imag[e],imag[g])
            bc=symmetric_product(imag[f],imag[g])
            pair[0][0]+=c*(x*y+x*z+y*z)
            pair[1][0]+=c*x*y*z
            for j in range(3):
                pair[0][j+1]+=c*(ab[j]+ac[j]-bc[j])
                pair[1][j+1]+=c*(z*ab[j]+y*ac[j]-x*bc[j])
        rows.extend(pair)
    return rows


def remove_common(ctx,values):
    common=ctx.constant(0)
    for v in values:
        if v: common=common.gcd(v) if common else v
    assert common
    k,_=ctx.gens()
    assert all(v in (k,k-1,k+1) for v,_ in common.factor()[1]), common
    return [v//common for v in values]


def certificate_equations(ctx):
    labels,real,imag=geometry(ctx)
    raw=real_triangle_rows(ctx,labels,real,imag)
    rows=[remove_common(ctx,raw[j]) for j in (0,1,4,8)]
    w=remove_common(ctx,[(-1)**j*determinant([row[:j]+row[j+1:] for row in rows[:3]])
                         for j in range(4)])
    equations=[w[1]*w[3]-w[2]**2]
    for i in (0,2):
        for power in (3,4,5):
            total=0
            for e,x in real.items():
                p,q=imag[e]
                a_ye=p*w[1]+q*w[2]
                ye_squared=p**2*w[1]+2*p*q*w[2]+q**2*w[3]
                if power==3:
                    val=a_ye*(3*x**2*w[0]-ye_squared)
                elif power==4:
                    val=x*a_ye*(x**2*w[0]-ye_squared)
                else:
                    val=a_ye*(5*x**4*w[0]**2-10*x**2*w[0]*ye_squared+ye_squared**2)
                total+=(int(i in e)-int(4 in e))*val
            equations.append(total)
    return determinant(rows), equations, w


def project(ctx):
    k,v=ctx.gens()
    determinant0, equations, w=certificate_equations(ctx)
    product=ctx.constant(1)
    for f,_ in determinant0.factor()[1]:
        if f.degrees()[1]==0:
            if f not in (k,k-1,k+1):product*=f
            continue
        gcd=ctx.constant(0)
        for polynomial in equations:
            r=f.resultant(polynomial,1)
            if r:gcd=gcd.gcd(r) if gcd else r
        assert gcd
        for factor,_ in gcd.factor()[1]:
            if factor not in (k,k-1,k+1):product*=factor
    p=ctx.constant(1)
    for factor,_ in product.factor()[1]:
        p*=factor
    degree=p.degrees()[0]
    reciprocal=sum(int(c)*k**(degree-i[0]) for i,c in p.to_dict().items())
    common=p.gcd(reciprocal)
    assert common.total_degree()==4 and not common % (k**4-3*k**2+1)
    return degree,common,determinant0,equations,w


def golden(ctx,determinant0,equations,w):
    import sympy as sp
    k,v=ctx.gens()
    t=sp.Symbol("t")
    field=sp.QQ.algebraic_field(sp.sqrt(5))
    phi=(1+sp.sqrt(5))/2
    def special(poly):
        reduced=poly % (k*k-k-1)
        coefficients={}
        for (i,j),c in reduced.to_dict().items():
            coefficients[j,]=coefficients.get((j,),0)+int(c)*phi**i
        return sp.Poly.from_dict(coefficients,(t,),domain=field)
    gcd=special(determinant0)
    for f in equations:
        gcd=gcd.gcd(special(f))
    witness=sp.Poly(t-(5-sp.sqrt(5))/4,t,domain=field)
    assert gcd.monic()==witness
    remainders=[special(f).rem(witness) for f in w]
    assert remainders[0].is_zero and not remainders[1].is_zero
    print("golden-ratio fibre",witness.as_expr(),"Cramer residues",[f.as_expr() for f in remainders])


if __name__=='__main__':
    ctx=fmpz_mpoly_ctx.get(["k","v"])
    degree,reciprocal,det,eqs,cof=project(ctx)
    print("rational projection degree",degree,"reciprocal gcd",reciprocal)
    golden(ctx,det,eqs,cof)

if __name__=='__main__':
    for characteristic in (101,109,139,149):
        ctx=fmpz_mod_mpoly_ctx.get(["k","v"],characteristic)
        degree,common,det,eqs,w=project(ctx)
        k,v=ctx.gens()
        root=next(z for z in range(characteristic) if (z*z-z-1)%characteristic==0)
        g=det%(k-root)
        for f in eqs:
            g=g.gcd(f%(k-root))
        assert g.total_degree()==1
        assert not (w[0]%(k-root))%g and (w[1]%(k-root))%g
        print("characteristic",characteristic,"projected degree",degree,"excluded final fibre",g)
