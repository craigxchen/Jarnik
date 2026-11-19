"""Full affine-gauge strict ideal of the split degree-one ten-pair system.

python p3_global_saturation.py --probe       # bounded global first-stage runs
python p3_global_saturation.py --full 101    # exact SymPy fallback, expensive
python p3_global_saturation.py --singular 101 > linear101.sing
Singular -q linear101.sing                    # full saturation, dimension, degree

Characteristic 0 means QQ. The chart u_01=0,v_01=1 is global on the strict
locus over every algebraically closed field: the needed simultaneous affine
source map is unique. It is the split gauge, not a real parametrization.
"""
import argparse
from itertools import combinations
import resource
import time

from flint import fmpz_mpoly_ctx, fmpz_mpoly_vec, fmpz_mod_mpoly_ctx

NAMES=[f"r{j}" for j in range(1,10)]+[f"b{j}" for j in range(5,10)]
EDGES=list(combinations(range(5),2))

def system(ctx):
    one=ctx.constant(1)
    generators=ctx.gens()[-len(NAMES):]
    u=[0*one,*generators[:9]]
    b5,b6,b7,b8,b9=generators[9:]
    # b = u-v.  These five affine functions solve exactly the four first
    # moment equations and b[0]=-1 (the pivots are a unimodular minor).
    b=[-one,
       -one+b5+b6-b7-b8,
       2*one-b5+b7+2*b8+b9,
       2*one-b6+2*b7+b8+b9,
       3*one-b5-b6+2*b7+2*b8+2*b9,
       b5,b6,b7,b8,b9]
    v=[r-d for r,d in zip(u,b)]
    weights=[[int(i in e)-int(4 in e) for e in EDGES] for i in range(4)]
    assert all(not sum(a[j]*b[j] for j in range(10)) for a in weights)
    equations=[]
    for power in range(2,6):
        differences=[r**power-s**power for r,s in zip(u,v)]
        equations.extend(sum(a[j]*differences[j] for j in range(10)) for a in weights)
    roots=u+v
    collisions=[]
    for i,j in combinations(range(20),2):
        d=roots[i]-roots[j]
        if d.total_degree()>0:collisions.append(d)
        else:assert d
    return equations,collisions


def singular(p):
    ctx=fmpz_mpoly_ctx.get(NAMES,"degrevlex")
    equations,collisions=system(ctx)
    print('LIB "elim.lib";')
    print(f'ring r = {p},({",".join(NAMES)}),dp;')
    print("ideal I = "+",\n".join(map(str,equations))+";")
    print("ideal diffs = "+",\n".join(map(str,collisions))+";")
    print("int j;")
    print("for (j=1;j<=size(diffs);j++) { I=sat(I,ideal(diffs[j]))[1]; }")
    print("ideal G=std(I);")
    # Singular dim for the unit ideal is -1. Explicitly use degree zero there.
    print('if (reduce(1,G)==0) { print("empty; dim -1; degree 0"); }')
    print('else { print("dimension"); dim(G); print("degree"); mult(G); }')
    print("quit;")


def probe():
    resource.setrlimit(resource.RLIMIT_AS,(1200*1024**2,1200*1024**2))
    ctx=fmpz_mpoly_ctx.get(NAMES,"degrevlex")
    eq,collisions=system(ctx)
    print(f"variables={len(NAMES)} equations={len(eq)} collision factors={len(collisions)}",
          flush=True)
    print("input terms",list(map(len,eq)),flush=True)
    # Start on the full, un-specialized unsaturated ideal; an exact SAT
    # workflow may use this preliminary basis before saturating. No partial
    # basis or initial monomial set below is used to report dimension.
    for limits in ((48,12000,2048),(96,12000,2048)):
        t=time.monotonic()
        basis,success=fmpz_mpoly_vec(eq,ctx).buchberger_naive(limits=limits)
        print("QQ",limits,"generated",len(basis),"longest",max(map(len,basis)),
              "done",success,"seconds",round(time.monotonic()-t,3),flush=True)
    for p in (101,103,107):
        ctx=fmpz_mod_mpoly_ctx.get(["inv",*NAMES],p,"lex")
        equations,collisions=system(ctx)
        equations.append(ctx.gen(0)*collisions[0]-1)
        # Actual first localization step; eliminate inv after lex GB. A
        # total-degree-priority pair queue is used within the lex order.
        modular_first_stage(ctx,p,equations,48,3000)
        if p==101:
            modular_first_stage(ctx,p,equations,256,20000)


def modular_first_stage(ctx,p,inputs,basis_limit,terms_limit):
    from heapq import heappush,heappop
    basis=[]
    leading=[]
    queue=[]
    def insert(f):
        f=f*pow(int(f.leading_coefficient()),-1,p)
        lead=f.monoms()[0]
        index=len(basis)
        for other,previous in enumerate(leading):
            if not any(a and b for a,b in zip(lead,previous)):
                continue
            lcm=tuple(max(a,b) for a,b in zip(lead,previous))
            heappush(queue,(sum(lcm),other,index,lcm))
        basis.append(f)
        leading.append(lead)
    for f in inputs:
        insert(f)
    processed=0
    longest=max(map(len,basis))
    started=time.monotonic()
    reason="complete"
    while queue:
        _,i,j,lcm=heappop(queue)
        ex=tuple(c-a for c,a in zip(lcm,leading[i]))
        ey=tuple(c-b for c,b in zip(lcm,leading[j]))
        f=ctx.from_dict({ex:1})*basis[i]-ctx.from_dict({ey:1})*basis[j]
        processed+=1
        while f:
            if len(f)>terms_limit:break
            old=f
            for g in basis:
                f=f%g
                if not f or len(f)>terms_limit:break
            if f==old:break
        longest=max(longest,len(f))
        if len(f)>terms_limit:
            reason="term limit"
            break
        if f:
            insert(f)
            if len(basis)>=basis_limit:
                reason="generator limit"
                break
    print(f"GF({p}) first saturation: inputs={len(inputs)} pairs={processed} "
          f"basis={len(basis)} longest={longest} pending={len(queue)} "
          f"stop={reason} seconds={time.monotonic()-started:.3f}",flush=True)


def dimension_degree(monomials, n):
    """Hilbert series of a monomial ideal via the lcm-lattice numerator."""
    from collections import defaultdict
    from math import comb
    zero=(0,)*n
    if zero in monomials:
        return -1,0
    # Remove the nonminimal leading monomials before inclusion-exclusion.
    ordered=sorted(set(monomials),key=sum)
    minimal=[]
    for u in ordered:
        if not any(all(a<=b for a,b in zip(v,u)) for v in minimal):
            minimal.append(u)
    coeff={zero:1}
    for u in minimal:
        update=defaultdict(int,coeff)
        for v,m in coeff.items():
            lcm=tuple(max(a,b) for a,b in zip(u,v))
            update[lcm]-=m
        coeff={v:m for v,m in update.items() if m}
    for order in range(n+1):
        val=sum(m*comb(sum(v),order) for v,m in coeff.items() if sum(v)>=order)
        if val:
            return n-order, (-1)**order*val
    raise ArithmeticError("zero Hilbert numerator for a non-unit ideal")


def full(p):
    """Exact reference calculation using only SymPy with python-flint."""
    import sympy as sp
    ctx=fmpz_mpoly_ctx.get(NAMES,"degrevlex")
    equations,collisions=system(ctx)
    gens=sp.symbols(" ".join(NAMES))
    inverse=sp.Symbol("z_inverse")
    domain=sp.QQ if p==0 else sp.GF(p)
    def expr(f):
        return sp.Poly.from_dict({k:int(c) for k,c in f.to_dict().items()},
                                 gens,domain=domain).as_expr()
    polys=[expr(f) for f in equations]
    for i,g in enumerate(collisions):
        gb=sp.groebner(polys+[inverse*expr(g)-1],inverse,*gens,domain=domain,
                       order="lex",method="f5b")
        polys=[h.as_expr() for h in gb.polys if h.degree(inverse)==0]
        if polys==[1]:
            print("empty; dimension -1; degree 0")
            return
        print("saturated",i+1,"of",len(collisions),flush=True)
    gb=sp.groebner(polys,*gens,domain=domain,order="grevlex",method="f5b")
    leading=[f.LM(order="grevlex").exponents for f in gb.polys]
    print("dimension and degree",dimension_degree(leading,len(gens)))


if __name__=="__main__":
    ap=argparse.ArgumentParser(description=__doc__)
    actions=ap.add_mutually_exclusive_group(required=True)
    actions.add_argument("--probe",action="store_true")
    actions.add_argument("--singular",type=int,metavar="P")
    actions.add_argument("--full",type=int,metavar="P")
    args=ap.parse_args()
    if args.probe:probe()
    elif args.full is not None:full(args.full)
    else:singular(args.singular)
