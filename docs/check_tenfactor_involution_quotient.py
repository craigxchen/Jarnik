"""Sparse exact polynomial identities for the ten-factor quadratic quotient."""

from fractions import Fraction as Q

from check_tenfactor_real_outercut_certificate import A,C,E,R,S,X0,mv

NVAR=6  # t,b,c,d,e,x
ZERO=(0,)*NVAR


def const(c):
    return {} if c==0 else {ZERO:Q(c)}


def var(i):
    e=[0]*NVAR;e[i]=1
    return {tuple(e):Q(1)}


def add(*ps):
    out={}
    for p in ps:
        for m,c in p.items():out[m]=out.get(m,0)+c
    return {m:c for m,c in out.items() if c}


def scale(p,c):
    return {m:a*c for m,a in p.items() if a*c}


def mul(p,q):
    out={}
    for m,c in p.items():
        for n,d in q.items():
            key=tuple(a+b for a,b in zip(m,n))
            out[key]=out.get(key,0)+c*d
    return {m:c for m,c in out.items() if c}


def power(p,n):
    out=const(1)
    for _ in range(n):out=mul(out,p)
    return out


def neg(p):return scale(p,-1)


def main():
    t,b,c,d,e,x=map(var,range(6))
    u=add(t,neg(b),c)
    v=add(t,neg(b),neg(d),e)
    w=add(b,neg(c),d,neg(e))
    q=add(scale(t,2),neg(b),neg(d))
    a=(u,v,w,q,x,b,c,add(t,neg(x)),d,e)
    assert all(add(*(scale(a[j],A[i][j]) for j in range(10)))=={} for i in range(4))
    # Four pivot coordinates have independent coefficients in the linear
    # system, so the displayed six-dimensional kernel is the entire kernel.
    from check_tenfactor_real_outercut_certificate import inverse
    piv=tuple(tuple(A[i][j] for j in range(4)) for i in range(4))
    inv=inverse(piv)
    assert all(mv(inv,tuple(row[j] for row in piv))==tuple(Q(i==j) for i in range(4))
               for j in range(4))
    cubes=tuple(power(z,3) for z in a)
    F=tuple(add(*(scale(cubes[j],A[i][j]) for j in range(10))) for i in range(4))
    H1=add(power(u,3),neg(power(v,3)),neg(power(c,3)),neg(power(d,3)),power(e,3))
    H2=add(power(u,3),power(v,3),power(w,3),neg(power(q,3)),neg(power(d,3)))
    H3=add(power(b,3),neg(power(c,3)),power(d,3),neg(power(e,3)),neg(power(w,3)))
    bigB=add(power(v,3),neg(power(q,3)),power(t,3),neg(power(e,3)))
    relation=add(bigB,scale(mul(t,add(power(x,2),neg(mul(t,x)))),3))
    assert F[0]==relation
    assert F[1]==add(F[0],H1)
    assert F[2]==H2
    assert F[3]==add(H2,H3)
    for h in (H1,H2,H3):
        assert all(sum(m)==3 and m[-1]==0 for m in h)
    assert all(m[0]==0 for m in H3)
    assert max(m[0] for m in H1)==2
    leading={tuple([0]+list(m[1:])):z for m,z in H1.items() if m[0]==2}
    assert leading==scale(add(c,d,neg(e)),3)
    h=add(c,d,neg(e));k=add(c,neg(d),e);y=add(t,neg(b))
    assert h==add(u,neg(v)) and k==add(b,neg(w))
    assert H1==scale(add(mul(h,mul(y,add(y,k))),neg(mul(mul(d,e),add(d,neg(e))))),3)
    assert H3==scale(add(mul(k,mul(b,add(b,neg(k)))),
                        neg(mul(mul(add(c,e),add(d,neg(e))),add(c,neg(d))))),3)
    # The coefficient exchange fixes the entire sign-row matrix.
    assert all(row[4]==row[7] for row in S)
    a0=tuple(c0+z for c0,z in zip(C,mv(E,X0)))
    lips=tuple(sum(map(abs,row)) for row in E)
    assert a0[4]+a0[7]-(lips[4]+lips[7])*R>Q(3,4)
    assert abs(a0[4]-a0[7])-(lips[4]+lips[7])*R>Q(17,20)
    print("PASS: complete linear kernel and four exact sparse cubic identities")
    print("PASS: three homogeneous quotient cubics, quadratic in t, and identical-column involution")
    print("PASS: t>0.75 and |2x-t|>0.85 throughout the certified real box")


if __name__=='__main__':
    main()
