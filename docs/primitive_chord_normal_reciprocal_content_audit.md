# Primitive-chord normal: the endpoint square and its reciprocal content

This is an exact arithmetic audit of the integer normal to an actual
endpoint chord.  It refines the [primitive-chord strip transform](endpoint_primitive_chord_transform.md)
by recording a divisibility hidden in the endpoint data and by distinguishing
the content of two endpoint rows from the content of the whole transformed
tuple.  The resulting square equation does not give a small-root Runge
relation or a radius descent for a general arc.

Let `P,Q,z_i` be distinct integer points of squared norm `N>0` on a minor
circle arc of length at most `C N^(1/4)`, with `P,Q` its extreme points.
Write

```text
Q-P=g w,     w=u+i v primitive,     s=u^2+v^2,
G=g s,       B=det(w,P),            T=G/2.
```

The endpoint norm equality gives `w dot P=-T`, so `G` is even.  For any
source row `z_i=P+h_i`, use the exact index-`s` coordinates

```text
p_i=w dot h_i,        q_i=det(w,h_i).
```

The inverse frame is

```text
h_i=(u p_i-v q_i, v p_i+u q_i)/s,                    (1)
```

and integrality is equivalent to both numerator coordinates being divisible
by `s`.  These are the same coordinates as in the cited strip transform.

## 1. The endpoint divisibility makes the first residue test automatic

Both `u` and `v` are units modulo `s` (for `s=1`, read the congruences as
vacuous): primitivity gives `gcd(u,s)=gcd(v,s)=1`.  Applying (1) to `P`
itself, with dot coordinate `-G/2` and determinant coordinate `B`, and
multiplying by two gives

```text
s | 2B.                                                (2)
```

Set `t=2B/s`, an integer.  The endpoint norm and Gaussian product become

```text
4N=s(g^2+t^2),
P=w(-g+i t)/2,       Q=w(g+i t)/2.                    (3)
```

The displayed halves are Gaussian integers because they are the actual
endpoint rows.  Formula (3) keeps any common Gaussian factors and their
orientations; `s` is not automatically an ordinary scalar content of the
whole tuple.

For `s>1`, put `k=v u^(-1) mod s`.  Then (1) reduces to

```text
p_i=k q_i mod s,       k^2+1=0 mod s.                 (4)
```

The exact circle conic is

```text
p_i^2+q_i^2-G p_i+2B q_i=0.                          (5)
```

Reducing (5) modulo `s` after (4) gives `2B q_i=0 mod s`, which holds
for **every** integer `q_i` by (2).  Thus the tempting congruence does not
force the small strip coordinate into a proper residue class.  Higher
congruences and the full integral conic can still constrain actual rows;
only this first-order reduction is redundant.

## 2. The two equal projections give a nonzero square constant

Let the primitive normal be `v_normal=(-v,u)` and put

```text
x_i=v_normal dot z_i=B+q_i,
y_i=det(v_normal,z_i)=T-p_i.
```

The endpoint rows have `(x,y)=(B,T)` and `(B,-T)`.  Every row satisfies

```text
x_i^2+y_i^2=sN,
y_i^2=T^2-2B q_i-q_i^2.                          (6)
```

The proposed common-projection model is therefore exactly (6), with
`K=T^2`, `A=-B`, and signed shift `q_i`.  Its reciprocal factorization is

```text
(y_i-T)(y_i+T)=-q_i(2B+q_i),
p_i(G-p_i)=q_i(2B+q_i).                              (7)
```

The second line is (5) rearranged: the two endpoint square roots yield no
independent divisor modulus.  The strip shifts satisfy
`|q_i|<=W sqrt(s)`, where the actual sagitta `W<=C^2/8`, but the linear
coefficient is large:

```text
|B|^2=sN-T^2,
T^2=g^2s^2/4<=C^4 N/(4g^2),
|B|=sqrt(sN) sqrt(1-|Q-P|^2/(4N)).                 (8)
```

In the difficult primitive-chord scale `g=1`, `s` can be comparable to
`sqrt(N)`.  The geometric allowance for `|q_i|` is then of order
`N^(1/4)`, while `|B|` has order `N^(3/4)`; the corresponding raw
monic-root allowance for `q_i(2B+q_i)=T^2-y_i^2` is of order `N`.
The alternative small roots `q_i` do **not** occur as factors
`y_i^2=q_i(Z-q_i)`: the nonzero constant `T^2` remains.  An even
squareclass relation among the `q_i` alone therefore does not turn (6)
into the small-root monic relation used by the
[elementary Runge bound](runge_even_product_square_bound.md).

## 3. Endpoint content is not whole-tuple content

The ordinary content of the two projected endpoint rows is

```text
D_end=gcd(|B|,T)=(s/2) gcd(g,|t|).                   (9)
```

The right side is an integer even when `s` is odd, since then `g,t`
are even.  If all projected rows are retained, their common ordinary
content is exactly

```text
D_all=gcd(|B|,T, {p_i,q_i : i not P}).              (10)
```

Every genuine interior row has `q_i!=0`, because the chord line meets
the circle only at `P,Q`.  Hence `D_all<=|q_i|<=W sqrt(s)` for each
interior row.  Integer **scalar** normalization of this particular
projected frame has squared norm

```text
N_scalar=sN/D_all^2 >= N/W^2 >=64N/C^4.           (11)
```

For example, the actual `C<2` arc on `N=85`

```text
P=(7,6),   I=(6,7),   Q=(2,9),
w=(-5,3),  s=34,     g=1,       B=-51, t=-3,
(p_I,q_I)=(8,-2),    T=17,
```

has `D_end=17` but `D_all=1`.  Its interior monic root is
`q_I(2B+q_I)=208`, coprime to the square parameter `T^2=289`.
Thus even endpoint content cannot be divided out of both parameter and
interior roots in this literal short arc.

There is a separate Gaussian normalization to keep in view.  The
projected row is just a common Gaussian similarity:

```text
x_i+i y_i=-i conjugate(w) z_i.                     (12)
```

Consequently the Gaussian gcd of the entire projected tuple is
`conjugate(w)` times the Gaussian gcd of the source tuple, up to a unit.
Removing that common Gaussian factor returns the original tuple up to
a quarter-turn.  Inequality (11) concerns ordinary scalar content in
the displayed normal frame, not a new least physical radius after full
Gaussian reduction.

The exact remaining bridge would have to control simultaneous
higher-modulus residues or the squareclasses of the large factors
`2B+q_i`, or construct a different small-root equation that removes
`T^2` while retaining integer parameter height.  Endpoint equality and
the index-`s` congruence alone supply none of these.  The
[checker](check_primitive_chord_normal_reciprocal_content.py) verifies
the identities on parametrized equal-norm endpoint pairs and the
literal three-row short-arc fixture.
