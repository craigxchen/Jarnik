# A rotated cubic family that splits at every place above 17

The fixed-basis obstruction in
[the mod-17 note](cubic_norm_descent_mod17_nonsplitting.md) does not extend
to every rational orthogonal rotation of that basis. For the explicit
rotation below, every nonzero rational parameter `r` with `v_17(r)>=3` makes
the residual cubic split over `Q_17` at all six embeddings of `K(i)`.
This is a local statement. It neither proves splitting over `K(i)` nor
constructs a full Boolean profile.

Let `alpha^3-3alpha+1=0`, `K=Q(alpha)`, and start with

```text
w=2-alpha^2,       v=alpha^2+2alpha-2,       delta=3.
```

Take the rational rotation

```text
p=-13/14,       q=-3/14,       p^2+3q^2=1,
W=pw+qv=(5alpha^2-3alpha-10)/7,
V=-3qw+pv=(22-13alpha-11alpha^2)/7,
L=1+W=(5alpha^2-3alpha-3)/7.
```

It preserves the trace-form basis relations. In particular,
`W^2+V^2/3=4`, just as `w^2+v^2/3=4`.
Use `W,V` in the quartic construction of
[the cubic norm note](cubic_norm_descent_root_incidence.md), retaining
`delta=3`. Write its monic residual cubic as `F=P/(t-i)`.

## The analytic coefficient formulas

Put

```text
u=2r/(1+delta*r^2),
s=(1-delta*r^2)/(1+delta*r^2),
h=delta*u^2,       A=2h/s,       B=2u/s.
```

Direct polynomial division gives

```text
F=t^3+(b/L)t^2+(c/L)t+d/L,
b=iL-AL+2uV,
c=1+W(1-2h)-iLA+VB(-h-1+is),
d=-2A+i[1+W(1-2h)-VB(h+1)].                         (1)
```

These expressions are defined for every nonzero rational `r`. At `r=0`
they have a removable singularity and give

```text
F_0=(t-i)(t+i)^2.
```

At a place above 17 put `e=v_17(L)`. The quantities `W,V,i` are integral,
and `u,s,h,A,B` belong to `Z_17[[r]]`, with `s` a unit series.
Formula (1) therefore proves the coefficientwise inclusion

```text
F-F_0 in r * 17^(-e) Z_17[[r]][t].                   (2)
```

## The six places and their leading discriminants

The three simple roots of the defining cubic modulo 17 are `7,13,14`;
the two simple roots of `x^2+1` are `4,13`. Hensel lifting gives all six
embeddings `K(i) -> Q_17`. For either choice of `i`, one has

```text
alpha mod 17             7       13       14
alpha mod 17^2          75      132       82
L mod 17                 0        3        0
v_17(L)                  1        0        1
L/17 mod 17, if e=1     11        -        6.          (3)
```

The leading tangent quadratic, obtained by putting `t=-i+rz`, is

```text
F(-i+rz)/r^2 at r=0
  =-2iz^2-4(1+i)Vz/L-16delta/L.
```

Its discriminant and the leading coefficient of the full cubic
discriminant are respectively

```text
D=32i[V^2-4delta L]/L^2,
C=-4D=-128i[V^2-12L]/L^2,
Disc(F)=r^2 C+O(r^3).                                (4)
```

The last coefficient can also be derived by observing that the third
root tends to `i`: its two squared differences from the roots tending
to `-i` contribute `(-2i)^4`, and the tangent quadratic has leading
coefficient `-2i`, giving the factor `-4` in (4).

For `i=4 mod 17`, the numerator `-128i[V^2-12L]` has residues
`16,4,16` at the three places in (3). All are nonzero squares.
For `i=13 mod 17` the residues are their negatives, also nonzero
squares because `-1` is a square modulo 17. Consequently, at every
place,

```text
v_17(C)=-2e,       C is a nonzero square in Q_17.       (5)
```

These are statements about the actual 17-adic coefficients. The
residue computations do not assert that `C` equals its residue.

## A whole parameter disk splits

The discriminant of a monic cubic `t^3+bt^2+ct+d` is

```text
b^2 c^2-4c^3-4b^3d-27d^2+18bcd.
```

It has total degree at most four in its three coefficients. Since
all coefficient series of (1) lie in `17^(-e) Z_17[[r]]`, every
coefficient of its discriminant series lies in `17^(-4e) Z_17`.
Using the zero constant and linear terms and the quadratic term (4)
therefore gives the rigorous remainder bound

```text
Disc(F)/r^2-C in r * 17^(-4e) Z_17[[r]].              (6)
```

Suppose `r!=0` and `v_17(r)=k>=3`. Since `e` is either zero or one,
`k-4e>-2e`. Equations (5)-(6) imply that `Disc(F)/r^2` is `C`
times a unit congruent to one modulo 17. Such a unit is a square
in `Q_17`, so `Disc(F)` is a nonzero square.

Equation (2) also gives `F=F_0 mod 17`. The simple root `i` of `F_0`
lifts by Hensel's lemma to a root of `F` in `Q_17`. Factoring out
that root leaves a quadratic whose discriminant differs from
`Disc(F)` by a nonzero square. That quadratic therefore splits too.

Thus `F` splits completely at every one of the six places for every
nonzero rational `r` with `v_17(r)>=3`; `r=17^3` is one explicit
parameter. A test at the fixed prime 17 cannot exclude every rotated
member of this family.
Other primes or global arithmetic may still obstruct splitting over
`K(i)`.

## A separate split-prime fixture

The rational rotation and parameter

```text
p=-971/973,       q=36/973,       r=102
```

give a residual cubic with three distinct roots modulo 197 at every
combination of `alpha=34,169,191` and `i=14,183`. The denominators of
the monic coefficients are units at these places. Hensel's lemma thus
gives splitting over `Q_197` at all six places. Since both 17 and 197
are congruent to 17 modulo 36, a uniform exclusion at split primes in
that residue class cannot follow just from the unit-parameter table at
17. This fixture makes no assertion of global splitting.

The dependency-free
[checker](check_cubic_rotated_mod17_local_escape.py) verifies the exact
rotation identities, Hensel residues, valuations in (3), and the square
residues used in (5). It also checks the discriminant valuations
`4,6,4` for `r=17^3` at both choices of `i`, and the six simple split
reductions at 197. The analytic argument above, rather than a finite
sample of parameters, proves the assertion for the whole parameter disk.
