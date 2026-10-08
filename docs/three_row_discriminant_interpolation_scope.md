# Three-row constant-imaginary interpolation and its arithmetic scope

This note isolates an exact three-row consequence of a common oriented
Gaussian block. It gives a finite full-factor fixture and identifies two
limits of a discriminant or Vieta argument. It proves no asymptotic
exclusion or radius-uniform count.
The finite identities are checked by
[check_three_row_discriminant_interpolation_scope.py](check_three_row_discriminant_interpolation_scope.py).

The two fixtures below have a ramified factor above 2 in some rows.
Their ordinary coordinate gcds are one, but those rows are not coprime
to their Gaussian conjugates. The actual primitive-anchor reduction
requires this stronger condition. A separate
[conjugate-primitive seven-block fixture](finite_full_support_conjugate_primitive_configuration.md)
satisfies it exactly. None of these finite fixtures is an
asymptotically balanced endpoint family.

Let `Y_1=Y_2=Y_3=1` and

```text
P_i=X_i+i,       X_i=r+g a_i,
g>0,             g | r^2+1,
s=(r^2+1)/g,    f(t)=g t^2+2r t+s.
```

Then `Norm(P_i)=g f(a_i)` and `disc(f)=-4`. Suppose the three pair
differences are pairwise coprime and the pair-private norm for `ij`
is `|a_j-a_i|`, so that the signed pair residue is a unit.
The pair factors incident to row `i` imply

```text
q_i=f(a_i) / product_(j!=i)(a_i-a_j) in Z.
```

Lagrange interpolation gives the *whole* quadratic identity

```text
f(t)=sum_i q_i product_(j!=i)(t-a_j).                 (1)
```

Order the nodes `a_1<a_2<a_3`, and put `A=a_2-a_1`,
`B=a_3-a_2`, `C=A+B`. Positivity of `f` makes
`(q_1,q_2,q_3)=(u,-v,w)` with `u,v,w>0`. Recenter the
interpolation variable at `a_1`; its linear half-coefficient is
`r'=r+g a_1`. Comparing coefficients in (1) gives

```text
g=u-v+w,
2r'=w B-u A-g C,
(r')^2+1=g u A C.                                       (2)
```

The last equation is precisely `disc(f)=-4`. Viewed as a quadratic
equation for `v`, its other rational root is

```text
v'=2(u B+w A)/C-v.                                      (3)
```

There is no general integral Vieta mutation: the numerator in (3)
need not be divisible by `C`. Gauss reduction of the positive
unimodular Gram form `[[g,r],[r,s]]` likewise changes the marked
nodes, so the reduced form's bounded coefficients do not bound the
original `a_i`.

## An exact finite factor fixture

Take

```text
g=17, r=4, s=1,      (a_1,a_2,a_3)=(1385,1482,1747),
(A,B,C)=(97,265,362),
(u,v,w)=(929,1453,541),
(X_1,X_2,X_3)=(23549,25198,29703).
```

Direct substitution gives

```text
Norm(P_1)=554555402=17·97·362·929,
Norm(P_2)=634939205=17·97·265·1453,
Norm(P_3)=882268210=17·265·362·541.
```

The seven displayed norm factors are pairwise coprime, and they have
the primitive representations

```text
17=4^2+1^2,   97=9^2+4^2,   265=16^2+3^2,
362=19^2+1^2, 929=23^2+20^2,
1453=38^2+3^2, 541=21^2+10^2.
```

In fact the pairwise gcds of the three
row norms are exactly `17·97`, `17·265`, and `17·362`, in the
corresponding row order. Each of these pair factors also divides
`X_j-X_i`, so the incident rows have the same real coordinate modulo
every common prime power. Since every `P_i=X_i+i` is primitive in
`Z[i]`, each odd split factor chooses a unique orientation, and this
congruence makes the orientation agree in its incident rows. Thus
this is a genuine finite
three-row oriented-factor configuration, with the ramified factor
`2` in the `13` pair. Formula (3) gives
`v'=197+12/181`, explicitly failing integrality.

The fixture is not an asymptotically balanced family: its seven
factors have different fixed sizes, and it supplies no sequence
with all logarithms `w+o(w)`.

A second finite fixture has somewhat closer common and pair scales:

```text
g=109, r=-33, s=10,      (a_1,a_2,a_3)=(134,291,304),
(A,B,C)=(157,13,170),   (u,v,w)=(73,4513,4549),
(X_1,X_2,X_3)=(14573,31686,33103),
f(t)=(10t-3)^2+(3t-1)^2.
```

Its row norms are `212372330`, `1004002597`, and `1095808610`;
the pairwise norm gcds are exactly `109·157`, `109·13`, and
`109·170`. The seven odd split-core norms after peeling the
ramified factor `2` from `C` are
`109,157,13,85,73,4513,4549`, pairwise coprime. Primitive
Gaussian representations are respectively

```text
109=10^2+3^2,  157=11^2+6^2, 13=3^2+2^2,
85=9^2+2^2,    73=8^2+3^2,
4513=48^2+47^2, 4549=65^2+18^2.
```

There is also a direct orientation check. With `H=10+3i`, one has
`r+i=H(-3+i)`, and the three Gaussian row gcds, up to units,
are `H(-11+6i)`, `H(3-2i)`, and `H(-11-7i)` in the `12,23,13`
order. Their norms are the preceding pair gcds. This remains a
finite fixture; fixed ratios among its seven norms give no
balanced infinite family.

## Two scale cautions

If all three positive differences `A,B,C=A+B` are primitive
Gaussian norms and pairwise coprime, then `A,B` are odd and
`1 mod 4`, while `C` is `2 mod 4`. Indeed a primitive odd norm is
`1 mod 4` and a primitive even norm is `2 mod 4`; placing the
even norm at `A` or `B` would make the remaining odd norm
`3 mod 4`. Therefore a strict model in which all three pair
norms have only odd split-prime support and all pair residues
are exactly one is impossible. The unavoidable factor `2` can
instead be a bounded ramified residual, as in the fixture.

Also, centering `r` by `|r|<=g/2` does not imply `s` has the same
scale as `g`. For every integer `t>=2`,

```text
g=2t^2+2t+1=Norm((t+1)+it),
r=-(2t+1),          s=2,
r+i=(-1+i)((t+1)+it),
r^2+1=g s.
```

Here `g` tends to infinity and `|r|<=g/2`, but `s` is constant.
The row scale `f(a_i)~g^3` remains possible when `|a_i|~g`.
