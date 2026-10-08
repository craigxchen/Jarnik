# Central circle moments as a fixed-residue curve

The central moment equations have a direct interpretation as residues of
a rational differential. This identifies relevant literature without
supplying the missing uniform arithmetic height estimate.

## 1. Exact algebraic dictionary

Take `n=2h` distinct nonzero nodes `q_i`, with `q_0=1`, and nonzero weights
`u_i` satisfying

```text
sum_i u_i q_i^j=0,              -(h-1)<=j<=h-1.
```

Let `P(z)=product_i(z-q_i)`. Partial fractions give

```text
omega = sum_i u_i dz/(z-q_i)
      = kappa z^(h-1) dz/P(z),             kappa !=0.       (1)
```

Here is a direct verification. The nonnegative moments force the numerator
of the partial fraction sum to have degree at most `h-1`. The negative
moments force a zero of order at least `h-1` at zero. Nonzero residues
exclude the zero differential, so the numerator is exactly as stated.
The differential has simple poles at the `q_i`, residues `u_i`, and two
zeros of order `h-1`, at zero and infinity.

Conversely, a differential with these two ordered zeros and labeled simple
poles can be put in (1) by a Mobius transformation. Scaling the coordinate
then puts its first pole at one. Comparing Laurent expansions recovers
the displayed moment equations. Thus fixing `u` identifies our complex
moment curve with a fixed-residue fiber in the stratum

```text
H(h-1,h-1; [-1]^(2h)).                               (2)
```

The zero labels matter: inversion exchanges them and sends the normalized
nodes to `1/q_i`. They are not being identified in (2).
For six labels it fixes the rowspace
`span(1,(q+q^(-1))/2,(q-q^(-1))/(2i))`, hence preserves each isotropic
ruling. It acts inside each genus-five component as its covering
involution `d -> -d`; it does not exchange the two components. For
`n=6`, pullback changes `kappa` to `-kappa/product_i q_i` and leaves
the labeled residues unchanged.

For the rational circle problem, all `u_i` are integers and all `q_i` lie
in `Q(i)` with `q_i conjugate(q_i)=1`. Those arithmetic and real-locus
conditions are additional to the complex identification (2).

## 2. Relevant primary source

Chen, Gendron, Prado and Tahar's
[Isoresidual curves, arXiv:2412.16810v3](https://arxiv.org/html/2412.16810v3)
studies these fibers. Definition 1.1 calls the proper-subsum hyperplanes
the resonance arrangement. Theorems 1.5 and 1.7 give their Euler
characteristic and connected components outside that arrangement. For
`H(2,2;[-1]^6)` their formulas give two components of genus five, in
agreement with our independent pencil and simple-branch calculation.
Theorem 1.2 supplies a meromorphic differential on each fiber using
relative periods between the two zeros; its saddle-connection periods
are residue subsums. These are geometric statements, not uniform bounds
for the heights of rational points. The current source was checked on
September 13, 2026.

Their Section 2.1.2 uses period residues, which are `2pi i` times our
coefficient residues: use `lambda_i=2pi i u_i` for (1), or first replace
`omega` by `omega/(2pi i)` to use `lambda=u`. This fixed scaling leaves
the fibers and genus statements unchanged.

For reference, substituting `a1=a2=h-1`, `a=2h-2` in their formula gives

```text
-chi = (2h-2)! * (2h-6+2/h).
```

For `h=2,3,4`, this means respectively one genus-zero component, two
genus-five components, and one genus-901 component. The six-label case
already has the [self-contained proof](six_point_simple_branch_genus.md).

## 3. An exact logarithmic function and its limitation

For ordinary integer weights, the expression

```text
F(q)=product_i q_i^(u_i)                              (3)
```

is a rational function on the normalized moment curve. In a simply
connected coordinate neighborhood avoiding the poles, an antiderivative
of (1) is `sum_i u_i log(z-q_i)`. Its difference between the two ordered
zeros, zero and infinity, equals `-sum_i u_i log(-q_i)`, modulo constant
periods. The term at infinity is zero because `sum u_i=0`.
Differentiating along the fixed-weight fiber consequently gives

```text
d(relative period) = -sum_i u_i dq_i/q_i = -dF/F.      (4)
```

This local identity follows directly from partial fractions; it does not
depend on selecting global logarithm branches. On the circle locus,
`q_i=exp(i theta_i)` locally, so `F=exp(i sum_i u_i theta_i)`.

The integrality of residues and of their subsums does not make
`sum_i u_i theta_i` an integer or an integer multiple of `2pi`. The poles
are Gaussian rational numbers, not generally roots of unity. Nor does
integrality quantize the period from an arbitrary point of the fiber to
a boundary point. Any argument using (4) must retain the height and
degree of (3), which depend on the actual central weights. No
radius-independent estimate follows just from the period structure.

The [Belyi period-map proof](six_point_belyi_period_map.md) now gives the
exact logarithmic differential, all ramification indices, and the map's
degree by a separate algebraic calculation. It also descends the symmetric
trace to a rational function on the isotropic base. These results retain
their dependence on the moving central weights.
