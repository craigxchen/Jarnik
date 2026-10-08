# Two continuous local bounds for the standalone joint involution

For the old cuts `16` and `25`, the six transformed points, without the
original anchor, have common Gaussian matching scalar contribution at least

```text
(e/2-2E_p) log p.                                       (1)
```

Here `e` is the old core exponent and `E_p` is the total correction budget
from [joint_involution_five_cut_inflation.md](joint_involution_five_cut_inflation.md),
equation (5). This is a continuous valuation argument: it permits arbitrary
extra cancellations and does not round valuations to a finite grid.

## 1. The standalone sorted formula

Sort the six signed phase valuations as `d_1<=...<=d_6`. With the anchor
removed, their local radius width is `d_6-d_1`. The common matching scalar
is at least

```text
Gamma = ((d_6-d_5)+(d_6-d_4)+(d_3-d_1)+(d_2-d_1))/2.      (2)
```

Thus four phases at most `a` and one at least `b` give
`Gamma>=(b-a)/2`; four phases at least `b` and one at most `a` give the
same bound. Another useful form is

```text
2 Gamma=(d_6-d_1)+(d_6-d_5)+(d_2-d_1)+(d_3-d_4).         (3)
```

If two low phases are at most zero, two middle phases are equal, and two
high phases are at least `e`, (3) gives `Gamma>=e/2`.

All brackets, forms, and Gaussian lifts below use the exact conventions
and two binomial identities in the five-cut note.

## 2. Cut `16`, with no correction error

The fixed phases are `(e,0,0)`. Put `v=v_p(V)`. The alternate expressions
for `Q_4,Q_5` have the common pattern

```text
Q_i=V P_i-R_i P_1,                                     (4)
```

where both scalar coefficients `R_i` are units at this cut. If `v=0`,
the first term has first Gaussian valuation zero and the second valuation
`e`. Thus `d_4,d_5<=0`, and (2) applies to four low phases.

If `v>0`, both `H` and `F` are units: use `V-H=2A` and
`V-F=2(A+B)`, where `A` and `A+B` are unit matching products. In the
primary expression for `Q_6`, the two summands have Gaussian orders
`(e,0)` and `(e,e)`, so `d_6>=e`.

The two summands in (4) have orders `(v,v)` and `(e,0)`.
Consequently:

* if `0<v<e`, then `d_4=d_5=v`;
* if `v>e`, then `d_4=d_5=e`;
* if `v=e`, then `d_4,d_5>=e`.

The first case uses (3), and the latter cases use four high phases in
(2). In every case `Gamma>=e/2`.

## 3. Cut `25`, with no correction error

The fixed phases are `(0,e,0)`. Put `h=v_p(H)`. The primary expressions
for `Q_4,Q_6` have the pattern

```text
Q_i=H P_i+/-R_i P_2,                                  (5)
```

where both `R_i` are units. If `h=0`, these give `d_4,d_6<=0`, and the
old rows `1,3` supply the other two low phases.

If `h>0`, the identities `H-F=2B` and `V-H=2A` show that `F,V` are
units. The alternate expression

```text
Q_5=V P_5-2 Delta_35 Delta_46 Delta_52 P_1              (6)
```

has summand orders `(e,0)` and `(e,e)`, giving `d_5>=e`. The two phases
from (5) equal `h` when `0<h<e`, equal `e` when `h>e`, and are both at
least `e` when `h=e`. Equations (2)--(3) again prove `Gamma>=e/2`.

## 4. A uniform error transfer

Write `E=E_p`. The following proof applies to both cuts, with `h=v_p(V)`
for `16` and `h=v_p(H)` for `25`. Call the two moving outside rows the
middle rows, and the remaining moving row the inside row.

Allow the more general local data

```text
v_pi(P_i)=e 1_(i in S)+alpha_i,   v_barpi(P_i)=beta_i,
v_p(Delta_ij)=e 1_(i,j in S)+delta_ij,
alpha_i,beta_i,delta_ij>=0,
E=sum_i(alpha_i+beta_i)+sum_ij delta_ij+v_p(2).
```

In particular both orientation corrections may be positive; local row
primitivity is unnecessary. An old row inside the cut has signed phase
in `[e-E,e+E]`, and an outside old row has phase in `[-E,E]`. A scalar
unit in the idealized proof has valuation at most `E`, and a scalar with
idealized order `e` has order in `[e,e+E]`. The total budget bounds
the sum of the relevant bracket and row corrections, not merely each
one separately. These hypotheses are preserved, with exactly the same
`E`, by the complement reflection in
[joint_involution_standalone_local_symmetry.md](joint_involution_standalone_local_symmetry.md).

Assume `e>4E`. If `h<=2E`, the first-orientation order of the first term
in each middle-row binomial is at most `3E`, whereas the other is at least
`e`. Hence the two middle phases are at most `3E`. Together with the two
old outside phases, these give four lows at most `3E` and a high at least
`e-E`. Thus `Gamma>=e/2-2E`.

Now suppose `h>2E`. The two complementary forms used above have
valuation at most `E`, by their differences from unit matching products.
The moving inside phase is at least `e-E`: its conjugate first-term
order is a complementary-form order plus `beta_inside`, whose sum is
at most `E`, whereas the second term has order at least `e`. For each
middle row the conjugate orientation is `rho_i+beta_reference<=E`,
where `rho_i` is the scalar coefficient's correction valuation. It
cannot cancel with the first term, whose conjugate order exceeds `2E`.

* If `h<e-2E`, the first term strictly dominates in the first orientation.
  Both middle phases lie in `[h-E,h+E]`. They are strictly above the two
  old outside phases and below both inside phases. In (3), the first
  difference is at least `e-2E`, the middle gap `d_3-d_4` is at least
  `-2E`, and the other two differences are nonnegative. Thus
  `Gamma>=e/2-2E`.
* If `h>e+2E`, the second term strictly dominates in both orientations.
  Both middle phases equal the same old inside phase, which is at least
  `e-E`. Four phases are at least `e-E`, and an old outside phase is at
  most `E`. Hence `Gamma>=e/2-E`.
* If `|h-e|<=2E`, possible first-orientation cancellation only increases
  its order. Each middle phase is at least `e-3E`. There are therefore
  four phases at least `e-3E` and an old outside phase at most `E`, giving
  `Gamma>=e/2-2E`.

Finally, if `e<=4E`, the nonnegative bound `Gamma>=0` already gives (1).
All prime powers and arbitrarily high cancellation depths are covered.

## Scope and audit

The root agent independently supplied the idealized `16` proof; the
uniformity-audit agent supplied `25`, the common error estimate, and this
write-up. The fresh-algebraic and symmetry agents independently checked
the error proof and its extension to both orientation corrections,
required for the exact complement symmetry. The full combination with
the other fourteen cut bounds and the actual radius normalization is
[joint_involution_standalone_inflation.md](joint_involution_standalone_inflation.md).
It proves that this particular joint map loses endpoint scale even when
the old anchor is omitted. Neither local statement assumes that transformed
residues remain small or that the output retains a uniform profile.
