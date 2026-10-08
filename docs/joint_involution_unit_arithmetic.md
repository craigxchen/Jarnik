# Arithmetic constraints on the joint six-point involution

This independently audits the map in
[joint_segre_involution.md](joint_segre_involution.md). It proves an exact
projective-unit invariance and a radius inequality for any circle realization.
Including the exceptional boundary contents gives the actual output-height
range `42w-o(w) <= h_out <= 48w+o(w)`. The sharper generic value `48w`
is not asserted. No uniform circle-point bound is proved.

## 1. Arithmetic domain and five preserved characters

In the chart `(infinity,0,1,a,b,c)`, put

```text
A=(a-1)(c-b), B=(c-1)(b-a), C=c-b, D=b-a, E=b(a-1),
L=A+B+C+D+E=b(c-1),
F=c(a+1-b)-a,             G=c(a+b-1)-a,
H=c(b+1-a)+a-2b,          U=c(b+1-a)-a,
V=c(a+b-1)-a(2b-1).
```

The joint involution has the exact formulas

```text
a'=aH/V, b'=bF/V, c'=cH/F,
a'-1=-(a-1)F/V, b'-1=-(b-1)U/V, c'-1=(c-1)U/F,
a'-b'=-(a-b)G/V,
a'-c'=(a-c)HG/(FV), b'-c'=(b-c)UG/(FV).                  (1)
```

Independently, the dual-coordinate deck formula factors as
`w'_0=w_0(G/F)^2`. Every exceptional form is a sum or difference of
two individual old perfect-matching monomials:

```text
F=C-B, H=B+C,
U=(B+D)+(C+D), V=(A+B)+(A+C), G=(L-B)+(L-C).               (2)
```

Here `B+D=c(b-a)`, `C+D=c-a`, `A+B=(b-1)(c-a)`,
`A+C=a(c-b)`, `L-B=a(c-1)`, and `L-C=c(b-1)` are themselves
matching monomials. Each paired pair of matchings has no common edge.
Its ratio is an alternating six-cycle projective unit, of height
`12w+o(w)` in the old full `64`-cut profile. Hence it cannot equal
`1` or `-1` for large `w`. Thus all five forms in (2) are nonzero,
and (1) gives six distinct rational output directions.

This does not give a quantitative real lower bound for the forms. In the
ordered chart `1<a<b<c`, the old matching coordinates are positive,
so only `F` among these forms can have real cancellation.

For the standard units in the order
`a,b,c,a-1,b-1,c-1,a-b,a-c,b-c`, the exponents of their additional
factors after (1) are

```text
F:  0  1 -1  1  0 -1  0 -1 -1
G:  0  0  0  0  0  0  1  1  1
H:  1  0  1  0  0  0  0  1  0
U:  0  0  0  0  1  1  0  0  1
V: -1 -1  0 -1 -1  0 -1 -1 -1.                           (3)
```

This integer matrix has rank four and a five-dimensional integer kernel.
An integral kernel basis gives the following exact transformations:

| Character | Transformed value |
|---|---|
| `bc/a` | unchanged |
| `(a-1)/b` | negative of the old value |
| `b(c-1)/(b-1)` | negative of the old value |
| `b(a-c)/(a(a-b))` | negative of the old value |
| `b(b-c)/((b-1)(a-b))` | unchanged |

These five rational functions are multiplicatively independent, as follows
from their exponent vectors in the nine distinct polynomial factors.
Algebraic independence is not asserted. Every integer combination has
exactly the same rational height before and after the map.

In particular,

```text
X=a(c-1)/(c(b-1)),       X'=-X,
X=Delta_42 Delta_63 Delta_15 /
                 (Delta_14 Delta_62 Delta_53).            (4)
```

The bracket word has zero signed degree at each row. Its internal cut
sum is `1` on twelve of the `64` subsets, `-1` on twelve, and zero on
forty. The full-cut height formula therefore gives

```text
h(X)=12w+o(w).                                            (5)
```

If the output also admits a full nearuniform, negligible-residue profile
of block scale `w'`, equations (4)--(5) imply `w'/w -> 1`. This excludes
fixed-power scale reduction **under those additional output-profile
hypotheses**, without assuming that the map preserves them.

## 2. A profile-independent radius bound

Let six distinct Gaussian lattice points lie on a circle of radius `R`
in an arc of length at most `C sqrt(R)`. For each perfect matching form
the Gaussian integer `P_M=product_(ij in M)(z_i-z_j)`.
Writing `z_i=z_center+R exp(i theta_i)` gives

```text
z_i-z_j=2iR exp(i(theta_i+theta_j)/2)
                         sin((theta_i-theta_j)/2).
```

Every matching uses each angle once. Thus the matching products have
one common complex phase up to sign, and their ratios belong to
`Q(i) intersect R=Q`. Write `P_M=gamma x_M`, where the `x_M` are
integers with gcd one. An integer Bezout combination gives
`gamma in Z[i]`, so `|gamma|>=1`. Consequently

```text
h([x_M])<=log max_M|P_M|
          <=3log(C sqrt(R))=(3/2)log R+3log C.             (6)
```

The matching point is unchanged by real projective changes of circle
parameters or individual row rescalings. Therefore (6) applies to any
circle realization of the output moduli point. It needs no common unit
class, primitive-row assumption, or restriction on the circle center.

## 3. Actual output height, including exceptional contents

Clearing common factors in (1) gives the raw matching map

```text
[-A F U G : -B V U G : C V U G : -D F V G : -E F^3].      (7)
```

Let `m` and `M` be the least and greatest absolute values among all
fifteen old matching monomials. Every pair of exceptional forms has a
sum or difference equal to twice one old matching monomial. Thus at most
one form has modulus smaller than `m`. If it is `F`, use the second
or third coordinate of (7); if it is `U,V`, or `G`, use the last.
The form `H` is absent from (7). Together with (2), this proves

```text
m^4 <= max|raw output| <=8 M^4.                           (8)
```

Use the actual central bracket factorization
`|Delta_e|=c_e b_e`, with

```text
c_e=product_(S containing both endpoints of e) n_S,
b_e=N(L_e)|t_e| in Z_(>0),
log b_e<=beta=2s+log T,
(1-eta)w<=log n_S<=(1+eta)w.
```

In particular every old matching has logarithmic size
between `48(1-eta)w` and `48(1+eta)w+3beta`.

The exact finite valuation certificate is
[check_joint_segre_local_content.py](check_joint_segre_local_content.py),
with its full table in
[joint_segre_local_content_audit.json](joint_segre_local_content_audit.json).
It uses all five individual decompositions (2) and all ten pair identities.
At a core prime, write `a_j` for the five old matching orders and `z_i`
for the orders of `F,H,U,V,G`. Each individual unequal-branch
decomposition fixes the smaller order. Each pair identity with matching
order `q` requires `min(z_i,z_j)<=q`, with equality when the two orders
differ. These are necessary valuation conditions, not assumptions of
generic position.

For all `64` cuts, the sum of the generic common orders in (7) is `144`.
Maximizing over the necessary conditions instead gives sum `150`:
there is one additional possible order at cuts `34`, `25`, `16`
and their complements. For example, at `a=1` one has
`A=E=0`, `V=H`, and `U=G=bc-1`. The extra relation `bc=1`
can occur while the remaining old directions are distinct. It is thus
incorrect to use `144w+o(w)` as a universal common-gcd upper bound.

Here is an explicit transfer from the finite table to actual integer
valuations. At a rational prime `p`, put

```text
B_p=sum_(all15 edges) v_p(b_e)+v_p(2).
```

If the core exponent at this prime is `e>0` and `B_p<e`, set

```text
u_i=ceil((v_p(form_i)-B_p)/e).                            (9)
```

The lower bounds and individual unequal-branch constraints pass to `u_i`.
For a pair with unequal `u_i`, the actual form valuations are unequal,
so their minimum is `q e+delta`, with `0<=delta<=B_p`.
Its shifted ceiling in (9) is exactly `q`. Equal `u_i` satisfy the
required upper bound as well. Therefore `u_i` is a feasible integer
word in the finite certificate. Since each coordinate of (7) uses one
old matching and three form factors, its actual common-gcd exponent
is at most `e k_S+4B_p`, where `k_S` is the table maximum.

If `B_p>=e`, pair identities show that at most one form has valuation
above `3e+B_p`. A coordinate omitting it has valuation at most
`12e+4B_p<=16B_p`. At primes outside the core the same argument gives
`4B_p`. Summing these bounds proves

```text
log gcd(raw output)<=150(1+eta)w+240beta+16log2.           (10)
```

The generic divisor still supplies the lower bound
`log gcd(raw output)>=144(1-eta)w`. Combining with (8) yields

```text
42w-342eta w-240beta-16log2 <= h_out
                  <=48w+336eta w+12beta+log8.             (11)
```

In particular the actual output height is at least `42w-o(w)`.
Equation (6) then gives, for bounded normalized output arc constants,

```text
log R_out>=28w-o(w).
```

Relative to the inherited input radius `log R_in=32w+o(w)`, this
projective-height estimate alone excludes radius descents stronger than
exponent `7/8`. The remaining interval cannot be removed by replacing
the six exceptional orders with their generic values. The actual common
Gaussian chord content supplies a stronger obstruction, as recorded below.

## Audit and scope

The stronger, actual Gaussian-circle statement for the operation retaining
the old anchor is proved in
[joint_involution_five_cut_inflation.md](joint_involution_five_cut_inflation.md):
five explicit cuts force `log|gamma|>=6.5w-o(w)`, and hence
`R_new>=R_old^(193/192-o(1))`. Its normalized arc constant diverges.
That argument uses common chord content in addition to the projective
height bound above; removing the original anchor changes its local width.
The two continuous standalone-cut estimates are recorded separately in
[joint_involution_two_pair_cuts.md](joint_involution_two_pair_cuts.md).

The completed standalone argument is
[joint_involution_standalone_inflation.md](joint_involution_standalone_inflation.md).
Sixteen directly controlled cuts force `log|gamma|>=8w-o(w)` without
retaining the old anchor. Consequently the arc-preserving ordering has
`R_out>=R_in^(49/48-o(1))`, and every labeling of this fixed-three pullback
has normalized arc length at least `exp(w/3-o(w))`. The latter assertion
does not require the transformed points to stay in the old arc. These
statements allow common rotations, but not arbitrary subsequent projective
changes of the output directions.

The root and fresh-algebraic agents independently checked the map formulas.
The symmetry agent supplied the exact local-content enumeration, including
the six exceptional cases. The uniformity-audit agent checked the integer
kernel, all `64` cut counts for (5), the phase and Bezout proof of (6),
and the shifted-ceiling error transfer and finite constants (10)--(11).
The calculations retain actual residues and prime powers. They do not
assert that all six exceptional contents occur globally at once.
