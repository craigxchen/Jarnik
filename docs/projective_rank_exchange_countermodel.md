# Rank, height gaps, and additive exchanges do not determine the balanced cut norm

## Status

This note gives an actual rational `M_0,k` family for every fixed `k`.
It has full multiplicative rank, a uniform height gap for every
integer projective-unit word, comparable growing heights for every
cross-ratio, and all the additive Plucker exchanges. Thus these
properties alone cannot bound `k`.

The family also matches the **entire** balanced-cut projective-unit
height norm on four and five points. Six points supply the first
distinction: a product of two explicit cross-ratios must lose an
additional `W/32+o(W)` of height in the balanced circle profile. This
isolates information a successful additive argument must retain.

The construction is not an endpoint lattice cluster and does not
have the full balanced Gaussian cut profile. It is a countermodel
only to arguments using the rank/height-gap conclusion, individual
cross-ratio heights, and universal additive identities without that
additional profile or metric information.

## 1. Distinct linear collision factors

Fix `k>=4`, set `q=k+1`, and for a positive integer parameter `T` put

```text
x_i(T)=i+q^i/T,       0<=i<k,
L_ij(T)=(j-i)T+q^j-q^i,       i<j.
```

The points are distinct, increasing rational numbers and

```text
x_j(T)-x_i(T)=L_ij(T)/T.                        (1)
```

All linear polynomials `L_ij` have distinct roots. Equivalently the
secant slopes

```text
alpha_ij=(q^j-q^i)/(j-i)
```

are distinct. For a fixed upper index `j`, strict convexity of `q^i`
makes these slopes strictly increase with the lower index. Slopes
with upper index `j` are all larger than every slope with upper index
less than `j`: their minimum is at least

```text
(q^j-q^(j-1))/(k-1)>q^(j-1),
```

whereas the latter slopes are less than `q^(j-1)`.

For different edges `e,f`, let `R_ef` be the nonzero integer resultant
of their linear polynomials. Evaluating the elementary linear
elimination identity gives the exact integer divisibility

```text
gcd(L_e(T),L_f(T)) divides |R_ef|.               (2)
```

Every pair gcd is therefore bounded independently of `T`.

## 2. Uniform height control for unbounded integer exponents

Let `c` be any integer edge vector with zero vertex degrees, as in
[coupled_crossratio_height_norm.md](coupled_crossratio_height_norm.md).
Then `sum_e c_e=0`, so its projective unit is

```text
Lambda_c(T)=product_e L_e(T)^c_e.                (3)
```

Write `n=||c||_1/2`, and form the positive integers

```text
A=product_(c_e>0)L_e(T)^c_e,
B=product_(c_e<0)L_e(T)^(-c_e).
```

Both unreduced products have `n` linear factors, counted with
multiplicity. If

```text
D_k=product_(e<f)|R_ef|,
B_k=max_(i<j)((j-i)+(q^j-q^i)),
```

then (2) implies

```text
gcd(A,B) <= D_k^(||c||_infinity).                (4)
```

For completeness, at each rational prime the common valuation of
the two products is at most the sum, over a positive edge and a
negative edge, of the minimum of their two weighted valuations.
Each such minimum is at most `||c||_infinity` times the unweighted
minimum, which is bounded by the valuation of `R_ef`. This retains
all prime powers and is uniform in the exponents `c_e`.

For `T>=1`, every factor satisfies `T<=L_e(T)<=B_k T`. Reducing
`A/B` therefore proves

```text
n log T-||c||_infinity log D_k
 <= h(Lambda_c(T))
 <= n log T+n log B_k.                          (5)
```

In particular, uniformly over **all** integer zero-degree vectors,

```text
h(Lambda_c(T))=(||c||_1/2)log T+O_k(||c||_1).    (6)
```

For example, once `log T>4 log D_k`, every nonzero vector satisfies

```text
h(Lambda_c(T)) >= (||c||_1/4)log T>0.            (7)
```

The evaluated projective-unit group consequently has rank exactly
`k(k-3)/2`. Equation (7) also supplies a Euclidean height gap at
least as strong as the form used in the power-size rank theorem,
after increasing `T` as needed. No relation with very large integer
coefficients can escape (5).

Every ordered quadruple cross-ratio has four nonzero coefficients
of absolute value one. Thus

```text
h(lambda_ijkl(T))=2 log T+O_k(1)                 (8)
```

simultaneously for all such quadruples. Their real values converge
to the finite, nondegenerate cross-ratios of `0,1,...,k-1`.
All Plucker and cluster exchange identities hold exactly because
the `x_i(T)` are actual rational points on a projective line.

## 3. Why four- and five-point unit heights do not distinguish it

For comparison set the formal height scale `W=16 log T`, so (8)
becomes `h(lambda)=W/8+O_k(1)`. This is a comparison scale, not a
claim about a lattice-circle radius of the constructed points.

For the balanced-cut norm use

```text
C_S(c)=sum_(i<j, i,j in S)c_ij,
ell(c)=(1/2)sum_unoriented_cuts |C_S(c)|,
w=W/2^(k-1).
```

Empty and singleton cuts contribute zero. On four rows the three
nontrivial cuts are the three pair matchings; opposite edge
coefficients agree by the zero vertex degrees. Hence

```text
ell(c)=||c||_1/4              when k=4.
```

Since then `w=2 log T`, the balanced prediction `w ell(c)` is
exactly the leading value in (6).

On five rows every nontrivial unoriented cut has a unique two-row
side, so its coefficient is that side's edge coefficient. Thus

```text
ell(c)=||c||_1/2              when k=5.
```

Now `w=log T`, giving the same agreement. Therefore full
multiplicative rank, every universal additive exchange, and the
entire leading projective-unit height norm on at most five rows
remain compatible with this explicit rational family.

## 4. The first six-row distinction is explicit

For six rows form

```text
u=(x_0-x_2)(x_3-x_4)/((x_0-x_4)(x_2-x_3)),
v=(x_0-x_5)(x_1-x_3)/((x_0-x_3)(x_1-x_5)).
```

Their product has the following nonzero edge coefficients:

```text
c_02=c_05=c_13=c_34=1,
c_03=c_04=c_15=c_23=-1.                        (9)
```

Thus `||c||_1=8`. The 15 cuts with a two-row side contribute
eight to the sum of absolute cut coefficients. The ten balanced
three-against-three cuts contribute six. Therefore

```text
ell(c)=(8+6)/2=7.
```

With `W=16 log T` and `w=W/32`, the balanced profile predicts

```text
h(uv)=7W/32+o(W).
```

But the actual rational family has, by (6),

```text
h(uv)=4 log T+O_k(1)=W/4+O_k(1).                (10)
```

The difference is `W/32`, one full balanced block weight.
Both individual cross-ratios have height `W/8+o(W)` in either
model, so this discrepancy cannot be detected by their individual
height estimates.

More precisely, the cut coefficients for `u` and `v` have opposite
signs on exactly one unoriented six-row cut:

```text
{1,2,3}|{0,4,5}:       C_S(u)=-1, C_S(v)=1.     (11)
```

The full balanced core therefore supplies a common cancellation of
logarithmic size `W/32+o(W)` when the two rational units are
multiplied. The generic linear-factor family has only bounded
pair-resultant cancellation, so this extra factor is absent.

## 5. Scope for an additive or S-unit continuation

The construction shows that a count argument cannot use only the
rank conclusion, the Euclidean height lower bound, the individual
cross-ratio heights, and their universal additive equations. These
are simultaneously realized here for every fixed number of rows.
It also shows that all five-point unit-height data miss the first
new balanced-cut interaction.

Equation (11) identifies a concrete six-point quantity to retain in
a further attempt: the large, specifically located cancellation
between two cross-ratios that each participate in many additive
exchanges. However, that cancellation is already supplied by the
corresponding core block. No incompatibility among all these
six-point cancellations has been proved.

In particular, this note does not bound the actual primitive
residues, and its rational family is not a small-residue realization
of the balanced Gaussian system. The remaining question is whether
the complete collection of such large shared factors can coexist
with the actual Gaussian phase and small-residue constraints.

## 6. The first nonlinear six-point exchange, with its denominator retained

Keep the units `u,v` from Section 4 and put `F=1-uv`. In the actual
balanced six-row system write the reduced product as `uv=p/q`, with
`q>0`. The conductor calculation has seven positive and seven negative
cut coefficients, all of absolute value one. The correction factor in
the simultaneous height theorem has height `o(w)`, where `w=W/32`.
It follows separately for the numerator and denominator that

```text
log |p|=7w+o(w),       log q=7w+o(w).                       (12)
```

This uses the positive and negative conductor products separately;
the height of their maximum alone would not suffice. Multiplication
by a rational factor of height `o(w)` changes each of these two
logarithms by at most `o(w)`.

The exact nonlinear exchange is

```text
F=(q-p)/q,       gcd(q-p,q)=1.                              (13)
```

The full unit-height lower bound implies `uv!=1`, so `q-p` is a
nonzero integer. Thus its presently available separation is precisely

```text
|1-uv| >= 1/q = exp(-7w+o(w)).                              (14)
```

The common block already cancelled in forming `uv` does not divide
the reduced numerator and denominator of (13) a second time. To turn
this particular exchange into a contradiction, an additional argument
would have to prove `|1-uv|<1/q`. Neither the shared cut nor the fact
that the six original angles occupy a short interval implies that
estimate: projective cross-ratios retain the relative positions inside
that interval.

One can also check exactly what the generic collision orders of this
non-unit function contribute. With the cut side excluding zero, its
seven pole cuts are

```text
23, 15, 125, 235, 1235, 145, 1245.
```

Each has order `-1`. At every other nontrivial cut the generic order
is zero, including the shared cancellation cut `123|045`. There are
no additional generic boundary zeros. This assertion concerns formal
collision orders, not a prohibition on extra prime divisibility at
the actual integer configuration. In particular the actual numerator
`q-p` could be much smaller than `exp(7w)`; no such smallness has been
proved here.

For an exact verification of the collision-order assertion, substitute
`x_i=101+epsilon*d_i` on the cut and `x_i=b_i` off the cut, using

```text
b=(2,7,13,23,37,53),       d=(3,11,19,29,43,61).
```

Expand the numerator

```text
(x0-x4)(x2-x3)(x0-x3)(x1-x5)
 -(x0-x2)(x3-x4)(x0-x5)(x1-x3)
```

and the first product as its denominator. On all 25 nontrivial cut
sides, their first nonzero integer coefficients give exactly the
orders above. The elementary edge orders first give the lower bound
`min(0,C_S(u)+C_S(v))`; a nonzero leading coefficient at this one
specialization certifies that no further generic cancellation occurs.
Thus the finite check certifies nonvanishing of the relevant formal
leading polynomials, rather than substituting numerical evidence for
an inequality about all arithmetic configurations.

The independent uniformity-audit agent checked Sections 1--5,
including the resultant bound for unbounded exponent vectors and the
complete four- and five-row norm comparisons.
