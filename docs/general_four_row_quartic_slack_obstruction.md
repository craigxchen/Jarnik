# General four-row slack has an unavoidable quartic defect

The Walsh four-row argument has the following exact extension to arbitrary
sign rows. Its error term need not be small even after removing rows:
an unbounded Paley family has error `(1-o(1))W` at **every** quadruple.
This obstructs extraction of small-defect quadruples from general
Hadamard profiles. It is not an endpoint configuration or a counterexample
to a possible general maximum-slack theorem. The four-row argument here
alone does not treat all characters; the subsequent
[prime-field theorem](paley_primefield_polynomial_character_gap.md)
does so at unbounded orders for capacity-two assignments.
Indeed the constructed profiles themselves have `u_max=O(W/M)`:
each off-diagonal Gram differs from `-b tau` by at most `4tau+1`.

## Exact identity for arbitrary weighted sign rows

Let physical columns have positive weights `w_j`, and put

```text
W = sum_j w_j,
G_xy = sum_j w_j s_xj s_yj,
T_1234 = sum_j w_j s_1j s_2j s_3j s_4j,
kappa = 4log C-log 4-D,             u_xy = kappa-G_xy.
```

For four distinct rows take `c=(1,1,-1,-1)` and its exact character height
`V_c = sum_j (w_j/2)|s_1j+s_2j-s_3j-s_4j|`. Directly on the sixteen
sign patterns,

```text
|z1+z2+z3+z4|/2 = 3/4 +(1/4)sum_(i<k) zi zk -(1/4)z1z2z3z4.
```

Substitute `z=(s1,s2,-s3,-s4)` to obtain

```text
V_c-W/2 = (W-T_1234)/4-kappa/2
          +(u13+u14+u23+u24-u12-u34)/4.                 (1)
```

Whenever this character is nonunit and the actual source phases satisfy
the endpoint arc hypothesis, its usual support-four lower bound is
`V_c >= (W+D)/2+log 8-2log(4C)`. Thus

```text
u12+u34 <= u13+u14+u23+u24 +(W-T_1234)+8log 2.          (2)
```

The three balanced pairings obey corresponding versions of (2), all with
the same nonnegative defect

```text
W-T_1234 = 2 sum_(j: s1j s2j s3j s4j=-1) w_j.          (3)
```

In a Walsh parallelogram the unflipped defect vanishes. Averaging its
flip errors is precisely what permits the stronger Walsh estimate.
For arbitrary sign profiles no such closure is supplied by (1).

## Unbounded Paley cores have no small-defect quadruples

Let `q` be a prime congruent to three modulo four, `M=q+1`, and `chi`
the quadratic character, with `chi(0)=0`. The normalized Paley matrix is

```text
H_(infinity,j)=H_(x,infinity)=1,
H_(x,j)=-chi(x-j) for x!=j,       H_(x,x)=-1.
```

It has `HH^T=MI`. For four distinct rows `Q`, let
`t_Q=sum_(j in F_q) product_(x in Q) H_(x,j)`, omitting the constant
column. If all four rows are finite, the summand away from their four
diagonal positions is the quadratic character of a squarefree quartic.
If one row is infinity, it is a constant sign times the character of a
squarefree cubic, away from three positions. Consequently

```text
|t_Q| <= 3sqrt(q)+4                                   (4)
```

uniformly for **every** quadruple. This uses Weil's multiplicative
character estimate `(d-1)sqrt(q)` for a squarefree degree-`d` polynomial;
the diagonal replacements cost at most `d`. A primary research source
stating the needed estimate is Kim--Yip--Yoo,
[Explicit constructions of Diophantine tuples over finite fields,
Lemma 2.1](https://link.springer.com/article/10.1007/s11139-024-00888-5).

Take a fixed `b>=5` copies of each nonconstant column. Assign one distinct
physical column to each row and flip its entry there. Such assignments
exist because `bq>=q+1`; no label-multiplicity bound is needed. Initially
give all physical columns weight `tau`. A quadruple encounters exactly
four assigned physical flips, so its quartic correlation changes by at
most `8tau`. Thus

```text
W-T_Q >= [b(q-3sqrt(q)-4)-8]tau = (1-o(1))W.            (5)
```

This holds after *any* row deletion leaving at least four rows, and under
any averaging with nonnegative weights over retained quadruples. In
particular it precludes a uniform `O(W/M+1)` defect bound by row removal
or choosing special quadruples in general Hadamard matrices.

## Actual prime weights and the scope of the obstruction

The nearby-prime clustering argument in
[the existing construction](strict_obtuse_prime_box_countermodels.md)
gives `r=bq` distinct rational primes congruent to one modulo four with

```text
tau <= log p_j < tau+1/r,        tau=4log M+O(1).
```

For completeness, partition `[M^4,2M^4]` into `r` equal intervals. The
fixed-modulus prime count gives at least `r` split primes in some interval
for large `M`; use the log of its left endpoint as `tau`. The log-width
is less than `1/r`. Hence `W=O(Mlog M)`. Every correlation changes by
less than one under this perturbation, and (5) remains valid with an
additional harmless absolute error of two.

Moreover, every pair has nominal Gram at most `(-b+4)tau`; actual Gram
is at most `(-b+4)tau+1`. Therefore all pair phase inequalities hold for
`C=1,D=0` and sufficiently large `q`.

Every balanced four-row certificate also satisfies the phase-height
lower bound. Before flips, its pair correlations are all `-b tau`, so
(1), or its Gram form, gives

```text
V_c-W/2 = (W-T_Q)/4+b tau/2.
```

Four flips lower `V_c` by at most `4tau`; actual prime perturbations
alter `V_c-W/2` by less than `5/2`. Combining this with (4), all three
balanced pairings have

```text
V_c-W/2 >= b(q-3sqrt(q)-4)tau/4+b tau/2-4tau-5/2 > 0
```

for large `q`. This is stronger than the support-four necessary bound
`V_c-W/2 >= -log 2` at `C=1,D=0`. These characters are nonunit because
their heights are positive. Arbitrary orientations of Gaussian factors
give literal primitive sign-profile rows with those exact heights, but
their angular span is not controlled here.

Thus even simultaneous pair and balanced four-row necessary phase-height
inequalities do not force small quartic defects. Their satisfaction does
not imply simultaneous actual source-phase localization. The subsequent
[prime-field theorem](paley_primefield_polynomial_character_gap.md)
extends this obstruction to every individual primitive half-height
inequality at unbounded orders. Different arithmetic arguments and
joint source-phase constraints remain unresolved.

The companion checker verifies the sixteen-sign identity, the exact
weighted identity (1), Hadamard orthogonality, and every quadruple at
orders `12,20,32,44`. It checks the character-sum comparison including
diagonal errors and the one-flip perturbation bound. The asymptotic
claim (4) relies on the stated Weil bound, not finite enumeration.
