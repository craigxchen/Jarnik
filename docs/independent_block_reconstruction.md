# Reconstructing endpoint clusters from independent block equations

## Status

The central-truncation theorem extracts independent Gaussian blocks
and small primitive equations from a hypothetical failure of the
uniform endpoint bound. This note proves a converse under explicit
height and marginal-balance conditions. Such block data reconstruct
actual primitive lattice-circle configurations, with a quantitative
bound on their normalized arc length.

This does not prove that the block systems exist or that they are
impossible. It identifies a sufficient arithmetic target whose
exclusion would settle the original goal. In particular, the
reconstruction retains integrality and all denominator costs.

## 1. Data and the quantitative reconstruction bound

Let `m>=1`. For each subset `S` of `[m]`, let `H_S` be a nonzero
Gaussian integer. Assume all the blocks and their conjugates have
the prescribed coprimality: each block is coprime to its conjugate,
distinct blocks are coprime, and distinct blocks are coprime to
each other's conjugates. Units are allowed unless a size condition
below excludes them. Define

```text
A_i = product_(S contains i) H_S,
h_i = K_i A_i,       K_i in Z[i] minus {0},
gcd_G(h_i,bar(h_i))=1,
t_i=Im h_i != 0.
```

Set

```text
P = product_(S nonempty) |H_S|,
E = sum_(i=1..m) log|K_i|,
a = min_i log|A_i|,
T = max_i |t_i|,
sigma = a - (1/2)log P,
beta = T exp(-a).
```

Assume the phase ratios `h_i/bar(h_i)` are distinct, and
`beta<=1/2`. Then there are `m+1` distinct Gaussian integers on
a primitive circle of radius `R`, occupying an arc of length at
most `C_data sqrt(R)`, where

```text
P <= R <= P exp(E),
C_data <= 8 T exp(E/2-sigma).                    (1)
```

The phase-distinctness condition is automatic if `m>=2` and

```text
max_i log|K_i| < min_(S nonempty) log|H_S|.       (2)
```

For `m=1`, nonzero `t_1` already distinguishes the one ratio from
the anchor. Condition (2) is a convenient sufficient condition,
not an extra requirement when distinctness is known directly.

## 2. Exact integrality and the least possible radius

Let `L` be a Gaussian least common multiple of `h_1,...,h_m`.
Define

```text
z_0=bar(L),
z_i=bar(L) h_i/bar(h_i),       1<=i<=m.          (3)
```

These are Gaussian integers because `bar(h_i)` divides `bar(L)`.
They all have modulus `R=|L|`.

Every nonempty block occurs in at least one `A_i`, hence divides
at least one `h_i`. By block coprimality their product divides
`L`. Conversely, every `h_i` divides
`(product_(S nonempty) H_S)(product_j K_j)`, so `L` divides that
same Gaussian integer. Taking moduli gives the radius bounds in
(1). No prime support has been discarded in this comparison.

This is the **least possible radius** for these phase ratios with
the anchor ratio `1` included. Indeed, let an integral anchor be
`b!=0`. If `b h_i/bar(h_i)` is integral, conjugate coprimality
implies `bar(h_i)|b`. Thus `bar(L)|b`, so `|b|>=|L|`.

The configuration (3) is already primitive. To see this, consider
a Gaussian prime dividing `bar(L)`, and choose `h_i` having
maximal valuation at its conjugate prime. Since `h_i` is coprime
to its conjugate, the valuation of `z_i` at the original prime is
zero. A prime not dividing `bar(L)` is already absent from `z_0`.
Consequently the common Gaussian gcd of (3) is a unit.

## 3. Distinctness and the arc estimate

For primitive Gaussian integers `h_i,h_j`, equality of their
conjugate ratios implies `h_i=+h_j` or `h_i=-h_j`. One can see
this either from uniqueness of a reduced Gaussian fraction or
by observing that their quotient is rational real and that both
coordinate pairs have integer gcd one.

If `i!=j`, choose a subset containing `i` and excluding `j`.
Equality up to sign would imply `H_S|K_j A_j`; block coprimality
then gives `H_S|K_j`, contrary to (2). Each ratio also differs
from the anchor `1`, since `t_i!=0` and a primitive real Gaussian
integer is a unit. This proves the asserted distinctness criterion.

For each `i`,

```text
|Im h_i|/|h_i| <= T/|A_i| <= beta.
```

Changing the sign of `h_i` does not change its conjugate ratio.
Choose that sign to make its real part positive. The argument of
`h_i` then lies in `[-arcsin(beta),arcsin(beta)]`; the arguments
of all ratios, including the anchor, lie in an interval of length
at most `4 arcsin(beta)`. For `beta<=1/2`, this is at most
`8 beta`, less than a full turn. Multiplying by `z_0` preserves
the angular span, so the normalized arc length is at most

```text
8 beta sqrt(R)
 <= 8 T exp[-a+(log P+E)/2]
 = 8 T exp(E/2-sigma).
```

This proves (1). The smallness of the original correcting factors'
arguments was never assumed.

## 4. The role of the empty block and marginal balance

Write

```text
w_S=log N(H_S),
W_c=sum_S w_S.
```

Then, exactly,

```text
log P=(W_c-w_empty)/2,
log|A_i|=(1/2)sum_(S contains i) w_S,
sigma=min_i[(1/2)sum_(S contains i)w_S-W_c/4]
       +w_empty/4.                              (4)
```

Thus the empty block supplies positive slack once the row marginals
are close enough to `W_c/2`. In the formal exactly balanced weight
model `w_S=w`, (4) gives `sigma=w/4`. Exact equality of distinct
positive block norms is not an arithmetic construction: coprime
nonunit blocks cannot have equal norms. The calculation describes
the weight model only.

For actual sequences it suffices to have

```text
w_empty=w(1+o(1)),       w -> infinity,
min_i log|A_i| = W_c/4+o(w),
E+log T=o(w).                                   (5)
```

These hypotheses imply `sigma=w/4+o(w)`. They also imply
`beta<=1/2` eventually, and (1) gives

```text
C_data <= exp[-w/4+o(w)] -> 0.                  (6)
```

If the number of rows tends to infinity and their ratios are
distinct, such sequences would contradict the original uniform
bound, even for one arbitrarily small fixed positive arc constant.

Near-uniform individual block weights alone do not imply the
second condition in (5) when the number of blocks grows: marginal
errors can accumulate across many blocks. The separate pair-moment
estimate in central truncation is what supplies that condition.

## 5. Central truncation supplies the full converse hypotheses

Use the notation of
[endpoint_central_truncation.md](endpoint_central_truncation.md).
There are `k=m+1` selected original rows, original scale
`W=log R_original^2`, correction scale `s=log D`, and core total
weight `W_c=W-2s`. The core products and primitive factors satisfy

```text
h_i=K_i A'_i,
log|K_i|<=s,
W/4-s <= log|A'_i| <= (1+zeta)W/4,
log T <= zeta W/4,
s<=kW/sqrt(M),       zeta=O(k/sqrt(M)).
```

Let `w=w'_empty`. All core blocks have weight comparable to
`W/2^(k-1)`, and their relative errors tend to zero. Formula (4)
and the displayed row bounds give the finite estimates

```text
w/4-s/2 <= sigma <= w/4+s/2+zeta W/4,
E <= (k-1)s,
log C_data <= log 8 - w/4 + ks/2 + zeta W/4.     (7)
```

Take `k=floor(c log_2 M)` with any fixed `0<c<1/2`. Then

```text
(ks+zeta W)/w = O(k^2 2^k/sqrt(M)) -> 0.
```

For a fixed original constant `C<sqrt(2)`, the strict pair gap
also gives `W` bounded below by a positive constant times `M`.
Consequently `w` tends to infinity. All assumptions in (5), the
distinctness criterion (2), and (6) therefore follow from the
audited central output.

For data actually obtained from a circle, reconstruction simply
recovers the selected phase ratios at their least primitive radius.
Its additional value is a converse for abstract candidate data:
one need not assume they were originally produced by a circle,
nor silently overlook a denominator cost. Equations (1)--(7)
show exactly when the candidates themselves produce the requested
kind of actual counterconfiguration.

## 6. Scope and audit

This yields an equivalent way to seek a contradiction: exclude
unbounded families of independent, primitive block equations with
the marginal slack and total correction budget above. Central
truncation would produce them from any failure of uniformity;
reconstruction would turn any such family back into actual
endpoint clusters.

No existence or exclusion of these families is proved here. The
uniform bound remains open in this project.

An independent audit checked the Gaussian lcm construction, its
minimality and common-gcd normalization, the distinctness criterion,
the angular constants, and the growing-row height calculation.
Exact arithmetic checks covered 1,000 independent Gaussian block
constructions with varying powers and correcting factors. The audit
is recorded in Section 9 of
[uniformity_audit_parallel.md](uniformity_audit_parallel.md).
No Lean formalization is claimed.
