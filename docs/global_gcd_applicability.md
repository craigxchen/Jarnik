# Prime-mass escape and global GCD height inequalities

## Status

This note tests global GCD inequalities on the actual normalized
anchor equations. The normalization makes the shifts fixed, so
moving targets are not the obstacle for this application. The
outside-prime contribution is the obstacle: it has a leading
height much larger than the common divisor being estimated.

The calculation does not exclude a different rational map, a
different boundary divisor, or a theorem treating the entire
incidence arrangement simultaneously. It supplies no uniform
endpoint bound.

## 1. Every fixed prime set loses the conductor mass

Consider a hypothetical sequence of endpoint clusters with
`M -> infinity`, fixed `C<=1`, and `W=log R^2`. The single-prime
profile bound in [inverse_valuation_profile.md](inverse_valuation_profile.md)
gives, for every rational prime power `p^e || R^2`,

```text
e log p <= 4W/M.                                 (1)
```

Indeed that note has `m_e log p<=W/M` and `m_e>=e/4`.
For any fixed finite rational prime set `S_f`, therefore,

```text
sum_(p in S_f) e_p log p <= 4|S_f| W/M = o(W).    (2)
```

The conclusion also holds for any varying set with cardinality
`o(M)`. For a primitive anchored numerator `h_i`, each rational
norm valuation is at most the corresponding `e_p`. Its norm
mass on a fixed prime set is consequently `o(W)` as well.

This is stronger than failure to extract a fixed supporting set:
every proposed fixed set carries an asymptotically negligible
fraction of the conductor. It uses the original cardinality
`M`, not merely the number of extracted rows.

## 2. Fixed shifts from the small-value equations

Fix `n>=2` nonanchor positions in the central extraction, and
use the equations

```text
h_i=K_i A_i=x_i+i t_i,   0!=t_i in Z,
log|h_i|=W/4+o(W),
sum_i log|K_i|+sum_i log(2|t_i|)=o(W).            (3)
```

The fixed-`n` projection of the full cut profile is uniform.
Thus a block belonging to all `n` selected rows has total norm
weight `W/2^n+o(W)`. Block coprimality and the small correcting
factors give

```text
log|gcd_G(h_1,...,h_n)|=W/2^(n+1)+o(W).          (4)
```

More explicitly, the common core divisor is the product of the
blocks containing all selected indices. It divides the gcd.
Any additional common valuation is bounded in total by the sum
of the logarithmic heights of the `K_i`, so the error in (4)
is `o(W)`.

Define a point over `Q(i)` by

```text
y_i=h_i/(2i t_i),    y_i-1=bar(h_i)/(2i t_i),
Q=[1:y_1:...:y_n].                              (5)
```

All target shifts are now exactly one. No theorem for moving
shifts is needed to formulate a GCD estimate at this point.
The denominators in (5) have total logarithmic height `o(W)`.
At the complex place, `log max(1,|y_i|)=W/4+o(W)`, and the
total finite-place height is at most the denominator height.
With absolute Weil heights normalized over `Q(i)`, this proves

```text
h(Q)=W/4+o(W).                                  (6)
```

## 3. Exact leading size of the boundary term

Let `S` be any fixed finite set of places of `Q(i)` containing
the complex place. Put `B=(X_0 X_1 ... X_n=0)` and

```text
B_out(Q)=sum_(v outside S) lambda_v(B,Q).
```

For the affine coordinates in (5), set
`m_v=log max(1,|y_1|_v,...,|y_n|_v)`. The standard local
height formula is

```text
lambda_v(B,Q)=(n+1)m_v-sum_i log|y_i|_v.         (7)
```

The sum of the finite `m_v` is `o(W)`. Moreover, (1)--(3)
show that the numerator and denominator mass at the fixed
finite places in `S` is `o(W)`. The product formula and (3)
then give

```text
B_out(Q)=sum_i log|y_i|+o(W)
        =nW/4+o(W).                             (8)
```

This handles denominator places and the normalization at split
primes explicitly. At a Gaussian prime of rational norm `p`,
an exponent contributes `(log p)/2`; the complex place uses
the usual modulus. Thus (4), (6), and (8) use the same absolute
height normalization.

The global exceptional-divisor height at the point
`[1:1:...:1]`, written `log GCD^+(y_1-1,...,y_n-1)`, satisfies

```text
log GCD^+(y_1-1,...,y_n-1)
   =W/2^(n+1)+o(W).                              (9)
```

At finite places away from the small denominators it is exactly
the log gcd in (4), with conjugation. The total change at the
denominators is `o(W)`. Its complex contribution is exactly zero
eventually: `y_i-1=-bar(y_i)`, and `max_i |y_i|>=1` for large
height. This proves (9).

## 4. Comparison with the available global inequalities

Yasufuku's [2024 paper, Theorem 5(ii)](https://link.springer.com/article/10.1007/s00605-024-02018-1)
gives, for fixed shifts over a number field and outside a proper
hypersurface, an upper bound of the form

```text
log GCD^+ < epsilon_2(n) h(Q) + B_out(Q)/(n-1),  (10)
```

where `epsilon_2(n)>0` can tend to zero as `n` grows. Applying
our calculations, its two sides have leading coefficients

```text
left:  1/2^(n+1),
right: epsilon_2(n)/4 + n/(4(n-1)).              (11)
```

Even replacing `epsilon_2(n)` by zero leaves the right coefficient
strictly larger for every `n>=2`. Thus (10) is compatible with
the hypothesized profiles even away from its exceptional set.
Controlling that exceptional set alone would not fix this
particular application.

Yasufuku's [2026 codimension-two result, Theorem 1.2](https://londmathsoc.onlinelibrary.wiley.com/doi/full/10.1112/blms.70349)
also has an outside-`S` boundary term, with coefficient one,
in addition to a degree-dependent height term. Its center must
be a codimension-two complete intersection avoiding every
coordinate point. These hypotheses and the boundary term must
both be retained when choosing two forms. Merely replacing an
`S`-unit theorem by this all-points theorem does not remove the
loss in (8).

The source statements in this section are distinct from the
original height calculation in Sections 1--3. No general
impossibility theorem for global GCD methods is asserted.
The untested possibility is a transformation or simultaneous
arrangement inequality whose useful divisor mass is large
relative to its correctly computed boundary cost.

## Verification boundary

The proof uses the previously audited single-prime profile and
central-truncation results, elementary Gaussian gcd valuations,
the product formula, and the displayed definitions of local
heights. An independent audit checked the normalization, both
denominator errors, the extra gcd valuation bound, and the two
original source statements. These calculations and
source-applicability checks are prose; no new Lean theorem is
claimed.
