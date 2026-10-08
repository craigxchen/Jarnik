# Conductor-adaptive CRT and the exact half-depth threshold

## Status

The finite-depth construction in
[power_class_saturation.md](power_class_saturation.md) includes a complete
proof for prescribed depths at every chosen core prime. Its correction
factors are fixed before the core exponents become large. It does not
treat congruence depths comparable to those exponents.

This note makes that boundary precise. Congruences covering slightly
more than half of a row's core norm already force the actual small-value
equation, provided coefficients and targets are small. At exactly half
depth, the implication can fail even for primitive Gaussian integers
and constant coefficients and targets. The same congruences also give
uniqueness of the small correction-to-target ratio.

These are exact arithmetic statements. No cross-row height lower bound
or uniform endpoint bound is proved here.

## 1. Branch congruences are ordinary integer congruences of larger modulus

Let `A` be a nonzero Gaussian integer coprime to its conjugate and
supported on odd split primes. Let `L|A`, and take `K in Z[i]` and
`t in Z`. Then

```text
bar(K) bar(A) = -2i t mod L

iff

N(L) divides Im(KA)-t in Z.                     (1)
```

Indeed, `KA=0 mod L`, so the first congruence is equivalent to

```text
2i(Im(KA)-t)=0 mod L.
```

Since `L` is coprime to `2i`, it divides the real integer
`Im(KA)-t`. Its conjugate also divides that integer. The two
Gaussian divisors are coprime, so their product `N(L)` divides
it. Conversely, integer divisibility by `N(L)` gives the displayed
Gaussian congruence. This proves both directions without discarding
any prime powers.

For independent core blocks, choose any divisors `J_S|H_S` and set

```text
A_i=product_(S contains i) H_S,
L_i=product_(S contains i) J_S.
```

The branch congruences modulo every incident `J_S` are equivalent
to the single congruence modulo `L_i`, since these factors are
pairwise coprime. Equation (1) then supplies the rational integer
modulus `N(L_i)`, rather than just its square root.

## 2. Slightly more than half depth forces the exact integer value

Under the preceding assumptions, if the branch congruence holds and

```text
N(L)>|K| |A|+|t|,                              (2)
```

then `Im(KA)=t` exactly. The difference in (1) is an integer of
absolute value strictly smaller than its nonzero modulus.

To express this as a depth threshold, write

```text
theta=log N(L)/log N(A).
```

If `|A|` tends to infinity and

```text
log max(1,|K|,|t|)=o(log|A|),
```

then every fixed `theta>1/2` satisfies (2) eventually. More
quantitatively, if `|t|<=|K||A|`, it suffices to have

```text
log N(L) > (1/2)log N(A)+log|K|+log 2.          (3)
```

Thus it is unnecessary to impose every prime power of the full
row conductor. Enough branch congruences to cover more than half
its logarithmic norm already recover the actual small-value equation.
The complete full-depth CRT problem, with the current negligible
coefficient and target heights, is therefore a faithful reformulation
of that equation, not a weaker local substitute.

## 3. An exact primitive example at the half-depth boundary

Let `u` be a positive even integer, and set

```text
L=u+i,
A=(u+i)^2,
K=i,
t=-2.
```

Then

```text
Im(KA)-t=u^2-1+2=u^2+1=N(L),
N(A)=(u^2+1)^2.
```

Consequently the branch congruence in (1) holds, but the actual
imaginary part is of order `|A|`, not bounded. Meanwhile

```text
log N(L)/log N(A) = 1/2.
```

All relevant primitivity assumptions hold. The coordinates of `L`
are coprime and have opposite parity, so `L` is coprime to its
conjugate. The same is true of its square `A` and of `KA`.
The unit correction `K=i` has logarithmic height zero. Repeated
prime powers are permitted throughout the model.

This proves that the strict half-depth threshold cannot be lowered
by the modulus-versus-size argument alone. It is a one-row example;
it is not asserted to have a growing uniform cut profile.

## 4. A single large incident block determines a small target ratio

Keep `A,L` fixed. Suppose two pairs `(K,t)` and `(K',t')` satisfy
the branch congruence modulo `L`. Multiplying the two congruences
by the other integer target and subtracting gives

```text
L divides t bar(K')-t' bar(K),
bar(L) divides t K'-t' K.                       (4)
```

We used that `bar(A)` is invertible modulo `L`, which follows
from `L|A` and `gcd(A,bar(A))=1`.

Therefore, if

```text
|t| |K'|+|t'| |K| < |L|,
```

the integer Gaussian numerator in (4) vanishes. In particular,
among solutions with nonzero targets and bounds

```text
|K|<=U,        0<|t|<=T,        2UT<|L|,
```

there is at most one rational Gaussian ratio `K/t`.

At the extracted scale, coefficients and targets have logarithmic
heights negligible relative to **each individual** core block.
Thus even a single incident block `H_S` determines the small
ratio `K_i/t_i` uniquely, once `A_i` is fixed. This is a genuine
rational-reconstruction uniqueness statement; it does not assert
that such a small ratio exists.

## 5. The concrete CRT counting attempt and its limitation

Write `A=x+iy`, with `gcd(x,y)=1`, and set `n=N(L)`.
The branch problem is exactly

```text
yp+xq-t=0 mod n,        K=p+iq.                 (5)
```

For a fixed target `t`, its solutions form a translate of a
rank-two lattice of index `n` in `Z^2`. Surjectivity of
`(p,q) -> yp+xq mod n` follows from `gcd(x,y)=1`.
If the target is included as a third integer variable, the kernel
in `Z^3` likewise has index `n`, since the coefficient of `t`
is minus one.

The large index by itself cannot exclude a very short solution.
For example, take any positive even integer `x`, set `A=x+i`
and `L=A`. Then the modulus is `x^2+1`, while the fixed vector
`(p,q,t)=(1,0,1)` solves (5), and indeed solves the exact integer
equation. Thus comparing the volume of a coefficient box with
the CRT index is not a valid deterministic nonexistence proof.
The uniqueness in Section 4 is compatible with exactly one such
short rational target.

For a common block `H_S`, the congruences for two incident rows
give

```text
H_S divides
 t_j bar(K_i)bar(A_i)-t_i bar(K_j)bar(A_j).      (6)
```

Once the actual small-value equations hold, the expression in
(6) is real, so `N(H_S)` divides it. This holds simultaneously
for all common core blocks. For primitive `h_i`, the exact gcd
normalization additionally accounts for the correcting factors
and gives the normalized determinant relation

```text
N(gcd_G(h_i,h_j)) t_ij=x_i t_j-x_j t_i,
h_i=x_i+i t_i,
```

with its correcting factors as described in
[all_anchor_primitive_transition.md](all_anchor_primitive_transition.md).
Eliminating real coordinates from several such relations gives
the known rank-two cocycles. This direct CRT elimination has not
produced a new height inequality beyond those relations.

The unresolved arithmetic question is therefore specific: can
the row lattices (5), whose coefficients come from independent
shared Gaussian blocks, simultaneously have these unique tiny
target ratios? Their individual indices, the uniqueness lemma,
and prescribed finite-depth congruences do not answer that
question. A successful conductor-adaptive CRT argument needs
additional control of the simultaneous shortest targets; none
is supplied by the present index calculation.

An independent audit checked the exact norm modulus, the strict
threshold and its primitive boundary example, the ratio uniqueness,
both lattice indices, and the point at which the cross-row expression
becomes real. These are prose proofs; no Lean formalization is claimed.
