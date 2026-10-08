# A prime-field polynomial gap for Paley characters

For the normalized Paley Hadamard matrix of prime order parameter
`q=3 mod 4`, put `M=q+1`. Every nonzero zero-sum ternary row vector `c`
of support `h<q` satisfies

```text
s := #{a:(H^T c)_a != 0} >= M/2+1-h,
L := ||H^T c||_1/M >= 2h(M+1-h)/(M(h+2)).                 (1)
```

At square-root support the second bound tends to two. This closes the
support interval left between the [fourth-moment bound](paley_sqrt_support_fourth_moment_gap.md)
and the [high-support capacity bound](hadamard_ternary_support_spectral_gap.md).
It uses reduction modulo the prime `q`, rather than a character-sum
estimate. In particular, for `M>=36`, five or more repeated columns and
**every** one-flip assignment with label capacity two have all nonzero
saturated characters above primitive half-height. Thus the individual
primitive half-height inequalities alone do not exclude these profiles.
This is not a short-arc realization or an endpoint theorem.

## The polynomial and parity argument

Index rows and columns by infinity and `F_q`. Write
`Q_xy=chi(x-y)`, so the finite block of `H` is `-I-Q`, with the infinity
row and column all ones. Thus `Q^T=-Q` and `Q^2=J-qI`.
Set `epsilon=c_infinity`, write `f` for the finite part of `c`, and put
`T=H^T c`. Since `sum c=0`,

```text
T_infinity=0,    sum f=-epsilon,
T_finite=epsilon*1-f+Qf.
```

Euler's criterion shows that the integer vector `Qf+epsilon*1`, reduced
modulo `q`, is the evaluation on `F_q` of

```text
P(X)=epsilon+sum_y f_y (X-y)^((q-1)/2) in F_q[X].          (2)
```

If `epsilon=0`, its degree is at most `(q-3)/2`; otherwise its degree is
at most `(q-1)/2`. This polynomial is nonzero. Indeed each integer
entry of `Qf+epsilon*1` has absolute value at most `h<q`, so an identically
zero reduction would make that integer vector zero. For nonzero epsilon
this contradicts its sum `q*epsilon`. For zero epsilon it gives `Qf=0`,
and `Q^2f=-qf` then contradicts `f!=0`.

A nonzero polynomial has at most its degree many roots. Since
`Qf+epsilon*1=T_finite+f`, its support is contained in the union of the
supports of `T` and `f`. The two cases give, respectively,

```text
s+h >= M/2+1,       s+(h-1) >= M/2.
```

These are the same first inequality in (1). Every nonzero transform
entry has even absolute value `2<=t<=h`: parity follows from `sum c=0`.
For such an entry `(t-2)(h-t)>=0`, whence

```text
Mh=sum T_a^2 <= (h+2) sum |T_a|-2hs.
```

Inserting the support bound proves the second inequality in (1).

## Equality characters at norm one

For prime Paley `M>12`, every zero-sum ternary `c` with `L=1` is a
row-coordinate pair difference, a half-sum or half-difference of two
distinct nonconstant columns, or a signed nonconstant column.

To prove this, equality in Hadamard inversion at each support row forces
every active transform entry to have magnitude `h`. Thus `T=h v` for a
sign vector on `s` coordinates, with `hs=M`. If `h=M`, inversion gives
a signed column. Otherwise `h<=M/2<q`, and (1) gives
`h+s>=M/2+1`. If both `h,s>=3`, then `hs=M` implies
`h+s<=M/3+3<M/2+1`, a contradiction. The remaining possibilities are
`h=2`, which gives a row pair, or `s=2`, which gives a half-sum or
half-difference by inversion. The constant transform coordinate is zero,
so the columns involved are nonconstant. Conversely each listed vector
has norm one. The separate order-twelve checker supplies its classification.

## Every capacity-two assignment at sufficiently large prime order

Take integer `b>=5` copies of each of the `M-1` nonconstant columns.
At each row flip one entry, using distinct physical columns, and assign
at most two rows to any original column label. There are `r=b(M-1)`
physical columns. For an integer zero-sum `c`, let

```text
B(c)=sum_j |c^T s_j|/2.
```

First suppose `c` is ternary. Any selected flip changes one absolute
coefficient by at most one, so `B>=bML/2-h`. At `4<=h<=2sqrt(M)`, (1)
therefore gives for `b=5`

```text
B-r/2 >= f_M(h)/(h+2),
f_M(h)=(5/2)(h-2)M-6h^2+(11/2)h+5.                     (3)
```

For `M>=36`, this is positive throughout the interval. The quadratic is
concave in `h`, and its endpoint values are

```text
f_M(4)=5M-69>0,
f_M(2sqrt(M))=5x^3-29x^2+11x+5>0,   x=sqrt(M)>=6.
```

The last polynomial is positive at six and has positive derivative for
`x>=6`. Also `2sqrt(M)<q` in this range, so (1) applies. Increasing `b`
increases the lower bound because `ML>=M>M-1`.

For `h>=2sqrt(M)`, the existing high-support proof applies to every
capacity-two assignment. Explicitly, if `L=1+delta` and `k` selected
flips lower coefficient magnitude, then

```text
B=bML/2+h-2k,   k<=2s,
s<=M/h+M*delta/2+M*delta/h,
B-r/2 >= b/2+h-4M/h+M*delta*(b/2-2-4/h)>0.
```

Here `h>=12`, making every nonconstant term nonnegative. At `h=2`,
`L=1` and `B-r/2>=b/2-2>0`.

Finally, for a zero-sum integer character of amplitude `A>=2`, inversion
and the flip bound give

```text
B >= (bM/2) A - sum |c_i| >= (b/2-1)MA > b(M-1)/2.
```

Consequently every nonzero integer character satisfies `B>r/2`.

There is an unflipped copy of every label. If `lambda` is rational,
`sum lambda=0`, and all physical coefficients `lambda^T s_j` are integers,
comparison of a flipped column with an unflipped copy forces every
`2lambda_i` to be integral. Conversely any zero-sum integer `c=2lambda`
gives integral coefficients because sign columns equal one modulo two.
Thus the result includes the full saturated rational-row character lattice.

Choose distinct split primes with weights `tau<=log p_j<tau+1/r`,
`tau>1`, as in the [nearby-prime construction](strict_obtuse_prime_box_countermodels.md).
The same dyadic pigeonhole argument applies to every fixed `b` and
unbounded prime Paley order: partition `[X,2X]` into `r` intervals
with `X` a fixed positive multiple of `M^4`. The fixed-modulus prime
number theorem gives more than `r(r-1)` split primes there for large
`M`, so one interval contains `r` of them. Its logarithmic width is
less than `1/r`, and `tau=4log M+O_b(1)`.
Writing `V=sum_j |c^T s_j|log p_j/2` and `W=sum_j log p_j`, nonnegative
weights give, uniformly even for unbounded character amplitude,

```text
V-W/2 >= tau*(B-r/2)-1/2 >0.
```

The last strict inequality uses integer `B` and `B-r/2>=1/2`.
These sign profiles have literal Gaussian-integer realizations. Choose
one Gaussian prime `pi_j` above each `p_j` and set

```text
z_i=product_j pi_j^((1+s_ij)/2) conjugate(pi_j)^((1-s_ij)/2).
```

All rows have norm `product_j p_j`. Every physical column has both
signs, so the common Gaussian gcd is a unit. The unflipped Hadamard
columns distinguish the rows, so the points are distinct. For fixed
`b`, the nearby-prime construction gives `W=O(M log M)`. Their arguments
are uncontrolled: neither endpoint localization nor an endpoint
counterexample follows from this construction. The
[canonical skew-Hadamard argument](paley_canonical_all_character_half_height.md)
can instead exploit its specific flip assignment; (1)--(3) are independent
of that assignment and require a prime field.

For the canonical assignment, the later
[compatible-phase construction](paley_compatible_phase_coordinate_gap.md)
also enforces both coordinate lower bounds for every character using
complex factors and endpoint-scale phases. It does not enforce Gaussian
integrality of those factors or source rows.

The [exact checker](check_paley_primefield_polynomial_character_gap.py)
enumerates small prime orders, checks the polynomial identity and support
bound, checks the norm-one classification there, and tests the margin and
physical coefficient inequalities on deterministic larger fixtures.
