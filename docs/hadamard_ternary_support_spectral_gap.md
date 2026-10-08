# A high-support gap for ternary Hadamard characters

For a repeated Paley one-flip profile with an arbitrary label-capacity-two
assignment, the earlier fourth-moment argument left ternary characters
above a square-root support range unresolved.
An elementary Hadamard argument excludes the opposite, high-support end
for **all** zero-sum ternary characters of support at least `2sqrt(M)`
when each old label receives at most two assigned flips. The first step
is a transform inequality:
for any normalized Hadamard matrix, every ternary vector `c` of support
`h` satisfies

```text
sum_a |(H^T c)_a| >= 2h,                               (1)
```

unless `c` is exactly a positive or negative full column of `H`. The
stronger high-support conclusion also covers those full columns. No
Paley character-sum estimate is used. This does not prove a growth bound
for actual endpoint configurations; it narrows the possible short
characters that might supply one.

## Exact transform inequality

Write `T=H^T c`, `P=max_a |T_a|`, and `h=|supp(c)|`. Hadamard
orthogonality gives `sum_a T_a^2=Mh`. If `P<=M/2`, then

```text
sum_a |T_a| >= (sum_a T_a^2)/P = Mh/P >= 2h.
```

If `M/2<P<M`, orient a maximizing column `a` so that `T_a=P`.
On the support of `c`, there are `d=(h-P)/2` places where `c` and
`H_a` disagree; off the support there are `k=M-h` places. For every
other column `b`, orthogonality of `H_a` and `H_b` gives

```text
|T_b| <= k+2d = M-P.                                   (2)
```

Parseval and (2) yield

```text
sum_(b!=a)|T_b|
 >= [Mh-P^2]/(M-P)
 = h+P(h-P)/(M-P).
```

Consequently

```text
sum_a|T_a|-2h
 >= (h-P)(2P-M)/(M-P) >=0.                       (3)
```

If `P=M`, equality in the maximum forces `h=M` and `c=+/-H_a`;
these are the sole exceptions. Thus (1) holds for every other ternary
vector, with or without zero row sum. In normalized notation
`L(c)=sum|T_a|/M`, it says `L(c)>=2h/M`.

## Repeated-column one-flip consequence

Now take `b>=5` physical copies of the `M-1` nonconstant Hadamard
columns, with one distinct physical one-entry flip assigned per row.
Assume at most two rows are assigned to each old column label. For a
nonzero zero-sum ternary vector `c`, let
`B(c)=sum_j |c^T s_j|/2`, summed over the `r=b(M-1)` physical columns.
Put `L=sum_a|T_a|/M=1+delta`, and let `s` be the number of labels with
`T_a!=0`. Hadamard inversion ensures `L>=1`. Let `k` count supported
rows whose assigned flip label has nonzero transform coefficient and
whose sign agrees with that coefficient. A supported flip lowers the
coefficient absolute value by one exactly in this case; every other
supported flip raises it by one. Hence the **exact** count is

```text
B(c)=(bM/2)L+h-2k,        k<=2s.                  (4)
```

Every nonzero `|T_a|` is an even integer between `2` and `h`.
The elementary inequality `(t-2)(h-t)>=0` gives
`t^2<=(h+2)t-2h` on this set. Sum over the `s` active labels and use
Parseval `sum T_a^2=Mh` to obtain

```text
s <= M/h + M delta/2 + M delta/h.                   (5)
```

Substituting (5) and `k<=2s` into (4) gives

```text
B(c)-r/2
 >= b/2+h-4M/h + M delta (b/2-2-4/h).            (6)
```

For `b>=5` and `h>=8`, the coefficient of `delta` is nonnegative.
If also `h>=2sqrt(M)`, then `h-4M/h>=0`, so (6) is at least `b/2>0`.
This includes a full-column character: it has `s=1`, `delta=0`,
and (4)--(6) still apply. The all-ones column cannot be zero-sum.
Thus every zero-sum ternary character of support `h>=2sqrt(M)` is
above the primitive half-height threshold for actual distinct nearby
split primes, provided `M>=12`. The perturbation of prime log weights
by less than `1/r` changes its half-height excess by more than `-1/2`,
while `B-r/2>=b/2` at this support, so the strict phase-height gap is
uniform at large prime scale.

Combining only (6) with the previously proved
[Paley square-root support gap](paley_sqrt_support_fourth_moment_gap.md),
any possible cheap saturated character in the five-copy Paley family
would have to be ternary with support roughly
`(1/3)sqrt(M)<h<2sqrt(M)`. The subsequent
[prime-field polynomial argument](paley_primefield_polynomial_character_gap.md)
closes this remaining interval for prime Paley orders `M>=36`, proving
that every character exceeds half-height for every capacity-two
assignment. This does not place the actual Paley points on a short arc.
Simultaneous source-phase constraints remain outside these individual
height tests.

For the specific diagonal Paley assignment, the
[skew-Hadamard theorem](paley_canonical_all_character_half_height.md)
goes further: it proves the exact minimum `B=bM/2-2` for **every**
saturated integer character at every Paley order. The support theorem
here remains useful for arbitrary Hadamard matrices and assignments
with label capacity two.

The [finite checker](check_hadamard_ternary_support_spectral_gap.py)
tests (1), (4), (5), and (6) on every zero-sum ternary vector at Walsh
order eight and Paley order twelve, with the high-support conclusion
checked for the physical one-flip assignment.
