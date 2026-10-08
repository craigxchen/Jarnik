# Multiquadratic near-squares: exact trace lattices and the maximal-order gap

This note examines a field whose degree grows with the number of
primitive norms. It proves an exact trace-lattice identity, explains
why the obvious conjugate products do not amplify the small real
residues, and records the substantial maximal-order issue separately.
The latter has a concrete local obstruction: the near-square defects
are units at shared conductor primes, so dividing their products by
the common factors used to normalize radical products is generally
not integral. These statements do not rule out every field argument
and do not prove the uniform endpoint bound.

## 1. Near-square defects and all their conjugates

Let

```text
q_i=x_i^2+t_i^2,       x_i>0,       t_i in Z minus {0},
E=Q(sqrt(q_1),...,sqrt(q_m)),       n=[E:Q].
```

Changing the sign of a Gaussian numerator if necessary makes `x_i`
positive without changing its phase ratio or norm. In the central
system, `log x_i=W/4+o(W)` and `log max(1,|t_i|)=o(W)`.
Assume each norm squareclass is nontrivial. Define the algebraic
integer

```text
delta_i=sqrt(q_i)-x_i
       =t_i^2/(sqrt(q_i)+x_i).
```

It has a very small positive value in the distinguished real
embedding. Nevertheless its quadratic conjugate is
`-sqrt(q_i)-x_i`, and exactly

```text
Norm_(E/Q)(delta_i)=(-t_i^2)^(n/2).              (1)
```

For any integer exponents `a_i`, multiplicativity gives the exact
norm of `product_i delta_i^a_i`; no additional smallness appears
when more rows are multiplied. In particular, when every `|t_i|=1`,
all the `delta_i` are algebraic units, however small their chosen
real embeddings are.

This compensation can be quantified without estimating a product
of `n` individual conjugates crudely. Write

```text
lambda_i=log((sqrt(q_i)+x_i)/|t_i|),
sigma(sqrt(q_i))=epsilon_i(sigma) sqrt(q_i).
```

Then

```text
log|sigma(delta_i)|=log|t_i|-epsilon_i(sigma)lambda_i.  (2)
```

The `epsilon_i` are the corresponding Galois characters. If the
norm classes are nontrivial and pairwise distinct, character
orthogonality gives

```text
(1/n) sum_sigma epsilon_i(sigma)=0,
(1/n) sum_sigma epsilon_i(sigma)epsilon_j(sigma)=0  (i!=j).
```

Consequently the mean of the logarithmic conjugates of
`product delta_i^a_i` is `sum a_i log|t_i|`, while their normalized
variance is exactly `sum a_i^2 lambda_i^2`. The increasing field
degree does not alter these normalized formulas. The small chosen
embedding is balanced by the other embeddings, rather than supplying
many simultaneous small embeddings for free.

## 2. An exact trace-lattice theorem

Let `F` be a downward-closed family of subsets of `[m]`. Assume
the squareclasses of `product_(i in A)q_i`, for `A in F`, are
pairwise distinct. Define

```text
r_A=product_(i in A) sqrt(q_i),
b_A=product_(i in A) delta_i.
```

Then the two integral lattices inside `E` are **identical**:

```text
sum_(A in F) Z b_A = sum_(A in F) Z r_A.         (3)
```

Indeed, expansion of each product `b_A=product(sqrt(q_i)-x_i)`
uses only `r_B` with `B subset A`, has leading coefficient one,
and has integer coefficients. Downward closure makes this an
integer triangular matrix with determinant one. The inverse
expansion `sqrt(q_i)=delta_i+x_i` is also integral.

The radical basis is orthogonal for the positive definite trace
inner product. For `A!=B`, the squareclass of the symmetric
difference is nontrivial, so an automorphism negates that radical
and its trace is zero. For `A=B`, its trace is
`n product_(i in A)q_i`. If `r=|F|` and
`c_i=#{A in F:i in A}`, the exact Gram determinant is therefore

```text
det(Tr_(E/Q)(b_A b_B))_(A,B in F)
   = n^r product_i q_i^c_i.                    (4)
```

This retains the actual norm-factor divisibility: an odd prime's
exponent on the right is exactly `sum_i c_i v_p(q_i)`. A shared
cut prime contributes with every one of its row incidences, with
no omitted square factors.

Six-wise independence permits

```text
F={A:|A|<=3},
r=sum_(j=0..min(3,m)) binom(m,j),
c_i=sum_(j=0..min(2,m-1)) binom(m-1,j).          (5)
```

If all `m` norm classes are independent, one may take every subset,
so `r=n=2^m` and

```text
Disc(Z[sqrt(q_1),...,sqrt(q_m)])
   =n^n product_i q_i^(n/2).                   (6)
```

Thus the determinant does not become smaller when its basis is
written as products of tiny near-square defects. This is an exact
statement even when `n` grows exponentially in `m`; it does not
rely on a fixed-degree auxiliary-polynomial estimate.

The same character calculation gives, for any set `A` of size at
most six under six-wise independence,

```text
Tr_(E/Q)(product_(i in A)delta_i)
 =n (-1)^|A| product_(i in A)x_i.               (7)
```

These elementary traces are large, and the norms are already fixed
by (1). Neither operation makes the distinguished small values
into an unaccounted small rational integer.

## 3. The maximal order is substantially larger

Equation (6) is the discriminant of the displayed order, **not**
the discriminant of the full ring of integers. This distinction is
large in precisely the shared-cut setting and cannot be dismissed.

For a clean exact calculation, suppose each `q_i` is squarefree,
all its prime factors are `1 mod 4`, and the norm classes are
independent. Let `P_ram` be the union of the prime supports of the
`q_i`; primes in an empty core block are absent from this union.
Every quadratic subfield has squarefree defining integer congruent
to one modulo four, so the compositum is unramified at two.
At each odd ramified prime its inertia has order two and is tame.
The tame different exponent is therefore `e-1=1`, which gives
the exact field discriminant

```text
|Disc(E)|=product_(p in P_ram) p^(n/2).          (8)
```

At a prime whose incident set of rows is `S_p`, the logarithm of
the maximal-order index obtained from (6)--(8) is

```text
log[O_E:Z[sqrt(q_1),...,sqrt(q_m)]]
 = (n/2)log n
   +(n/4)sum_(p in P_ram)(|S_p|-1)log p.        (9)
```

The shared-prime terms in (9) are real and large. For example, the
radical product over a subset `A` can be divided integrally by

```text
s_A=product_p p^floor(|A intersect S_p|/2),
r_A/s_A=sqrt(squarefree_part(product_(i in A)q_i)).
```

Those divisions remove the odd-prime part of the index. Remaining
powers of two account for the usual integral half-radical bases.
It would be incorrect to claim that (3) alone controls what happens
after these maximal-order normalizations.

## 4. Why the obvious normalized near-square products are not integral

There is an exact local obstruction to applying those same divisions
to the tiny defects. Suppose an odd prime `p` divides `q_i` and
does not divide `t_i`. In the arithmetic setting
`h_i=x_i+i t_i` is conjugate-primitive, so the latter condition is
automatic: `p|q_i,t_i` would force `p|x_i`, contradicting
primitivity. Thus `x_i` is also a `p`-unit.

At **every** prime ideal of `E` above `p`, `sqrt(q_i)` has positive
valuation, whereas `x_i` is a unit. Therefore

```text
delta_i=sqrt(q_i)-x_i is a unit at every prime above p. (10)
```

If `p` divides two incident norms `q_i,q_j`, the radical product
`sqrt(q_i q_j)/p` is integral, but

```text
delta_i delta_j/p is not integral.              (11)
```

Indeed the numerator in (11) is a unit at every prime above `p`.
This holds with the full relevant prime valuations, not merely
with their supports. More generally, products of defects from
rows all incident to a given cut prime are units at that prime.

The correct integral expression reintroduces the large coordinates:

```text
sqrt(q_i q_j)/p
 =(delta_i delta_j+x_i delta_j+x_j delta_i+x_i x_j)/p. (12)
```

The rational term `x_i x_j/p` has a genuine denominator `p`, since
both coordinates are units there. If one divides by the entire
common norm factor, none of its conductor primes cancels against
`x_i x_j`. In a uniform profile that denominator has logarithmic
size of order `W/4`, rather than negligible height.

Consequently the maximal-order index is a possible source of new
arithmetic information, but simply dividing products of small
defects by its shared-prime factors is invalid. A successful
maximal-order argument would have to construct different integral
combinations and control their remaining large-coordinate terms.
No such construction is supplied by (3) or (9).

For a finite example take `q_1=65=8^2+1` and
`q_2=85=9^2+2^2`. The field has degree four, its discriminant is
`(5*13*17)^2`, and the raw radical order has index `80` in its
maximal order. The radical `sqrt(65*85)/5=sqrt(221)` is integral.
But for `delta_1=sqrt(65)-8` and `delta_2=sqrt(85)-9`,
`Norm(delta_1 delta_2)=16`, so
`Norm(delta_1 delta_2/5)=16/625` is not an integer. This explicitly
exhibits the failure of the analogous small-defect normalization.

## 5. A concrete obstruction to field degree alone

For every fixed `m`, there are infinitely many tuples of genuine
conjugate-primitive Gaussian integers

```text
h_i=X+c_i+i,       t_i=1,
q_i=(X+c_i)^2+1,
[Q(sqrt(q_1),...,sqrt(q_m)):Q]=2^m,              (13)
```

with distinct fixed even integers `c_i` and arbitrarily large even
`X`. Thus exact unit-sized imaginary parts coexist with the largest
possible multiquadratic degree and balanced row heights.

To ensure the last assertion, choose distinct split primes `p_i`
larger than every `(c_i-c_j)^2+4`. Prescribe
`X+c_i` to be a root of `-1 mod p_i`, then lift the residue modulo
`p_i^2` so that `v_(p_i)(q_i)=1`. For `j!=i`, the inequality on
`p_i` ensures `p_i` does not divide `q_j`: a common root would
force `p_i|(c_i-c_j)^2+4` or `p_i|c_i-c_j`. Chinese remaindering
these prescriptions with `X=0 mod 2` gives an infinite arithmetic
progression. Each norm has a private odd-valuation prime, proving
independence of all the squareclasses.

This is not a full uniform-cut example. In fact, for `d=c_i-c_j`,

```text
gcd(q_i,q_j) divides d^2(d^2+4),
```

by eliminating `X` from the two quadratic norm polynomials.
Consequently pairwise common norm masses are bounded for fixed
`m`, and the leading conductor pattern is concentrated on singleton
rows. Its least-radius reconstruction including the anchor has
radius comparable to `X^m`; its angular distance from the anchor
is of order `1/X`, so it misses the endpoint when `m>=3`.

The conclusion is correspondingly precise. Large field degree,
six-wise independence, the natural field norms and traces, and the
unmodified near-square trace lattice do not by themselves produce
a new exponent. Shared-cut data survive in the large maximal-order
index, but exploiting that index requires integral combinations
beyond the invalid normalization (11). This is the remaining field
route, rather than a proof that all Galois-symmetry arguments fail.

There is now a valid mixed normalization in
[maximal_order_pair_defects.md](maximal_order_pair_defects.md):
adding the mixed terms and the small rational correction
`-t_i t_j` before dividing by the full shared norm factor produces
the integral primitive pair defect `sqrt(q_ij)-x_ij`. The new note
computes its lattice index and its matching size cost. Further
graph traces, three-radical spin defects, and fixed-character sign
differences are treated in
[graph_trace_and_spin_defects.md](graph_trace_and_spin_defects.md).
Thus the obstruction in (11) applies to the unmodified product,
not to all combinations with the same nominal denominator.

As finite verification, 15 exact integer trace-Gram determinants
were evaluated for downward-closed degree truncations through five
independent radicals, including full field degree 32. They agreed
with (4). The index-80 example above also satisfies the exact
discriminant identity. These checks supplement the proofs and are
not Lean formalizations.
