# Odd balanced characters: the shared quadratic field is at the wrong level

The normalized **half-characters** of eight common-unit Gaussian points
do share one rational quadratic twist. Their **quarter-root elements**,
which would be needed to repeat the even-layer height proof, generally do
not share a quadratic extension of `Q(i)`. An exact positive cut-profile
fixture has quarter-root compositum degree `2^21` over `Q(i)`.

Moreover, even a fixed common quadratic extension does not supply the
integer perpendicular-coordinate bound used in
[the even-layer theorem](even_column_balanced_character_height.md): a
Pell family makes that coordinate arbitrarily small, with the weaker
quadratic approximation exponent. These are precise obstructions to
this particular extension of the proof. They neither disprove a
simultaneous character inequality using additional relations nor give
an actual eight-point endpoint counterexample. The uniform arc problem
remains open.

## 1. Exact common twist, allowing arbitrary nested layers

Remove the complete common factor, choose a Gaussian prime `pi_p`
above each varying split prime, and write the allocation in row `i` as
`a_i(p)`. For a balanced `lambda in {+1,-1}^8`, put

```text
c_lambda(p) = sum_i lambda_i a_i(p),
A_lambda = product_p pi_p^max(c_lambda(p),0)
                       conjugate(pi_p)^max(-c_lambda(p),0).
```

This is conjugate-primitive. With one common row unit, exactly

```text
product_i z_i^lambda_i = A_lambda / conjugate(A_lambda).       (1)
```

Every `lambda_i` is odd, so the parity

```text
c_lambda(p) = sum_i a_i(p) (mod 2)
```

is independent of the character. Consequently, for the squarefree integer

```text
t = product_(p : sum_i a_i(p) odd) p,
```

there are positive integers `h_lambda` with

```text
Norm(A_lambda) = t h_lambda^2,
r_lambda = A_lambda / (h_lambda sqrt(t)),
r_lambda^2 = product_i z_i^lambda_i.                          (2)
```

Thus all `r_lambda` belong to `Q(i,sqrt(t))`. This conclusion is literal
and survives cancellation among overlapping prime layers. The field is
`Q(i)` when `t=1`; otherwise it is quadratic over `Q(i)`.

If the source arguments occupy an interval of width `Delta`, (1) gives

```text
dist(arg A_lambda, pi Z) <= 2 Delta.
```

The desired fourth-root analogue is `beta_lambda=sqrt(A_lambda)`, whose
argument is within `Delta` of `(pi/2)Z`. If `t=1`, every `c_lambda(p)` is
even and these square roots can be taken in `Z[i]`. If `t>1`, the
Gaussian squareclass of `A_lambda` matters, not just its norm squareclass.

For every prime dividing `t`, its contribution to `[A_lambda]` in
`Q(i)^*/Q(i)^{*2}` is `[pi_p]` when `c_lambda(p)>0` and
`[conjugate(pi_p)]` when `c_lambda(p)<0`. These differ by the nontrivial
rational squareclass `[p]`. All other prime contributions are squares.
Hence sharing the norm twist `t` does not make the Gaussian squareclasses
equal. A single quadratic extension containing all `sqrt(A_lambda)`
would require their nonzero squareclasses to span dimension at most one.

## 2. A positive formal profile with exact field degree `2^21`

Take every unoriented nonconstant cut on eight labelled rows once:

```text
minority size             1    2    3    4
number of cuts           8   28   56   35
```

For size four use the representative containing row one. Assign each
cut an independent Gaussian prime block. This is a literal Gaussian
factorization, but **no endpoint arc realization is asserted**.

At unit formal log weights, its total weight is `127` and every pair
of rows is separated in `64` columns. Thus it passes strict pair
separation `64>127/2`, and its sign matrix satisfies

```text
SS^t = 128 I - 1 1^t.
```

In particular the allocation rows have full affine rank. It also has
positive value for the previously targeted intrinsic expression:

```text
7 S_1 + 2 S_2 - S_3 - S_4
  = 7*8 + 2*28 - 56 - 35 = 21 > 0.                         (3)
```

Equal log weights here describe a formal profile, not independent primes
of equal norm. For a literal unbounded family approaching these weights,
fix 127 distinct split primes and replace each block by an odd power:
choose odd `e_p` with `e_p log p=L+O(log p)` as `L` tends to infinity.
Both strict inequalities persist and the squareclasses below are unchanged.
This also realizes canonical nested layers, by repeating each cut `e_p`
times.

Use the 35 characters with `lambda_1=+1`. Set `x_1=0` and
`x_i=(1-lambda_i)/2` for `i=2,...,8`; the remaining seven bits have
weight four. In a singleton column, the negative-orientation indicator
is `x_i`. In a triple column `{i,j,k}`, it is majority of those bits:

```text
m_ijk = x_i x_j + x_i x_k + x_j x_k over F_2.              (4)
```

Even cuts contribute zero to the squareclass vector. For each odd cut,
its two Gaussian-prime coordinates are `1+m_T` and `m_T`. Therefore
the squareclass matrix has the same rank as the evaluation matrix of
the constant function and the functions `m_T`.

Triples containing row one supply every quadratic `x_i x_j` on seven
bits. Conversely (4) expresses every triple through these quadratics,
and on the weight-four slice

```text
x_i = sum_(j != i) x_i x_j.
```

The 21 squarefree quadratic monomials have exactly one linear relation
on this slice, their total sum being `binom(4,2)=0` modulo two. To prove
there are no others, suppose `q(S)=sum_{i<j in S} b_ij` vanishes on all
four-subsets. Comparing `T union {u}` with `T union {v}` gives

```text
sum_(j in T) (b_uj+b_vj)=0
```

for every three-subset `T` of the other five indices. Comparing two such
triples makes all five summands equal; their three-term sum then forces
each to be zero. Thus `b_uj=b_vj` for distinct indices, and all edge
coefficients are equal. The same argument applied to a constant
evaluation shows that constant `1` is not in their span: equal edge
coefficients evaluate to zero.

The exact rank is consequently `20+1=21`. Unique Gaussian-prime
valuations make these 21 squareclasses independent in
`Q(i)^*/Q(i)^{*2}`. Repeated quadratic adjunction therefore gives

```text
[Q(i, sqrt(A_lambda) : lambda_1=+1) : Q(i)] = 2^21.       (5)
```

Already the 56 triple cuts alone give the same rank. The choice of one
representative per conjugate pair is specified in (5); no claim of
minimal degree under independently conjugating representatives is needed.

This degree is still bounded in terms of the fixed eight rows, independently
of the radius. Equation (5) is an obstruction to a **single quadratic
quarter-root field**, not to arbitrary bounded-degree-field methods. The
prime supports and hence field discriminants can vary with the profile;
no uniform lattice covolume follows from a bound on the degree alone.

## 3. Even one fixed twist loses the integer axis gap

Let `G=2+i`, of norm five. For any positive Pell solution

```text
u_n + y_n sqrt(5) = (9+4 sqrt(5))^n,
x_n = u_n - 2 y_n,
B_n = x_n+i y_n,
A_n = G B_n^2.
```

Direct expansion gives

```text
Im(A_n) = x_n^2+4 x_n y_n-y_n^2
        = u_n^2-5 y_n^2 = 1,
Norm(A_n) = 5 Norm(B_n)^2.                               (6)
```

In particular every `A_n` is conjugate-primitive, all have the same
Gaussian squareclass `[G]`, and all `beta_n=sqrt(G) B_n` lie in the
single quadratic extension `Q(i,sqrt(G))`. For `n>=1`, `Re(A_n)<0`, so
the closest target axis for `beta_n` is imaginary. Since
`|A_n|` tends to infinity,

```text
dist(arg beta_n,(pi/2)Z)
    = (1/2) arcsin(1/|A_n|)
    ~ 1/(2 |beta_n|^2).                                  (7)
```

Thus `|beta_n| sin(dist)` tends to zero, contradicting any positive
uniform perpendicular-coordinate gap in this fixed twist. The field
norm records the Pell cancellation and cannot supply the Gaussian-axis
bound `|beta| sin(dist)>=1`.

Taking 35 consecutive Pell indices gives 35 distinct simultaneous
examples in that same field, with bounded ratios of their heights and
the same failure. They are **not claimed to obey the multiplicative
relations of 35 balanced characters of eight source points**. Their
purpose is precisely to show that common-field membership, integrality,
primitivity, distinctness, and simultaneous small arguments alone do
not repair the proof. Additional source-character structure would have
to do actual work.

## Verification

[The exact checker](check_odd_layer_character_quadratic_twist_obstruction.py)
checks the complete cut census, Gram and separation identities, the
rank-21 squareclass certificate, all 35 literal common-twist identities
with independent Gaussian primes, and the first 35 Pell identities.
The rank computation supplements the combinatorial proof above.
