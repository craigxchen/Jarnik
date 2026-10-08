# Exact torus relations among balanced half-characters

The 35 primitive Gaussian half-characters of eight common-unit rows obey
many multiplicative relations. Their rational content can be computed
exactly, even when several nested threshold layers use the same split
prime. After this reduction, a zero row-character relation has primitive
factor one. A nonzero combination is precisely an ordinary integer row
character. Thus the multiplicative relations alone do not supply an
additional axis-gap height inequality. In particular, all resulting
scalar inequalities pass the positive full 127-cut formal profile.

## Exact primitive-product identity

Remove the complete common Gaussian factor. For each varying split prime
`p`, choose `pi_p` above `p`, and let `a_i(p)` be its allocation exponent
at row `i`. Let `Lambda` be any choice of 35 balanced sign vectors modulo
global negation. For `lambda in Lambda`, write

```text
c_lambda(p) = sum_i lambda_i a_i(p),
A_lambda = product_p pi_p^max(c_lambda(p),0)
                     conjugate(pi_p)^max(-c_lambda(p),0).
```

For an arbitrary integer vector `n=(n_lambda)`, put

```text
L = sum_lambda n_lambda lambda,
q_p = sum_i L_i a_i(p) = sum_lambda n_lambda c_lambda(p),
G_n = product_lambda A_lambda^max(n_lambda,0)
                     conjugate(A_lambda)^max(-n_lambda,0),
A_L = product_p pi_p^max(q_p,0)
                conjugate(pi_p)^max(-q_p,0).
```

Then the **literal Gaussian-integer identity** is

```text
G_n = C_n A_L,
C_n = product_p p^k_p,
k_p = (sum_lambda |n_lambda| |c_lambda(p)| - |q_p|)/2
    = min(v_pi_p(G_n),v_conjugate(pi_p)(G_n)) >= 0.       (1)
```

The exponents `k_p` are integers: the numerator is even because
`|x|=x (mod 2)`. No assumption of disjoint blocks is used; each `a_i(p)`
may be a sum of arbitrarily many nested layers. In particular

```text
product_lambda (A_lambda/conjugate(A_lambda))^n_lambda
  = A_L/conjugate(A_L).                                  (2)
```

If `L=0`, then `A_L=1` and `G_n=C_n` is a positive rational integer.
For two combinations `n,m` with the same `L`, (1) gives the exact
content-preserving exchange `C_m G_n=C_n G_m`.

For example, take a two-element set `U` and four other distinct labels
`a,b,c,d`. Let `lambda_X` have negative set `X`. The balanced four-subset
exchange

```text
lambda_(U union {a,b}) + lambda_(U union {c,d})
 = lambda_(U union {a,c}) + lambda_(U union {b,d})        (3)
```

therefore gives

```text
C_right A_(Uab) A_(Ucd) = C_left A_(Uac) A_(Ubd),        (4)
```

with `C_left` and `C_right` computed by (1). The content can differ.
At one split prime, use `U={0,1}`, `(a,b,c,d)=(2,3,4,5)` and allocation
`(a_0,...,a_7)=(0,0,5,4,1,0,0,0)`. Its threshold layers have minority
sets `{2,3,4}`, `{2,3}` repeated three times, and `{2}`. The four
signed exponents in (3) are `(-8,+8,-2,+2)`. Consequently

```text
A_(Uab) A_(Ucd) = p^8,
A_(Uac) A_(Ubd) = p^2,
C_left/C_right = p^6.                                     (5)
```

This example shows why an exchange cannot be treated as literal equality
of the unnormalized `A` products.

## The phase and height information left after reduction

With one common source unit and `sum_i L_i=0`, the row character is
exactly

```text
product_i z_i^L_i = A_L/conjugate(A_L).                 (6)
```

The common factor and common unit cancel. If all row arguments lie in
an interval of width `Delta`, then

```text
dist(arg A_L, pi Z) <= (||L||_1/4) Delta.                (7)
```

For a conjugate-primitive nonunit `A_L`, its imaginary coordinate after
axis rotation is a nonzero integer. If `Delta<=C/sqrt(R)` and
`W=log R^2`, (7) yields the usual scalar half-character bound

```text
log Norm(A_L) >= W/2 - 2 log(C ||L||_1/4).               (8)
```

For `L=0`, the phase is identically zero and there is no nonunit to
which an axis gap can be applied. For `L!=0`, the precise height is

```text
log Norm(A_L) = sum_p |sum_i L_i a_i(p)| log p.          (9)
```

Thus (8) is exactly the scalar row-character inequality. Taking more
torus relations does not create a new primitive Gaussian integer: (1)
reduces every monomial to one `A_L` before the integer-coordinate step.
This is a limitation of *multiplicative relations plus the scalar axis
gap*, not of arguments using additive Gaussian identities or other
arithmetic structure.

## Full 127-cut formal test

Take one independent split-prime block for every nonconstant unoriented
cut of eight rows, with one representative from each complementary pair.
There are 127 cuts. At equal formal log weight `h`, `W=127h` and for every
nonzero integer `L` with `sum_i L_i=0`, (9) is

```text
log Norm(A_L) = h sum_cuts |sum_(i in cut) L_i| >= 64h.  (10)
```

To prove the bound, sum first over all 256 subsets `T`. Complementary
subsets have equal absolute sums, and the empty/full pair contributes
zero, so the cut sum is half the subset sum. Choose any `i` with
`|L_i|>=1` and pair every subset excluding `i` with the subset obtained
by adding `i`. Each pair contributes at least `|L_i|` by the triangle
inequality. There are 128 pairs, giving (10). A row difference
`L=e_i-e_j` attains 64.

The scalar lower bound (8) asks for only `W/2=63.5h` at leading order.
This margin persists for literal distinct split primes: choose odd
powers so that every block weight is `h+O(1)` as `h` grows. Then both the
minimum and maximum block weights, divided by `h`, tend to one. In fact
`64 min(w)-127 max(w)/2` tends to infinity. Equation (10) therefore
implies `log Norm(A_L)>W/2` for all nonzero `L` simultaneously when
`h` is large. For every fixed `C`, it also passes the full finite bound
(8) for all such `L`: their `l1` norms are at least two, so the largest
possible additive term on its right is `2 log(2/C)`, independent of `L`.
This is a formal height test; no endpoint arc realization
is asserted. It shows that even the entire family of scalar inequalities
obtained from torus monomials does not exclude this positive profile.

The common quadratic twist and squareclass rank of this same profile
are treated in [the odd-layer obstruction](odd_layer_character_quadratic_twist_obstruction.md).
The separate even-layer theorem uses actual Gaussian quarter-root
integers, which are unavailable here.

[The exact checker](check_balanced_half_character_torus_relations.py)
checks 120 literal nested-allocation Gaussian content identities, the
displayed exchange, and all 78,124 nonzero row vectors whose first seven
entries lie in `[-2,2]` and whose last entry makes the sum zero. The
subset-pairing argument above proves the unrestricted scalar bound.
