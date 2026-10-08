# Comparable pair squareclasses pass every six-point Witt obstruction

For arbitrarily large `H`, there are fifteen distinct primes `B_ij` in
`[H,2H]`, all congruent to one modulo eight, such that

```text
< epsilon_i product_(j!=i) B_ij : i=1,...,6 >
```

is a hyperbolic quadratic form over `Q`, for any assignment of three positive
and three negative signs `epsilon_i`.

Thus comparable pair-norm squareclasses are compatible with all rational
Witt and Hilbert-symbol conditions imposed by the central six-point moment
identities. This is not a full fair integer moment model. It does not provide
singleton factors of the required height, the relation `sum u_i=0`, a
controlled isotropic basis, circle nodes, or a least-radius estimate.

## 1. The necessary quadratic-form condition

Let six distinct rational unit-circle nodes `q_i=x_i+i y_i` have a nonzero
ordinary central vector `u` satisfying

```text
sum u_i = sum u_i q_i = sum u_i q_i^2 = 0.
```

For the rational bilinear form `beta(v,w)=sum u_i v_i w_i`, the three vectors
`1,x,y` are pairwise orthogonal and individually isotropic. Indeed,
the linear moments supply the products involving `1`, and

```text
sum u_i x_i^2 + sum u_i y_i^2 = sum u_i = 0,
sum u_i x_i^2 - sum u_i y_i^2 = Re(sum u_i q_i^2) = 0,
2 sum u_i x_i y_i = Im(sum u_i q_i^2) = 0.
```

They are independent: otherwise all six nodes lie on one affine line, whereas
a line meets the unit circle in at most two points. Hence `diag(u)` contains
a three-dimensional totally isotropic space and is hyperbolic. Here all
`u_i` are nonzero, as follows from the five-row Laurent Vandermonde kernel.

On the ideal two-level six-label cut profile, discard rational square
factors and write

```text
u_i = epsilon_i A_i^2 b_i,
b_i = product_(j!=i) B_ij.                           (1)
```

Singleton cuts occur squared, pair cuts to the first power, and balanced
triple cuts are absent. Rational isometry removes the `A_i^2`, so the Witt
condition sees only the signs and the fifteen pair squareclasses.

## 2. Exact local edge condition

Assume first that the `B_ij` are distinct odd primes, and put `a_i=epsilon_i b_i`.
The determinant is

```text
product_i a_i = - (product_(i<j) B_ij)^2.             (2)
```

It agrees modulo squares with three hyperbolic planes `H^3=<1,-1>^3`.
For an odd prime `p=B_ij`, write `a_i=p c_i`, `a_j=p c_j`. All other
coefficients are units. In the Hasse invariant

```text
h_p(a) = product_(r<s) (a_r,a_s)_p,
```

unit-unit pairs contribute one. For each `k` outside `{i,j}`, the two mixed
pairs contribute `(a_k/p)^2=1`. Only the edge pair remains:

```text
h_p(a) = (p c_i,p c_j)_p
       = Legendre(-c_i c_j,p)
       = Legendre(-epsilon_i epsilon_j
                   product_(ell!=i,j) B_iell B_jell, p).   (3)
```

We use the convention `(a,b)_p` for the quadratic Hilbert symbol and
`Legendre(a,p)` for the ordinary residue symbol. Formula (3) follows from

```text
(p^alpha u,p^beta v)_p
 = Legendre(-1,p)^(alpha beta)
   Legendre(u,p)^beta Legendre(v,p)^alpha.
```

Since `h_p(H^3)=1` at every odd prime, local hyperbolicity is exactly (3)
being `+1`. In particular, for `p=1 mod4`, the condition is

```text
Legendre(product_(ell!=i,j) B_iell B_jell, B_ij)=+1.  (4)
```

There are eight factors in this numerator. This also gives the condition
for pairwise coprime odd composite labels at a prime dividing one label to
odd order; even valuations impose no condition, and the unit part of that
label occurs squared and cancels.

For completeness, the local classification used here says that dimension,
determinant modulo squares, and Hasse invariant determine a quadratic form
over `Q_p`; see [Milne, Class Field Theory, Theorem 6.10 in the section on quadratic forms](https://www.jmilne.org/math/CourseNotes/CFT.pdf).
In this particular case one can also check sufficiency directly: condition
(3) makes the two-dimensional `p`-part hyperbolic, while the remaining
four-dimensional unit form has square determinant and is hyperbolic by
finite-field classification and Hensel lifting.

## 3. Comparable labels satisfying all local conditions

Take the primes `p=1 mod8` in `[H,2H]`. Color the unordered pair `{p,q}`
by `Legendre(p,q)`. Quadratic reciprocity makes this an undirected two-color
graph. Once the number of available primes is at least `R(15,15)`, Ramsey's
theorem supplies fifteen primes whose mutual symbols all equal one common
sign `sigma`.

Only the elementary finite Ramsey bound is needed:
`R(15,15)<=binom(28,14)`. The usual recurrence
`R(a,b)<=R(a-1,b)+R(a,b-1)` follows by splitting the neighbors of one vertex
according to color; its induction yields this binomial bound.

There are enough primes in the fixed progression for every sufficiently
large `H`, since their count in `[H,2H]` is asymptotic to `H/(4 log H)`.
The needed fixed-modulus prime number theorem is stated in
[Sutherland's MIT 18.785 notes, Lecture 18](https://ocw.mit.edu/courses/18-785-number-theory-i-fall-2021/mit18_785f21_full_lec.pdf).

Assign the fifteen selected primes bijectively to the edges of `K6`.
Every condition (4) is then `sigma^8=1`. At all other odd primes the
coefficients are units, so the Hasse invariant is one. At two, each `b_i`
is one modulo eight and therefore a square in `Q_2`; the form is visibly
isometric to `<epsilon_1,...,epsilon_6>=H^3`. At infinity its signature
is `(3,3)`. Thus it is locally hyperbolic everywhere and hence hyperbolic
over `Q` by Hasse--Minkowski. For a primary expository statement of the
local-global theorem, see [Conrad, The Local--Global Principle, Theorem 3.3](https://kconrad.math.uconn.edu/blurbs/gradnumthy/localglobal.pdf).
Applying that theorem to an isotropic vector, splitting a hyperbolic plane,
and repeating proves the required three-plane conclusion.

The logarithms of all fifteen pair labels are `log H+O(1)`. Consequently
the Witt obstruction alone cannot force a gap among their logarithmic sizes.
The construction makes no assumption about whether the monochromatic sign
is positive or negative; its even exponent is essential.

## 4. What the construction deliberately leaves unresolved

Hyperbolicity supplies rational isotropic vectors and rational maximal
isotropic spaces. It does not supply the particular vector `1` in an
isotropic space for the coefficients in (1). To obtain the zeroth central
moment from these squareclasses, one still needs

```text
sum_i epsilon_i b_i A_i^2 = 0.                       (5)
```

Although the split form has nonzero rational solutions, clearing their
denominators can introduce factors and heights incompatible with the desired
singleton blocks. Prescribed comparable singleton factors, pairwise conductor
coprimality, all remaining moments, and the angular/radius scale are additional
arithmetic constraints. In particular, assigning extra formal labels for the
six singleton and ten balanced cuts does not construct actual circle points.

## 5. Exact finite witness

The following fifteen primes lie in `[10^6,2*10^6]`, are all one modulo eight,
and have all 105 mutual Legendre symbols equal to `+1`:

```text
1037329, 1584721, 1379513, 1501889, 1869649,
1278713, 1700513, 1209353, 1118009, 1377881,
1242641, 1379449, 1636121, 1207417, 1747169.
```

Run `python3 docs/check_six_point_witt_pair_squareclasses.py`. It checks their
primality, comparability, all mutual residue symbols, determinant, and all
local invariants against `H^3`. It also checks the sign-sensitive general
formula (3) on 600 edge instances without restricting primes modulo eight.
The infinite comparable construction follows from Ramsey and the prime
number theorem, not from the finite witness.
