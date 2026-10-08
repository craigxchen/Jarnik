# Reciprocal-minor lcm as an exact lattice-content invariant

This note sharpens the divisibility `B_k | N(A)` in
`luna_fresh_reciprocal_minor_symmetry.md`.  The quotient is exactly a lattice
index.  The identity keeps every rational prime power, including 2,
and isolates what a reciprocal-minor height argument would still have to
control.

Let

```text
w_j=-d+i t_j,                  A=lcm_G(w_0,...,w_{k-1}),
q_j=A/w_j,                     M=N(A),
n_j=d^2+t_j^2,
delta_ij=det(q_i,q_j)=M d(t_j-t_i)/(n_i n_j).
```

The `t_j` are distinct, so the `q_j` span a rank-two rational lattice.  Put

```text
h = gcd_(i<j) |delta_ij|,
Lambda = sum_j Z q_j  subset Z[i]=Z^2.
```

Then `h=[Z^2:Lambda]`: this is the usual gcd-of-maximal-minors formula for
the index of a full-rank integer column lattice.

## 1. Exact formula for the denominator lcm

Reduce

```text
d(t_j-t_i)/(n_i n_j)=a_ij/b_ij,       gcd(a_ij,b_ij)=1,
```

with `b_ij>0`, and set `B_k=lcm b_ij`.  Since
`delta_ij=M a_ij/b_ij` is integral, `b_ij | M`, and coprimality gives the
stronger equality

```text
gcd(M,delta_ij)=M/b_ij.                              (1)
```

Taking the gcd of (1) over all pairs and using the elementary duality
between gcd and lcm among divisors of `M` gives

```text
B_k = M/gcd(M,h).                                    (2)
```

Thus the loss from the Gaussian norm `M` to the reciprocal denominator lcm
is precisely the part of the complementary-quotient lattice content `h`
that is supported on `M`.  Formula (2) is prime-by-prime: for every rational
prime `p`,

```text
v_p(B_k)=v_p(M)-min(v_p(M),v_p(h)).                   (3)
```

There is also an intrinsic finite-lattice version.  The Gaussian gcd of the
`q_j` is a unit.  Indeed, at each Gaussian prime the lcm `A` attains the
maximum valuation of one of the `w_j`, so one complementary quotient has
valuation zero.  Consequently the gcd of all real and imaginary coordinates
of all `q_j` is one.  Applying the maximal-minor index formula to the columns

```text
q_0,...,q_{k-1}, M e_1, M e_2
```

therefore yields

```text
gcd(M,h) = [Z^2 : Lambda + M Z^2].                   (4)
```

Equations (2) and (4) express `B_k` without choosing Gaussian associates or
factoring any `w_j`.

For the four-point endpoint target, (2) rewrites the missing estimate as

```text
M >= c^2 gcd(M,h) t_0^4/d^2.                         (5)
```

The exact identity is useful bookkeeping, but gives no bound on the new
factor `gcd(M,h)` by itself.

## 2. Changing the common denominator does not descend

For a fixed rational cotangent tuple `(t_j/d)`, its primitive integral
presentation is unique up to a positive common multiplier.  Replacing

```text
(d,t_j) by (c d,c t_j)
```

replaces `A` by an associate of `cA`, while every `q_j` and hence `h` stays
unchanged.  Therefore

```text
M_c=c^2 M,
B_k(c)=c^2 M/gcd(c^2 M,h).                            (6)
```

The target scale `t_0^4/d^2` also grows by `c^2`.  After division by that
scale, (6) can only stay fixed or decrease as `c` grows.  Removing the common
factor is consequently the strongest presentation for a prospective lower
bound; varying `d` through nonprimitive presentations cannot create a
descent.

Likewise, passing from `Lambda` to its saturation removes the determinant
factor `h`, but the required change of coordinates has determinant `1/h`
and is generally not a Euclidean similarity.  It does not preserve complex
norms, multiplication, or the circle radius.  A primitive basis of the
saturated `Z`-lattice therefore cannot replace `M` by `M/h` in the Gaussian
configuration.

## 3. An exact primitive family with unbounded content

The preceding warning is visible in a simple four-point family.  Let `d` be
odd with `gcd(d,6)=1`, and take

```text
t_j=d+4j,       0<=j<=3.
```

This is a primitive common-real-part presentation.  Write

```text
u_j=2j+i(d+2j),       N_j=N(u_j)=d^2+4dj+8j^2.
```

Then `w_j=(1+i)u_j`.  The `u_j` are pairwise coprime in `Z[i]`: a common
Gaussian divisor of `u_i,u_j` divides `2(j-i)(1+i)`; every `N_j` is odd, and
the only remaining possible rational prime is 3, which is excluded by
`3` not dividing `d`.  Hence one may choose

```text
A=(1+i) product_j u_j,
q_j=product_(ell != j) u_ell,
M=2 product_j N_j.                                   (7)
```

Directly from the reciprocal determinant formula,

```text
|delta_ij|=2d(j-i) product_(ell != i,j) N_ell.        (8)
```

After division by `2d`, the six values in (8) are

```text
N_2 N_3,  2N_1 N_3,  3N_1 N_2,
d^2 N_3,  2d^2 N_2,  d^2 N_1.                        (9)
```

Their gcd is one.  To see this, suppose a prime `p` divides all six.  If
`p` does not divide `d`, the last three expressions force `p` to divide
`N_1,N_2,N_3` (the case `p=2` is already excluded by the first odd
expression).  But

```text
N_2-N_1=4(d+6),       N_3-N_2=4(d+10),
```

so an odd `p` would divide 16, a contradiction.  If `p` divides `d`, then
`p` is neither 2 nor 3 and the first expression is congruent to
`32*72` modulo `p`, again a contradiction.  It follows that

```text
h=2d,       gcd(M,h)=2d,
B_4=M/(2d)=(product_j N_j)/d.                         (10)
```

Thus Gaussian primitivity of the `q_j` does not imply saturation of their
integer span: the index `[Z^2:Lambda]=2d` is unbounded even in primitive
presentations.  This is an exact counterfamily to a descent based only on
choosing a primitive basis.  It is not a counterexample to the desired
four-point height estimate: here `t_0=d` and `B_4` is much larger than its
quadratic target.

The finite audit in `check_reciprocal_minor_lattice_content.py` verifies
(2)--(4), the scaling law (6), and all formulas in the family (7)--(10).

## 4. Consequence for the route

The reciprocal-minor lcm contains exactly the Gaussian norm after quotienting
by the `M`-primary part of the complementary-quotient lattice content.  Any
new height inequality from this route must therefore use arithmetic that
bounds or compensates `gcd(M,h)`.  Gaussian gcd one, primitive common
denominator, and saturation of the underlying integer lattice do not supply
that control on their own.
